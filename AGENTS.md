# HemoCell 红细胞间粘附作用改造方案

## 1. 当前状态与执行边界

- 本文件是后续源码改造的实施方案。当前阶段完成方案维护；进入实施阶段后只修改 `src_hemocell` 中的 HemoCell 源代码，不修改 `twoCellShear` 或其他具体算例的 `.cpp`、`config.xml`、`.pos` 文件。算例接入方式由用户后续指导。
- 改造目标是：在不同红细胞的膜节点之间加入“短程 LJ 排斥 + 中程 Morse 吸引”的中心力；同一红细胞内部节点不施加该作用力。
- 低层代码继续沿用现有 `cellId` 不同的筛选逻辑，不在粒子内核中硬编码名为 `RBC` 的 cell type。与旧排斥函数一致，新作用默认面向所有不同 `cellId` 的节点对；具体算例只注册 RBC 时即为 RBC-RBC 粘附。
- 细胞-边界排斥 `applyBoundaryRepulsionForce()` 不在本次范围内；细胞膜本构力、IBM 插值/铺展和粒子序列化格式也不改。
- 保留原有两参数 `setRepulsion(...)` 及旧排斥公式；新增五参数 `setAdhesion(...)`，代码链直接仿照旧排斥链实现，不引入复杂配置对象或额外框架。

## 2. 已确认的现有代码链

现有细胞间排斥参数和计算路径为：

`config.xml` → 算例中的 `Config::read<T>()` → `HemoCell::setRepulsion(...)` → `HemoCellFields::{repulsionConstant, repulsionCutoff}` → `HemoCell::iterate()`/`writeOutput()` → `HemoCellFields::applyRepulsionForce()` → `HemoCellParticleField::applyRepulsionForce()` → `#define inner_loop`。

新增粘附链仿照上述结构：

`HemoCell::setAdhesion(...)` → `HemoCellFields` 粘附参数 → `HemoCell::iterate()`/`writeOutput()` → `HemoCellFields::applyAdhesionForce()` → `HemoCellParticleField::applyAdhesionForce()` → 粘附节点对内核。

相关源文件和后续职责如下：

| 文件 | 后续改造职责 |
| --- | --- |
| `src_hemocell/hemocell.h` | 声明粘附模式的公共 setter，并说明输入单位 |
| `src_hemocell/core/hemoCell.cpp` | 完成 SI/微米到 LBM 单位转换，保存参数并启用粘附力调度 |
| `src_hemocell/core/hemoCellFields.h/.cpp` | 保存粘附参数，增加与旧排斥力平行的 processing functional 和调度入口 |
| `src_hemocell/core/hemoCellParticleField.h/.cpp` | 声明/实现粘附力计算，写入分段解析力，并将固定 stencil 改为随 cutoff 扩展的通用半邻域遍历 |

`src_hemocell/config/config.h/.cpp` 的 XML 访问是通用的，不需要为新标签修改解析器。本阶段也不修改任何实际 `config.xml`；后续算例约定从 `<Adhesion>` 段读取参数后调用 `setAdhesion(...)`。

## 3. 势能求导后的实际力

对节点 `i` 和 `j`，定义

```text
dv   = x_i - x_j
r    = |dv|
rhat = dv / r
sigma = r0 / 2^(1/6)
```

节点 `i` 所受力为 `F_i = -dU/dr * rhat`：

```text
F_i(r) = (24 epsilon/r) [2 (sigma/r)^12 - (sigma/r)^6] rhat,
                                                               0 < r < r0

         2 alpha D0 [exp(-2 alpha (r-r0))
                     - exp(-alpha (r-r0))] rhat,              r0 <= r < rc

         0,                                                    r >= rc
```

并对节点 `j` 同时施加 `F_j = -F_i`。

实现时必须保持以下物理含义：

- `0 < r < r0` 时，LJ 标量系数为正，力沿 `rhat`，表现为排斥。
- `r0 < r < rc` 时，Morse 标量系数为负，力沿 `-rhat`，表现为吸引。
- `r = r0` 时两侧力均为零，力连续；分支归入任意一侧均不改变结果，代码统一采用 Morse 分支。
- `r >= rc` 严格置零，忠实实现给定的硬截断。`rc` 左侧的 Morse 力通常非零，因此这里存在力跳变；第一阶段不得擅自改成 shifted-force 或平滑开关函数。
- 在 `r0` 处，两个势能分别为 `-epsilon` 和 `-D0`。只有 `epsilon == D0` 时势能本身连续，但力始终连续；因此不强制两个可调参数相等，只在文档和日志中提示这一点。
- 精确 `r = 0` 时方向无定义且 LJ 奇异，因此在除法前保留一个最小的零距离判断，避免 NaN/Inf；不增加复杂的异常和诊断框架。

建议在内核中复用中间量以减少昂贵运算：先计算 `r2`，再由 `sigma2/r2` 得到六次幂和十二次幂；Morse 分支只计算一次 `e = exp(-alpha*(r-r0))`，再用 `e*e` 表示二倍指数项。

## 4. 参数格式与单位约定

后续算例接入时采用独立顶层段，XML 标签区分大小写：

```xml
<Adhesion>
    <r0>...</r0>           <!-- equilibrium distance [um] -->
    <rc>...</rc>           <!-- cutoff distance [um] -->
    <epsilon>...</epsilon> <!-- LJ energy scale [J] -->
    <D0>...</D0>           <!-- Morse well depth [J] -->
    <alpha>...</alpha>     <!-- Morse width parameter [um^-1] -->
</Adhesion>
```

公共接口建议为：

```cpp
setAdhesion(r0Micrometer, rcMicrometer,
            epsilonJoule, D0Joule, alphaPerMicrometer)
```

在 `param::lbm_*_parameters(...)` 已完成、`param::dx` 和 `param::df` 有效之后，由 setter 一次性转换并存储 LBM 参数：

```text
r0_lbm      = r0_um * 1e-6 / param::dx
rc_lbm      = rc_um * 1e-6 / param::dx
sigma_lbm   = r0_lbm / 2^(1/6)
epsilon_lbm = epsilon_J / (param::df * param::dx)
D0_lbm      = D0_J / (param::df * param::dx)
alpha_lbm   = alpha_um^-1 * param::dx / 1e-6
```

其中 `param::df` 是一个 LBM 力单位对应的牛顿值，`param::df * param::dx` 是能量换算尺度。内核只使用上述 LBM 参数，不在每个节点对中重复做单位换算。

setter 只做与公式直接相关的基本检查：

```text
r0 > 0
rc > r0
epsilon > 0
D0 > 0
alpha > 0
param::dx > 0
param::df > 0
```

日志风格仿照 `setRepulsion(...)`，记录输入参数并说明启用 adhesion；不增加复杂的参数管理或错误恢复逻辑。

## 5. 与现有 API/数据结构的兼容策略

1. 保留 `HemoCell::setRepulsion(T constant, T cutoff)`、`repulsionEnabled`、`repulsionTimescale`、`applyRepulsionForce()` 和旧公式，确保旧算例接口及物理公式不变；动态 stencil 修复也用于旧排斥力，以修正大 cutoff 时的漏算。
2. 平行新增 `HemoCell::setAdhesion(...)`、`adhesionEnabled`、`adhesionTimescale`、`setAdhesionTimeScaleSeperation(...)` 和 `applyAdhesionForce()`；命名和调用方式直接对应旧排斥链。
3. 在 `HemoCellFields` 中新增简单的 LBM 标量字段：`adhesionR0`、`adhesionCutoff`、`adhesionSigma`、`adhesionEpsilon`、`adhesionD0`、`adhesionAlpha`。`sigma` 由 `r0` 派生，不从 XML 读取。
4. 增加 `HemoAdhesionForce` processing functional，结构直接复制并调整现有 `HemoRepulsionForce`，由其调用粒子场的 `applyAdhesionForce()`。
5. 粘附力仍写入现有 `force_repulsion` 槽，以免扩大到粒子序列化和输出格式重构。本任务约定同一运行只启用 `setRepulsion()` 或 `setAdhesion()` 之一，不设计两种模型叠加。
6. 不增加 interaction-mode 类、配置封装层、自动参数回退或其他与本次公式无关的保护性代码。

## 6. 粘附节点对内核实施细节

按以下顺序处理每个候选节点对：

1. 跳过同一对象。
2. 跳过相同 `cellId`，确保不把粘附力施加到同一红细胞内部。
3. 保证一个无序节点对只计算一次。同一个 particle-grid bin 内仅保留 `j > i`，不同 bin 使用半邻域偏移，避免双计数。
4. 计算 `r2`，先拦截零距离，再计算 `r`。
5. `r >= rc` 时尽早跳过。
6. `r < r0` 使用 LJ 排斥力，否则使用 Morse 吸引力。
7. 通过现有 `force_repulsion` 槽对两节点执行等量反向累加，保持总内力为零。

除为支持动态 cutoff 所必需的邻域枚举外，不改变以下行为：

- 不同 `cellId` 的判断；
- 周期边界、block envelope 和现有半邻域的总体遍历结构；
- 每次作用力计算前清零 `force_repulsion`；
- 与膜本构力相加后再通过 IBM 铺展到流体的路径。

## 7. 修复固定邻域 stencil

现有 `applyRepulsionForce()` 只遍历当前 particle-grid bin 及各方向相邻的一层 bin，每个 bin 尺度为 1 LU。新模型暂定 `r0 = 0.3 um`、`rc = 1.1 um`；当 `dx = 0.5 um` 时：

```text
r0_lbm = 0.6 LU
rc_lbm = 2.2 LU
neighborRange = ceil(rc_lbm) = 3
```

因此必须删除 `rc_lbm <= 1` 的限制思路，直接修复 HemoCell 的固定 stencil：

1. 根据当前作用力 cutoff 计算 `neighborRange = ceil(cutoff_lbm)`。
2. 对网格偏移 `[-neighborRange, neighborRange]` 做三重循环，并检查邻居网格索引是否在本地 particle-grid 范围内。
3. 只保留一个字典序半空间，例如 `ox > 0`，或 `ox == 0 && oy > 0`，或 `ox == 0 && oy == 0 && oz >= 0`；同 bin 内再使用 `j > i`。这里 `ox/oy/oz` 是网格索引偏移，避免与物理分辨率 `param::dx` 混淆。这样每个节点对只计算一次。
4. 仍由真实欧氏距离 `r < cutoff` 决定是否施力，立方体候选范围只负责不漏掉节点对。
5. 该 cutoff 驱动的邻域枚举同时供旧 `applyRepulsionForce()` 和新 `applyAdhesionForce()` 使用，从 HemoCell 层面修复固定一层 stencil，而不是只为某个算例打补丁。

现有常用 `particleEnvelope`（如 20、25、30 LU）大于暂定的搜索半径 3 LU；实现时沿用 HemoCell 现有 block/envelope 通信机制，不额外设计新的通信层。

## 8. 算例接入边界

本阶段不修改任何具体算例。后续由用户指导算例侧完成以下调用：从 `config.xml` 的 `<Adhesion>` 段读取 `r0`、`rc`、`epsilon`、`D0`、`alpha`，再传给 `hemocell.setAdhesion(...)`。HemoCell 内核和 setter 不直接访问某个算例的 XML。

## 9. 验证计划

### 9.1 解析公式测试

对固定参数计算若干代表距离，并与独立手算/高精度参考值比较：

- `0.8*r0`：力为正、沿 `rhat`；
- `r0`：力为零；
- `(r0+rc)/2`：力为负、沿 `-rhat`；
- `rc` 和 `1.1*rc`：力为零；
- `r0` 左右各取一个很小扰动：两侧趋近零；
- 任意有效节点对：`F_i + F_j` 在浮点容差内为零。

同时用有限差分 `-[U(r+h)-U(r-h)]/(2h)` 验证远离 `r0`、`rc` 的解析力；不要跨越分段点做中心差分。

### 9.2 参数、单位与邻域测试

- 验证微米到 LU、焦耳到 LBM 能量、`um^-1` 到 `LU^-1` 的转换结果。
- 检查基本正值关系和 `rc > r0`，并验证五个输入到六个内部值（含派生 `sigma`）的换算。
- 用 `cutoff = 0.8、1.0、1.4、2.2 LU` 检查动态 stencil，特别验证网格索引相差 2 和 3、但真实距离仍小于 cutoff 的节点对不会漏算。

### 9.3 集成与回归测试

- 单进程两节点/两细胞静态测试确认三个分段和同 bin 不重复计数。
- 两个不同 `cellId` 的节点受力相反；同一 `cellId` 的节点对无该作用力。
- 以 1 个和至少 2 个 MPI rank 运行相同小算例，比较力和早期轨迹，确认 block/envelope 边界没有漏算或重复计算。
- 编译默认 HemoCell 库，并用旧排斥入口做回归检查，确认 `setRepulsion()` 和边界排斥没有被新链破坏。本阶段不修改或提交具体算例代码。
- 监测所有粒子力、速度和位置的有限性。由于现有 `FORCE_LIMIT` 只限制膜本构 `sv.force`，不能假定它会替新 `force_repulsion` 分量兜底；参数和时间步应通过收敛/稳定性试验确定，第一阶段不添加会改变模型的任意力截断。

## 10. 完成判据

只有同时满足以下条件，后续代码改造才算完成：

- `setAdhesion(...)` 接收五个外部参数，HemoCell 源码中没有硬编码调参值；未来算例使用的 XML 段名约定为 `<Adhesion>`；
- 物理单位到 LBM 单位的转换集中在 setter，内核只使用 LBM 值；
- LJ 段只排斥、Morse 段只吸引、截断外严格为零；
- 不同细胞节点成对等量反向受力，同一细胞节点不受该细胞间作用；
- 同 bin 节点对没有双计数，零距离不会产生 NaN/Inf；
- 邻域 stencil 能按 cutoff 动态扩展，并完整覆盖暂定 `rc = 1.1 um`（在 `dx = 0.5 um` 时为 `2.2 LU`、搜索到偏移 3）；
- 旧两参数排斥算例仍可编译运行，边界排斥不受影响；
- 解析、单位、动态邻域、单进程和 MPI 验证均通过；
- 本阶段没有修改 `twoCellShear` 或其他具体算例文件；
- 实施提交不夹带与本任务无关的格式化或重构。

## 11. 推荐实施顺序

1. 先建立解析力和单位换算参考数据。
2. 增加与旧排斥链平行的 `setAdhesion(...)`、参数字段、timescale、processing functional 和粒子场入口。
3. 将固定一层 stencil 改成 cutoff 驱动的通用半邻域遍历，并让旧排斥力和新粘附力共用该遍历逻辑。
4. 在粘附节点对内核接入去重、最小零距离判断和三段力，不动边界排斥。
5. 完成解析、单位、动态邻域、单进程和 MPI 对比，不修改具体算例。
6. 对照本文件“完成判据”逐项复核；具体算例接入和物理参数扫描等待用户后续指导。

## 12. `twoCellShear` 算例目标与本轮边界

`twoCellShear` 当前由 `src_hemocell/examples/oneCellShear` 复制而来，后续要改造成“双红细胞先静态粘附预松弛、再施加单边 Couette 剪切”的独立算例。本节是算例实施方案；撰写本节时不修改 `twoCellShear.cpp`、`config.xml`、`RBC.pos`、`RBC.xml` 或构建文件。

目标状态如下：

- 盒子名义物理尺寸为 `16 um × 16 um × 20 um`，坐标约定为：`x` 是流向，`y` 是涡量方向，`z` 是速度梯度和上下壁法向。
- `x`、`y` 周期，`z` 非周期；下壁 `z=0` 始终静止，上壁沿 `+x` 方向运动。
- 两个相同 RBC 的盘面均平行于 `xy` 平面，沿 `z` 方向上下叠放，初始最近膜节点距离接近 `r0=0.3 um`。
- `RBC.pos` 第一行生成的下方 RBC（`cellId=0`）始终固定：所有膜节点位置保持初始值、速度为零，不随流场平移、转动或形变；第二行生成的上方 RBC（`cellId=1`）保持正常流固耦合，可以位移和形变。
- 粘附力从初始时刻起启用；预松弛阶段上下壁均静止，必须推进完整的细胞-流体耦合迭代，使上方自由 RBC 在粘附力和膜本构力共同作用下相对固定下方 RBC 达到近似稳定。
- 预松弛结束后再启动上壁；下壁不动。上壁目标速度由 `config.xml` 中的物理剪切率计算，不在源码中写死。
- 使用两个 MPI rank，规则分块固定为 `1 × 1 × 2`，即只沿 `z` 方向分块。
- 本算例只调用新增的 `setAdhesion(...)`，不同时调用旧的 `setRepulsion(...)`。

## 13. 计算域、离散尺寸与 MPI 分块

### 13.1 网格尺寸

建议在 `config.xml` 中显式保存三个名义尺寸（单位 `um`），由算例转换为格点数：

```xml
<domain>
    <Lx>16.0</Lx>
    <Ly>16.0</Ly>
    <Lz>20.0</Lz>
    ...
</domain>
```

转换使用：

```text
nx = round(Lx * 1e-6 / dx)
ny = round(Ly * 1e-6 / dx)
nz = round(Lz * 1e-6 / dx)
```

当 `dx=0.5 um` 时得到 `nx=32`、`ny=32`、`nz=40`。这沿用 HemoCell 现有算例按“格点数乘 `dx`”描述盒子尺寸的约定。需要注意：非周期 `z` 方向的壁面格点位于 `0` 和 `nz-1`，因此用于壁速计算的离散壁间距是

```text
Hwall = (nz-1) * dx = 19.5 um
```

而不是名义值 `20 um`。如果后续要求壁面几何距离严格等于 `20 um`，应单独改为 `nz=41`；第一版按当前 HemoCell 约定采用 `40` 个格点。

### 13.2 固定 `1 × 1 × 2` 分块

不能只依赖 `defaultMultiBlockPolicy3D()` 自动选择分解方向。算例应显式用

```text
createRegularDistribution3D(nx, ny, nz, 1, 1, 2)
```

建立 `SparseBlockStructure3D`，再用 `OneToOneThreadAttribution` 对两个 block 和两个 MPI rank 做一一分配。启动时检查 `global::mpi().getSize()==2`；不是两核时输出清楚的错误并退出，避免实际分块与方案不符。

继续使用 `fluidEnvelope=2` 和 `particleEnvelope=20`。对暂定 `rc=2.2 LU` 而言，动态邻域需要三层 particle-grid bin；现有 `particleEnvelope=20` 足够覆盖该搜索和 RBC 跨 block 的副本同步。本算例的 RBC 位于盒子下半部，因此两个 rank 的细胞计算负载不会完全均衡，这是指定放置方式带来的预期结果，不应误判为分块失败。

## 14. 单边 Couette 边界与两阶段时间推进

### 14.1 不复用现有对称 Couette 初始化

现有 `iniLatticeSquareCouette(...)` 会设置

```text
top    = -0.5 * gammaDot * H
bottom = +0.5 * gammaDot * H
```

上下壁都会运动，不符合本算例。第一版不修改这个共享 helper，而是在 `twoCellShear.cpp` 内直接设置：

- `x`、`y` 周期，`z` 非周期；
- top 和 bottom 均使用现有 velocity boundary；
- bottom 始终为 `(0,0,0)`；
- 预松弛阶段 top 也为 `(0,0,0)`；
- 剪切阶段 top 为 `(uTopLbm,0,0)`。

目标壁速为：

```text
gammaDotLbm = gammaDotPhysical * dt
uTopLbm     = gammaDotLbm * (nz-1)
uTopSI      = gammaDotPhysical * (nz-1) * dx
```

以当前 `gammaDot=111 s^-1`、`dt=0.5e-7 s`、`dx=0.5 um`、`nz=40` 为例：

```text
uTopLbm = 2.1645e-4
uTopSI  = 2.1645e-3 m/s
```

初始化流体为静止平衡态，不提前建立 Couette 速度剖面。这样预松弛阶段没有剪切；上壁启动后，速度剖面由流体动力学自然发展。

### 14.2 预松弛和上壁启动

在 `<sim>` 中增加：

```xml
<tRelax>...</tRelax> <!-- 完整粘附预松弛步数 -->
<tRamp>...</tRamp>   <!-- 上壁从 0 线性增至目标速度的步数；0 表示立即启动 -->
```

预松弛不能沿用当前只调用 `lattice->collideAndStream()` 的 fluid warmup，因为那不会更新 RBC。加载两个 RBC 并启用 adhesion 后，`tRelax` 内每一步都调用完整 `hemocell.iterate()`，但 top/bottom 速度均保持为零。

第一版保留一个可调线性 ramp，避免上壁速度瞬时跳变给刚形成的粘附双联体带来非物理冲击：

```text
factor = clamp((iter-tRelax)/tRamp, 0, 1)
uTop   = factor * uTopTarget
```

当 `tRamp=0` 时，`iter>=tRelax` 直接使用目标速度。`tmax` 继续表示包含预松弛阶段在内的总 HemoCell 迭代终点，因此实际恒定剪切段长度为 `tmax-tRelax-tRamp`。

checkpoint 重启时根据已恢复的 `hemocell.iter` 重新判断当前处于 relaxation、ramp 还是 constant-shear 阶段，并在下一次 `iterate()` 前恢复相应的上壁速度；不能在重启后无条件从零壁速重新开始。

## 15. 两个 RBC 的初始位置估算

`RBC.pos` 的位置单位为 `um`，后三列为角度制旋转角。继续采用 `oneCellShear` 已验证的 `90 0 0`，使两个 RBC 盘面平行于 `xy` 平面，并令两者角度完全相同。

对当前 `RBC_FROM_SPHERE` 形状，`RBC.xml` 中半径为 `R=3.91 um`。连续双凹解析轮廓的最大半厚度约为

```text
hmax ≈ 1.1467 um
```

若上下两膜的名义间隙取 `r0=0.3 um`，中心距估算为

```text
deltaZcenter ≈ 2*hmax + r0
             ≈ 2.5934 um
```

第一版取易读的 `2.6 um`。建议的 `RBC.pos` 基线为：

```text
2
8.0 8.0 5.0 90 0 0
8.0 8.0 7.6 90 0 0
```

由此得到：

- 两个 RBC 在 `x/y` 方向同心、盘面平行；
- 名义最近表面距离约 `2.6-2*1.1467=0.3066 um`，刚好接近 `r0` 且位于 `rc=1.1 um` 内；
- 下 RBC 初始最低膜面约为 `z=3.85 um`，不会紧贴静止下壁；
- 上 RBC 初始最高膜面约为 `z=8.75 um`，整个双联体位于盒子下半部；
- `x/y` 方向膜边缘约在 `4.09–11.91 um`，不会与自身周期像初始接触。

离散三角网格、mesh inflate 和节点分布会让真实最小节点距离与 `0.3066 um` 略有差异。初始时刻必须输出一次两细胞节点位置并检查真实最小跨细胞节点距离；若偏差较大，只调整第二行的 `z=7.6`，保持两个 RBC 的 `x/y` 和旋转完全一致。初始距离不能明显小于 `r0`，否则 LJ 段可能在第一步产生很大的排斥力。

`RBC.pos` 的行序在本算例中具有明确含义：第一行必须是固定的下方 RBC，第二行必须是自由的上方 RBC。HemoCell 按读取顺序分配 `cellId`，所以实现中固定 `cellId=0`，不再增加额外的固定对象配置层。

## 16. `config.xml` 参数方案

建议将当前精简配置补全为以下几组；数值仍全部由算例配置读取，不在 C++ 中写死：

```xml
<parameters>
    <warmup>0</warmup>
    <outputDirectory>output</outputDirectory>
    <checkpointDirectory>checkpoint</checkpointDirectory>
    <logDirectory>log</logDirectory>
    <logFile>twoCellShear.log</logFile>
</parameters>

<ibm>
    <radius>3.91e-6</radius>
    <stepMaterialEvery>1</stepMaterialEvery>
    <stepParticleEvery>1</stepParticleEvery>
</ibm>

<domain>
    <Lx>16.0</Lx>
    <Ly>16.0</Ly>
    <Lz>20.0</Lz>
    <shearrate>111.0</shearrate>
    <rhoP>1025</rhoP>
    <nuP>1.1e-6</nuP>
    <dx>0.5e-6</dx>
    <dt>0.5e-7</dt>
    <fluidEnvelope>2</fluidEnvelope>
    <particleEnvelope>20</particleEnvelope>
    <kBT>4.100531391e-21</kBT>
</domain>

<Adhesion>
    <r0>0.3</r0>
    <rc>1.1</rc>
    <epsilon>...</epsilon>
    <D0>...</D0>
    <alpha>...</alpha>
</Adhesion>

<sim>
    <tRelax>...</tRelax>
    <tRamp>...</tRamp>
    <tmax>...</tmax>
    <tmeas>...</tmeas>
    <tcheckpoint>...</tcheckpoint>
</sim>
```

`epsilon`、`D0`、`alpha` 和足够的 `tRelax` 尚无用户确认值，方案不替用户假定。两个细胞、短程强非线性力且计算量很小，初始基线将 material、particle velocity 和 adhesion 的 timescale separation 都设为 `1`；等稳定性和收敛性确认后再考虑加速，避免第一轮调参同时混入时间尺度分离误差。

`RBC.xml` 第一版保持现有高阶 RBC 本构参数和 `radius=3.91e-6 m` 不变；`config.xml` 中的 `<ibm><radius>` 必须与其保持一致。已有细胞-边界排斥第一版不自动启用，因为当前初始壁面间隙充足且用户尚未给出壁面排斥参数；运行中如出现膜节点接近或穿入壁面，再单独增加 `BRepCutoff/kBRep` 并调用现有 `enableBoundaryParticles(...)`，不能用 cell-cell Adhesion 参数代替壁面模型。

## 17. `twoCellShear.cpp` 实施顺序

后续实现按以下顺序进行：

1. 从配置读取 `Lx/Ly/Lz/dx/dt/shearrate`，计算 `32×32×40` 网格并调用 `param::lbm_shear_parameters(*cfg,nz)`；日志同时打印名义尺寸、实际壁间距和目标上壁速度。
2. 校验 MPI size 为 2，显式创建 `1×1×2` 的规则 block 管理，再构建 HemoCell lattice；不依赖默认自动分解。
3. 在算例内建立 `x/y` 周期、`z` 非周期的上下 velocity boundary，将两壁初始速度设为零，并将整个流体初始化为静止平衡态。
4. 初始化 cell field，注册 `RbcHighOrderModel("RBC", RBC_FROM_SPHERE)`，从 `<ibm>` 设置 material/particle timescale。
5. 增加一个最小的 `HemoCell::setCellFixed(plint cellId)` 接口并在算例中调用 `setCellFixed(0)`。底层只需保存一个固定 cell ID；在粒子速度插值中令该 cell 的所有节点速度为零，在 `advanceParticles()` 中跳过这些节点的位置更新。这样下 RBC 从子步层面就不会发生平移或形变，而不是在每个完整迭代结束后再把已经移动的网格强行复位。固定节点仍参与 LJ/Morse 节点对计算，其反作用力仍保存在 `force_repulsion` 槽并铺展给流体；不能把固定节点从 adhesion 内核中删掉。该接口在 `initializeCellfield()` 之后、加载 fresh particles 或 checkpoint 之前调用。
6. 从 `<Adhesion>` 读取五个参数并调用 `hemocell.setAdhesion(...)`；调用发生在 LBM 单位参数和 cell field 已初始化之后。adhesion timescale 第一版设为 1。
7. RBC 输出至少包含 `OUTPUT_POSITION`、`OUTPUT_TRIANGLES`、`OUTPUT_VELOCITY`、`OUTPUT_FORCE` 和 `OUTPUT_FORCE_REPULSION`。最后一项虽然沿用旧名称，但在 adhesion 模式下保存的正是复用 `force_repulsion` 槽的粘附/LJ 相互作用力。
8. fresh run 时先完成可选的静止无粒子 fluid warmup，再读取两个 RBC；在粒子加载后立即写一次初始输出。checkpoint run 则加载 checkpoint，不能重新读 `RBC.pos`。固定 `cellId=0` 的设置由算例每次启动时完成，不需要加入 checkpoint 粒子序列化格式。
9. 主循环每一步先根据 `iter/tRelax/tRamp` 设置本步 top 速度，再调用一次 `hemocell.iterate()`；按 `tmeas` 输出、按 `tcheckpoint` 保存 checkpoint。
10. 删除当前只针对 `cellId=0` 的单细胞 stretch 日志逻辑，改为遍历 `CellInformationFunctionals` 返回的两个 cell ID，记录阶段、迭代数、两个中心、包围盒、面积、体积、中心相对位移和当前上壁速度。日志中明确标记 `cellId=0` 为 fixed、`cellId=1` 为 free。
11. 日志统一改名为 `TwoCellShear`，结束时释放 boundary-condition 对象并正常关闭自定义日志文件。

预松弛“稳定”第一版由固定 `tRelax` 控制，而不是在代码中加入复杂的自动收敛判据。调参时根据两细胞中心距、最近节点距、面积/体积变化和粘附力随时间的平台区间选择足够的 `tRelax`。

除上述按 cell ID 固定粒子所需的最小 hook 外，算例实现优先直接参考和复用 `src_hemocell/examples` 中已有写法：`oneCellShear` 的剪切参数和输出流程、`cube` 的显式上下壁边界、`cellCollision_interior_viscosity` 的双细胞流程、`stretchCell`/`HemoCellStretch` 的轻量 processing-functional 组织方式，以及其他算例的 timescale、checkpoint 和 CMake 接入方式。不新建 interaction manager、自动收敛控制器、配置封装类或大量兜底分支；只保留参数直接读取、两核检查和现有 HemoCell 常规错误处理。

## 18. 构建与运行接入

`twoCellShear` 位于仓库根目录，而复制来的 `compile.sh` 假定算例位于 `src_hemocell/examples/...`，所以当前 `../../build` 路径不成立。实施时不把 `twoCellShear` 注册进 HemoCell 的 `examples/CMakeLists.txt`，也不把算例目录再复制一份；改为两个彼此清楚分离的构建阶段。

### 18.1 第一阶段：编译 HemoCell 静态库

HemoCell 的 configure 和所有库编译中间文件统一放在：

```text
/home/jxh/adhesion_rbc/src_hemocell/build
```

第一阶段等价于：

```bash
cmake -S /home/jxh/adhesion_rbc/src_hemocell \
      -B /home/jxh/adhesion_rbc/src_hemocell/build \
      -DBUILD_TESTING=OFF \
      -DCMAKE_C_COMPILER=/bin/mpicc \
      -DCMAKE_CXX_COMPILER=/bin/mpicxx \
      -DMPI_C_COMPILER=/bin/mpicc \
      -DMPI_CXX_COMPILER=/bin/mpicxx \
      -DMPI_CXX_SKIP_MPICXX=TRUE \
      -DCMAKE_CXX_FLAGS=-DOMPI_SKIP_MPICXX

cmake --build /home/jxh/adhesion_rbc/src_hemocell/build \
      --target hemocell --parallel
```

完成后必须存在：

```text
/home/jxh/adhesion_rbc/src_hemocell/build/libhemocell.a
```

这里采用正确目录名 `build`；用户描述中的 `buiild` 视为笔误。第一阶段只要求默认 `hemocell` 库，不额外构建 examples、tests、interior-viscosity 或 solidify 变体。

### 18.2 第二阶段：独立编译并链接 `twoCellShear`

`twoCellShear/CMakeLists.txt` 改成一个可独立 configure 的小型 CMake 工程，职责如下：

1. 定义可缓存的 `HEMOCELL_SOURCE_DIR`，默认值为 `${CMAKE_CURRENT_SOURCE_DIR}/../src_hemocell`。
2. 将 `${HEMOCELL_SOURCE_DIR}/build/libhemocell.a` 声明为 imported static library；configure 时如果该文件不存在，立即给出“请先运行 compile.sh 第一阶段”的错误。
3. 为算例添加 HemoCell、core、config、helper、io、mechanics、Palabos 和 bundled external-library 所需 include 目录。
4. 用 `find_package(MPI REQUIRED)` 和 `find_package(HDF5 REQUIRED COMPONENTS C HL)` 找到与第一阶段一致的 MPI/HDF5，并将 `twoCellShear` 链接到 `libhemocell.a`、MPI、HDF5 C/HL 及其必要系统链接项。
5. 为算例保持与 HemoCell 一致的 `PLB_MPI_PARALLEL`、`PLB_USE_POSIX`、`PLB_SMP_PARALLEL` 和 `OMPI_SKIP_MPICXX` 编译定义，避免头文件配置与静态库 ABI 不一致。
6. 用 `add_executable(twoCellShear twoCellShear.cpp)` 固定目标名，不依赖 `${PROJECT_NAME}` 推导目标。
7. 设置 `CMAKE_RUNTIME_OUTPUT_DIRECTORY=${CMAKE_CURRENT_SOURCE_DIR}`，因此最终可执行文件固定为：

```text
/home/jxh/adhesion_rbc/twoCellShear/twoCellShear
```

第二阶段自己的 CMake cache 和对象文件放在：

```text
/home/jxh/adhesion_rbc/twoCellShear/build
```

不能把这些中间文件直接散落在 `twoCellShear` 源码目录；只有最终可执行文件输出到该目录。

### 18.3 `compile.sh` 的两阶段流程

`compile.sh` 必须根据脚本自身位置计算绝对路径，不能依赖用户从哪个目录调用。脚本顺序固定为：

1. source `src_hemocell/loadHemoCell.sh`；
2. configure `src_hemocell/build`；
3. 只构建 `hemocell` target，并检查 `src_hemocell/build/libhemocell.a`；
4. configure `twoCellShear/build`，显式传入 `HEMOCELL_SOURCE_DIR`；
5. 构建 `twoCellShear` target，并检查 `twoCellShear/twoCellShear`；
6. 任意一步失败立即退出，不能在旧静态库或旧可执行文件存在时误报成功。

第一次执行会完整编译 HemoCell；以后 HemoCell 源码未变化时，`cmake --build` 使用增量编译，只重新构建发生变化的对象。修改 `twoCellShear.cpp` 时，第二阶段只需重新编译算例并重新链接已有 `libhemocell.a`。

期望运行方式为：

```bash
source /home/jxh/adhesion_rbc/src_hemocell/loadHemoCell.sh
cd /home/jxh/adhesion_rbc/twoCellShear
./compile.sh
mpirun -np 2 ./twoCellShear config.xml
```

`src_hemocell/build`、`twoCellShear/build`、最终可执行文件和 `output[_N]/tmp[_N]` 结果均由仓库根目录 `.gitignore` 排除。

## 19. 算例验证计划

### 19.1 几何与初始化

- 启动日志确认域为 `32×32×40 LU`、block 为 `1×1×2`、MPI rank 数为 2。
- 初始输出确认恰好存在两个不同 `cellId`，中心接近 `(8,8,5.0)` 和 `(8,8,7.6) um`，盘面平行且未被边界加载逻辑删点。
- 计算初始最小跨细胞节点距离，目标约为 `0.3 um`，并确认没有节点距离小于 `r0` 很多。
- 确认 `r0=0.6 LU`、`rc=2.2 LU`、动态 stencil range 为 3，`particleEnvelope=20` 足够。

### 19.2 静止粘附预松弛

- relaxation 全程 top/bottom 速度均为零，fluid 无外加体力。
- 下 RBC 的所有节点位置、中心、面积、体积和包围盒必须保持初始值，节点输出速度为零；上 RBC 可以在粘附力下移动和形变。
- 最近节点距离向 `r0` 附近调整，跨细胞作用力仍按节点对等量反向；固定下 RBC 所受反力作为锚固反力保留并铺展到流体，不能要求双联体质心守恒漂移为零。
- 上 RBC 的面积、体积、包围盒、中心距和相对速度最终进入平台区；据此选择 `tRelax`，不能只看单个节点距离。

### 19.3 上壁启动与剪切阶段

- bottom 始终为零，top 按 ramp 到达 `uTopTarget`，日志中的 LBM/SI 壁速与配置剪切率一致。
- 无细胞基准或远离细胞区域的速度剖面最终接近 `u_x(z)=gammaDot*z`。
- 观察双联体是否保持粘附、滑移、滚动或解离，并同时检查 `OUTPUT_FORCE_REPULSION`、两细胞相对位移及形变。
- 剪切全过程再次确认下 RBC 无位移、无转动、无形变，上 RBC 的响应不受固定逻辑误伤。
- checkpoint 分别在 relaxation、ramp 和恒定剪切阶段重启一次，确认壁速阶段不会重置。

### 19.4 稳定性和参数扫描

- 首先以较短 `tmax`、高频 `tmeas` 做试运行，再增加总时长。
- 对 `epsilon/D0/alpha` 和 `dt` 做独立扫描；不要同时改变网格、timescale separation 和粘附参数。
- 检查 `uTopLbm` 保持低 Mach 数；当前示例值约 `2.16e-4`，远低于常用稳定性上限。
- 检查 RBC 是否接近上下壁；若需要壁面排斥，单独建立无 adhesion 的墙面接触基准后再与粘附模型组合。

## 20. `twoCellShear` 完成判据

只有同时满足以下条件，算例改造才算完成：

- `RBC.pos` 含两个平行 RBC，初始真实最近节点距接近 `0.3 um`，且没有明显初始重叠；
- 下方 `cellId=0` 的所有节点在 relaxation、ramp、恒定剪切和 checkpoint 重启后都保持初始位置与零速度，不发生位移、转动或形变；上方 `cellId=1` 正常位移和形变；
- 域、周期性、单边壁速和 `1×1×2` MPI 分块与本方案一致；
- 五个 adhesion 参数全部由 `<Adhesion>` 读取并传给 `setAdhesion(...)`，无算例内硬编码物理调参值；
- relaxation 推进完整 HemoCell 耦合过程，两壁静止，粘附双联体达到可辨识的平台状态；
- 剪切阶段只有上壁运动，壁速严格由剪切率、`dt` 和离散壁间距计算；
- fresh run 和三个阶段的 checkpoint restart 均行为一致；
- 两个 cell ID 都进入日志与输出，adhesion 槽力可供后处理；
- 两核运行无跨 block 漏力、重复力、NaN/Inf 或粒子被意外删除；
- `compile.sh` 能按“`src_hemocell/build/libhemocell.a` → 独立链接 twoCellShear”的两阶段流程工作，最终可执行文件位于 `twoCellShear/twoCellShear`，构建和结果文件均不进入 Git。
