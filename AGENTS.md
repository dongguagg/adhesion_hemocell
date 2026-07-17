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
