# HemoCell 红细胞间及红细胞-壁面粘附作用改造方案

## 1. 当前状态、策略变更与执行边界

- 当前源码已经实现不同 `cellId` 膜节点之间的“短程 LJ 排斥 + 中程 Morse 吸引”，并在 `twoCellShear` 中通过 `HemoCell::setCellFixed(0)`、粒子速度清零和跳过位置推进来固定下方 RBC。
- 新策略废弃上述按 `cellId` 强制固定方式。下方 RBC 改回正常 RBC：参与流体速度插值、粒子推进、膜本构、细胞间粘附、IBM 力铺展，可以平移、转动和变形。
- 为使下方 RBC 被物理地锚定并继续阻碍流体，在膜节点与固体壁面格点之间新增与细胞间粘附同形式、但参数独立的“LJ 排斥 + Morse 吸引”。下方 RBC 初始放在下壁附近，由细胞-壁面粘附限制其运动，而不是由运动学 hook 锁死。
- 本轮只维护本方案并写入 `AGENTS.md`，不修改任何 HemoCell 源码、`twoCellShear.cpp`、`config.xml` 或 `RBC.pos`。经用户确认后再进入实施阶段。
- 实施阶段允许修改完成本策略所需的 `src_hemocell` 源码和 `twoCellShear` 算例文件；不修改其他算例，也不夹带无关格式化、框架化重构或物理模型变更。
- 保留旧 `setRepulsion(...)`、`enableBoundaryParticles(...)`、`applyBoundaryRepulsionForce()` 及其公式，保证旧接口可用；新增壁面粘附链，不用新模型替换或删除旧边界排斥。
- 细胞膜本构、IBM 核函数、流体碰撞推进、粒子序列化和现有输出格式不改。两类粘附力继续复用 `force_repulsion` 槽，并与膜本构力相加后铺展到流体。

## 2. 已确认的代码链与后续职责

现有细胞间粘附链为：

`config.xml` → `twoCellShear` 中的 `Config::read<T>()` → `HemoCell::setAdhesion(...)` → `HemoCellFields::{adhesionR0,...}` → `HemoCell::iterate()`/`writeOutput()` → `HemoCellFields::applyAdhesionForce()` → `HemoCellParticleField::applyAdhesionForce()` → 不同 `cellId` 节点对内核。

现有边界排斥链为：

`HemoCell::enableBoundaryParticles(...)` → `HemoCellFields::populateBoundaryParticles()` → `HemoCellParticleField::populateBoundaryParticles()` → `HemoCell::iterate()`/`writeOutput()` → `HemoCellFields::applyBoundaryRepulsionForce()` → `HemoCellParticleField::applyBoundaryRepulsionForce()`。

新增细胞-壁面粘附链仿照边界排斥链实现：

`<cellWallAdhesion>` → `twoCellShear` 读取五个参数 → `HemoCell::setBoundaryAdhesion(...)` → LBM 参数换算并复用 `populateBoundaryParticles()` → `HemoCell::iterate()`/`writeOutput()` → `HemoCellFields::applyBoundaryAdhesionForce()` → `HemoCellParticleField::applyBoundaryAdhesionForce()` → 膜节点与壁面格点的 LJ/Morse 内核。

相关文件职责如下：

| 文件 | 后续改造职责 |
| --- | --- |
| `src_hemocell/hemocell.h` | 删除 `setCellFixed(...)`；声明 `setBoundaryAdhesion(...)` 和启用标志，说明五个输入的单位 |
| `src_hemocell/core/hemoCell.cpp` | 删除固定 cell setter；换算并保存壁面粘附参数；在 `iterate()`/`writeOutput()` 中按正确的清零和累加顺序调度两类粘附力 |
| `src_hemocell/core/hemoCellFields.h/.cpp` | 删除固定 cell 状态和判断；保存壁面粘附 LBM 标量；新增清零入口、`HemoBoundaryAdhesionForce` functional 和 `applyBoundaryAdhesionForce()` |
| `src_hemocell/core/hemoCellParticleField.h/.cpp` | 恢复所有 RBC 的正常插值和推进；实现壁面 LJ/Morse 力及 cutoff 驱动的壁面邻域搜索 |
| `twoCellShear/twoCellShear.cpp` | 删除固定 cell 调用及 fixed/free 语义；读取两个独立粘附段并启用细胞间和细胞-壁面粘附 |
| `twoCellShear/config.xml` | 将 `<Adhesion>` 改为 `<cellCellAdhesion>`，新增 `<cellWallAdhesion>`；两段各保存五个独立参数 |
| `twoCellShear/RBC.pos` | 将两个 RBC 整体下移，使下 RBC 靠近下壁，同时保持两 RBC 的细胞间初始间隙 |

`src_hemocell/config/config.h/.cpp` 的 XML 读取能力是通用的，不需要针对新标签修改解析器。

## 3. 两类粘附共用的解析力

对受力膜节点 `i` 和相互作用对象 `j`，定义

```text
dv    = x_i - x_j
r     = |dv|
rhat  = dv / r
sigma = r0 / 2^(1/6)
```

膜节点 `i` 所受力为：

```text
F_i(r) = (24 epsilon/r) [2 (sigma/r)^12 - (sigma/r)^6] rhat,
                                                               0 < r < r0

         2 alpha D0 [exp(-2 alpha (r-r0))
                     - exp(-alpha (r-r0))] rhat,              r0 <= r < rc

         0,                                                    r >= rc
```

- 对细胞间粘附，`j` 是另一细胞膜节点，并同时累加 `F_j=-F_i`。
- 对细胞-壁面粘附，`j` 是 `populateBoundaryParticles()` 生成的静态固体边界格点；只向膜节点累加 `F_i`，不创建可移动的壁面粒子，也不向壁面格点保存一个粒子力。固定速度壁面的反力由流体边界条件承担。
- 壁面作用中 `rhat` 从壁面格点指向膜节点。因此 `0<r<r0` 时 LJ 力把膜节点推离壁面，`r0<r<rc` 时 Morse 力把膜节点拉向壁面。
- `r=r0` 时两侧力均为零，统一进入 Morse 分支；`r>=rc` 严格为零，保留硬截断及其在 `rc` 处的力跳变，不擅自增加 shifted-force 或平滑开关。
- 两套参数分别允许 `epsilon!=D0`。此时势能在 `r0` 处不连续，但力仍连续；不强制两者相等，只在文档和日志中提示。
- 精确 `r=0` 时方向无定义且 LJ 奇异，必须在除法和开方后的归一化之前拦截，避免 NaN/Inf；不添加任意力截断。
- 内核复用 `r2`、`sigma2/r2`、六次幂和十二次幂；Morse 分支只计算一次 `e=exp(-alpha*(r-r0))`，用 `e*e` 表示二倍指数项。

## 4. XML 标签、公共接口与单位

XML 标签区分大小写，实施后固定使用以下两个独立顶层段：

```xml
<cellCellAdhesion>
    <r0>...</r0>           <!-- equilibrium distance [um] -->
    <rc>...</rc>           <!-- cutoff distance [um] -->
    <epsilon>...</epsilon> <!-- LJ energy scale [J] -->
    <D0>...</D0>           <!-- Morse well depth [J] -->
    <alpha>...</alpha>     <!-- Morse width parameter [um^-1] -->
</cellCellAdhesion>

<cellWallAdhesion>
    <r0>...</r0>           <!-- equilibrium distance [um] -->
    <rc>...</rc>           <!-- cutoff distance [um] -->
    <epsilon>...</epsilon> <!-- LJ energy scale [J] -->
    <D0>...</D0>           <!-- Morse well depth [J] -->
    <alpha>...</alpha>     <!-- Morse width parameter [um^-1] -->
</cellWallAdhesion>
```

细胞间接口继续保留现有名称以避免不必要的 API 破坏；壁面接口新增为：

```cpp
setAdhesion(cellCellR0Micrometer, cellCellRcMicrometer,
            cellCellEpsilonJoule, cellCellD0Joule,
            cellCellAlphaPerMicrometer)

setBoundaryAdhesion(cellWallR0Micrometer, cellWallRcMicrometer,
                    cellWallEpsilonJoule, cellWallD0Joule,
                    cellWallAlphaPerMicrometer)
```

两者都在 `param::lbm_*_parameters(...)` 完成、`param::dx` 和 `param::df` 有效之后由 setter 一次换算：

```text
r0_lbm      = r0_um * 1e-6 / param::dx
rc_lbm      = rc_um * 1e-6 / param::dx
sigma_lbm   = r0_lbm / 2^(1/6)
epsilon_lbm = epsilon_J / (param::df * param::dx)
D0_lbm      = D0_J / (param::df * param::dx)
alpha_lbm   = alpha_um^-1 * param::dx / 1e-6
```

其中 `param::df` 是一个 LBM 力单位对应的牛顿值，`param::df*param::dx` 是能量换算尺度。两个粒子内核只使用 LBM 参数，不在相互作用循环中做 SI 换算。

两个 setter 分别检查 `r0>0`、`rc>r0`、`epsilon>0`、`D0>0`、`alpha>0`、`param::dx>0` 和 `param::df>0`，并分别打印 cell-cell/cell-wall 输入值、换算值和启用日志。五个 cell-wall 参数完全独立于 cell-cell 参数，不复制、不回退、不共享默认值。

## 5. API、字段和调度兼容策略

1. 保留现有 `setAdhesion(...)`、`adhesionEnabled`、`adhesionTimescale` 和 `applyAdhesionForce()`，它们明确表示 cell-cell adhesion；只把算例 XML 来源从 `<Adhesion>` 改成 `<cellCellAdhesion>`。
2. 新增 `boundaryAdhesionEnabled` 和六个 LBM 字段：`boundaryAdhesionR0`、`boundaryAdhesionCutoff`、`boundaryAdhesionSigma`、`boundaryAdhesionEpsilon`、`boundaryAdhesionD0`、`boundaryAdhesionAlpha`。
3. `setBoundaryAdhesion(...)` 仿照 `enableBoundaryParticles(...)` 调用 `populateBoundaryParticles()`，但保存的是五参数 LJ/Morse 模型并启用 `boundaryAdhesionEnabled`。调用发生在上下壁 dynamics 已建立且 `initializeCellfield()` 已完成之后。
4. 第一版让 cell-cell 和 cell-wall adhesion 共用现有 `adhesionTimescale`，`twoCellShear` 固定将其设为 `1`。不新增两个可独立错开的粘附刷新频率，避免共享 `force_repulsion` 槽时一个模型刷新而另一个模型留下陈旧分量。
5. 在一次 adhesion 更新中，先用一个轻量 processing functional 将所有粒子的 `force_repulsion` 清零一次，再调用 `applyAdhesionForce()` 累加 cell-cell 力，最后调用 `applyBoundaryAdhesionForce()` 累加 cell-wall 力。`writeOutput()` 必须使用完全相同的清零和重算顺序。
6. 当前 `applyAdhesionForce()` 内部的清零逻辑移到上述统一入口；`applyBoundaryAdhesionForce()` 只做累加，不能再次清零，否则会覆盖 cell-cell adhesion。只有壁面粘附启用时也由统一入口先清零，所以不会逐步累积旧力。
7. 两类粘附仍共同写入 `sv.force_repulsion`，`OUTPUT_FORCE_REPULSION` 因而输出两者之和。第一版不增加新粒子字段、序列化字段或单独输出通道；后处理若要区分两类贡献，需另立后续任务。
8. 本算例不同时启用旧 `setRepulsion(...)` 或旧 `enableBoundaryParticles(...)`。旧 cell-cell repulsion 与旧 boundary repulsion 的 API、公式和单独运行行为保持不变。
9. 不增加 interaction manager、配置封装类、自动参数回退或通用势函数框架。

## 6. 细胞间粘附节点对内核

现有内核的物理和去重规则继续保留：

1. 跳过同一对象和相同 `cellId`，不对同一 RBC 内部节点施加该力。
2. 同一个 particle-grid bin 仅处理 `j>i`；不同 bin 使用字典序半邻域，只计算一次无序节点对。
3. 先计算 `r2`，拦截零距离和 `r>=rc`，再进入 LJ 或 Morse 分支。
4. 对两个膜节点等量反向累加到 `force_repulsion`。
5. 邻域范围继续由 `neighborRange=ceil(cellCellRcLbm)` 动态确定，真实欧氏距离决定是否施力。

本轮对该内核唯一必要的结构调整，是移除其私有清零步骤，让统一粘附调度先清零、cell-cell 与 cell-wall 两个内核随后安全叠加。

## 7. 细胞-壁面粘附内核与邻域搜索

`applyBoundaryAdhesionForce()` 直接参考 `applyBoundaryRepulsionForce()` 的壁面表示和 block/envelope 机制，但使用以下规则：

1. 若 particle grid 未更新，先调用现有 `update_pg()`。
2. 使用 `populateBoundaryParticles()` 已生成的静态固体-流体界面格点。第一版对所有被该函数识别的固体壁面使用同一套 `<cellWallAdhesion>` 参数；`twoCellShear` 通过几何位置使下壁是实际发生粘附的壁面，不增加只选择 `z=0` 的专用内核开关。
3. 根据 `boundaryAdhesionCutoff` 计算 `neighborRange=ceil(cutoffLbm)`，对每个壁面格点搜索三轴偏移 `[-neighborRange,neighborRange]` 内合法的 particle-grid bin，不能沿用旧边界排斥固定的 `±1` 搜索层。
4. 这是“膜节点—壁面格点”的二部配对，不使用 cell-cell 半邻域。每个壁面格点只遍历一次、每个候选 bin 中的膜节点只遍历一次，同一膜节点可以与 cutoff 内多个不同壁面格点发生作用；这是现有离散壁面粒子模型的含义，不改成最近点到连续平面的距离。
5. 使用全局坐标计算 `dv=membranePosition-(boundaryPoint+atomicLatticeLocation)` 和真实 `r2`。先拦截 `r2<=0` 与 `r2>=rc^2`，再按共享公式计算 LJ/Morse 力。
6. 不做 `cellId` 筛选：任何 cell type 的膜节点接近固体边界时都会受该壁面力，与现有 boundary repulsion 的低层语义一致。`twoCellShear` 只注册 RBC，因此实际是 RBC-wall adhesion。
7. 只向膜节点的 `force_repulsion` 加力；壁面不推进。该力与膜本构力一起经 `spreadParticleForce()` 铺展到流体，使正常演化的下方 RBC 对流体产生阻碍和反馈。
8. `particleEnvelope` 必须大于两套 cutoff 所需的搜索范围；两 MPI rank 下重点检查下壁所在 block 不漏算、不重复算。

## 8. 撤销固定 cellId 策略

实施时完整移除本项目为固定下 RBC 增加的最小 hook：

- 从 `hemocell.h` 删除 `setCellFixed(plint cellId)` 声明，从 `hemoCell.cpp` 删除其实现。
- 从 `HemoCellFields` 删除 `fixedCellEnabled`、`fixedCellId` 和 `isCellFixed(...)`。
- 从 `HemoCellParticleField::interpolateFluidVelocity()` 删除对固定 cell 的速度清零和 `continue`，恢复所有粒子从流体插值速度。
- 从 `HemoCellParticleField::advanceParticles()` 删除对固定 cell 跳过 `particle.advance()` 的分支，恢复正常位置更新和边界穿入检查。
- 从 `twoCellShear.cpp` 删除 `hemocell.setCellFixed(0)`、相关注释、fixed/free 日志标签和 checkpoint 重启时的固定语义。
- 不通过每步复位坐标、速度钳制、额外弹簧或其他隐式固定方法替代。下 RBC 是否保持在下壁附近只由 cell-wall adhesion、膜力、流体力和初始几何决定。

## 9. 源码级验证计划

### 9.1 解析式、单位和符号

- 分别用 cell-cell 和 cell-wall 两套独立参数检查 `0.8*r0` 为排斥、`r0` 为零、`(r0+rc)/2` 为吸引、`rc` 与 `1.1*rc` 为零。
- 在不跨越 `r0`/`rc` 的距离上，用势能中心差分验证解析力。
- 验证微米到 LU、焦耳到 LBM 能量、`um^-1` 到 `LU^-1` 的换算；检查五个输入得到六个内部量（含派生 `sigma`）。
- cell-cell 节点对检查 `F_i+F_j=0`；cell-wall 检查 LJ 指向离壁、Morse 指向壁面。
- 零距离必须跳过，所有输出力、速度和位置保持有限。

### 9.2 邻域、清零与累加

- 用 `cutoff=0.8、1.0、1.4、2.2 LU` 验证 cell-cell 与 cell-wall 搜索，特别覆盖索引差 2、3 且真实距离小于 cutoff 的组合。
- 单独启用 cell-cell、单独启用 cell-wall、同时启用两者三种模式，确认 `force_repulsion` 不跨时间步累积，也不会因第二个内核清零而丢失第一个内核贡献。
- 同一次输出重算与正常迭代得到相同 interaction force，`writeOutput()` 不改变下一步的力槽状态。
- 对平面壁的多个壁面格点贡献建立小测试，确认同一“膜节点—壁面格点”配对没有重复，而不同壁面格点的贡献按模型正常相加。

### 9.3 MPI 与回归

- 以 1 个和至少 2 个 MPI rank 对比靠近 block/envelope 边界的膜节点壁面力，确认没有漏算、重复计算或 rank 相关差异。
- 编译默认 `hemocell` 静态库，回归旧 `setRepulsion()`、旧 `enableBoundaryParticles()` 和边界排斥公式。
- 验证撤销固定 hook 后两个 cellId 都执行插值和推进，checkpoint 不依赖任何未序列化的 fixed-cell 状态。
- `FORCE_LIMIT` 只限制膜本构 `sv.force`，不能假定它替 `force_repulsion` 中两类粘附力兜底；稳定性通过参数和 `dt` 扫描确定，不添加任意钳制。

## 10. 源码改造完成判据

- `<cellCellAdhesion>` 与 `<cellWallAdhesion>` 各自提供五个参数，源码没有硬编码调参值或跨段回退。
- `setAdhesion(...)` 和 `setBoundaryAdhesion(...)` 在 setter 中完成所有单位换算，内核只使用各自的 LBM 参数。
- `applyBoundaryAdhesionForce()` 复用离散边界格点，按 wall cutoff 动态扩展搜索，严格实现 LJ/Morse/硬截断和零距离保护。
- cell-cell 与 cell-wall 力在每次更新前只清零一次并正确相加；`OUTPUT_FORCE_REPULSION` 是两类 interaction force 的总和。
- 下 RBC 不再有任何按 cellId 的速度清零、跳过推进或坐标复位逻辑，两个 RBC 都是正常可变形、可移动的 IBM 细胞。
- 壁面粘附力随膜节点总力铺展到流体，正常 RBC 对流体的阻碍作用得到保留。
- 旧 cell-cell repulsion、旧 boundary repulsion、膜本构、IBM、序列化和输出格式不被破坏。
- 解析式、单位、动态邻域、力槽累加、单进程、MPI 和回归验证通过，无 NaN/Inf 或粒子意外删除。

## 11. 推荐实施顺序

1. 先建立两套参数的解析力与单位换算参考数据。
2. 删除 `setCellFixed`、固定状态字段以及插值/推进中的固定分支，先恢复两个 RBC 的标准 HemoCell 行为。
3. 增加 `setBoundaryAdhesion(...)`、LBM 参数字段、启用标志、processing functional 和粒子场入口。
4. 将 adhesion 力槽清零移到统一调度入口，按“清零一次 → cell-cell → cell-wall”接入 `iterate()` 和 `writeOutput()`。
5. 实现 cutoff 驱动的壁面邻域遍历和 LJ/Morse 内核，保留旧 `applyBoundaryRepulsionForce()`。
6. 更新 `twoCellShear` 的 XML 标签、setter 调用、初始位置和日志语义。
7. 完成源码级测试、两核算例测试、checkpoint 与旧 API 回归，再逐项核对本文件完成判据。

## 12. `twoCellShear` 算例目标与本轮边界

`twoCellShear` 保持“双红细胞先静态粘附预松弛、再施加单边 Couette 剪切”的总体流程，但锚定机制从“强制固定 cellId=0”改为“下方 RBC 与静止下壁粘附”。本节是待用户确认的实施方案；撰写本节时不修改 `twoCellShear.cpp`、`config.xml`、`RBC.pos`、`RBC.xml` 或构建文件。

目标状态如下：

- 盒子名义物理尺寸为 `16 um × 16 um × 20 um`，坐标约定为：`x` 是流向，`y` 是涡量方向，`z` 是速度梯度和上下壁法向。
- `x`、`y` 周期，`z` 非周期；下壁 `z=0` 始终静止，上壁沿 `+x` 方向运动。
- 两个相同 RBC 的盘面均平行于 `xy` 平面并沿 `z` 方向上下叠放。初始跨细胞膜节点距离接近 `cellCellAdhesion/r0`，下 RBC 到离散下壁格点的最近距离接近 `cellWallAdhesion/r0`。
- `RBC.pos` 第一行仍生成下方 RBC（`cellId=0`），第二行生成上方 RBC（`cellId=1`），但行序只用于日志和几何识别，不再触发任何底层固定逻辑。
- 两个 RBC 都完整参与流固耦合，均可平移、转动和变形。下 RBC 预期由下壁粘附限制在壁面附近，但不要求所有节点、中心、面积或体积严格保持初值。
- cell-cell 和 cell-wall adhesion 从初始时刻同时启用；预松弛阶段上下壁均静止，必须推进完整的细胞-流体耦合，使下 RBC 与下壁、两个 RBC 之间的粘附及膜形变共同达到近似稳定。
- 预松弛结束后再启动上壁；下壁不动。上壁目标速度由 `config.xml` 中的物理剪切率计算，不在源码中写死。
- 使用两个 MPI rank，规则分块固定为 `1 × 1 × 2`，即只沿 `z` 方向分块。
- 本算例调用 `setAdhesion(...)` 和新增的 `setBoundaryAdhesion(...)`，不调用旧 `setRepulsion(...)` 或 `enableBoundaryParticles(...)`。

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

预松弛不能沿用只调用 `lattice->collideAndStream()` 的 fluid warmup，因为那不会更新 RBC。加载两个 RBC 并同时启用 cell-cell/cell-wall adhesion 后，`tRelax` 内每一步都调用完整 `hemocell.iterate()`，但 top/bottom 速度均保持为零。

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

记 cell-cell 与 cell-wall 的平衡距离分别为 `r0_cc` 和 `r0_cw`。忽略离散网格横向偏移时，中心高度的第一轮估算为：

```text
zLower       ≈ hmax + r0_cw
deltaZcenter ≈ 2*hmax + r0_cc
zUpper       ≈ zLower + deltaZcenter
```

当 `r0_cc=0.3 um` 时，两个中心仍相差约 `2.5934 um`，可取易读的 `2.6 um`。只有在用户确认 `cellWallAdhesion/r0` 后才能确定最终 `zLower`；例如若 `r0_cw=0.3 um`，连续几何给出 `zLower≈1.4467 um`、`zUpper≈4.0401 um`，可先试用：

```text
2
8.0 8.0 1.45 90 0 0
8.0 8.0 4.05 90 0 0
```

这个数值只是在 `r0_cw=0.3 um` 条件下的待验证基线，不作为未确认壁面参数的硬编码结论。几何调整遵循：

- 两个 RBC 在 `x/y` 方向同心、盘面平行；
- 上下两 RBC 的中心距保持约 `2*hmax+r0_cc`，因此移动下 RBC 时上 RBC 同步下移，不能破坏已有 cell-cell 初始间隙；
- 下 RBC 的名义最低膜面距 `z=0` 约为 `r0_cw`，确保一开始位于 wall cutoff 内但不深落入 LJ 强排斥段；
- 上 RBC 初始最高面仍远离上壁，避免上壁 cell-wall adhesion 在初始阶段参与；
- `x/y` 方向膜边缘约在 `4.09–11.91 um`，不会与自身周期像初始接触。

壁面内核使用膜节点到离散壁面格点的三维欧氏距离，而不是只用节点的 `z` 坐标；离散三角网格、mesh inflate、节点的 `x/y` 相位以及壁面格点间距都会使真实最近 wall-pair 距离偏离上述连续估算。因此 fresh run 在第一次迭代前必须计算并记录：

- 最小跨细胞节点距离，与 `r0_cc`、`rc_cc` 比较；
- 下 RBC 膜节点到下壁离散格点的最小距离，与 `r0_cw`、`rc_cw` 比较；
- 分别处于 cell-cell 和 cell-wall 的 LJ/Morse/cutoff 三个区间内的候选对数量。

若 wall 距离不合适，先整体同步调整两个 RBC 的 `z`；若 cell-cell 距离不合适，再微调第二行相对第一行的 `z`。两个初始最小距离都不能明显小于各自的 `r0`，否则第一步可能产生很大的 LJ 排斥并把膜节点推入边界删除逻辑。

`RBC.pos` 第一行仍必须是靠壁的下 RBC，第二行是上 RBC，以保持 `cellId=0/1` 的日志语义和后处理约定；该行序不再改变任何粒子的动力学处理。

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

<cellCellAdhesion>
    <r0>0.3</r0>
    <rc>1.1</rc>
    <epsilon>...</epsilon>
    <D0>...</D0>
    <alpha>...</alpha>
</cellCellAdhesion>

<cellWallAdhesion>
    <r0>...</r0>
    <rc>...</rc>
    <epsilon>...</epsilon>
    <D0>...</D0>
    <alpha>...</alpha>
</cellWallAdhesion>

<sim>
    <tRelax>...</tRelax>
    <tRamp>...</tRamp>
    <tmax>...</tmax>
    <tmeas>...</tmeas>
    <tcheckpoint>...</tcheckpoint>
</sim>
```

两段参数互相独立。现有 cell-cell 基线可以迁移到新标签，但 cell-wall 的 `r0`、`rc`、`epsilon`、`D0`、`alpha` 尚无用户确认值，本方案不替用户假定；也不因函数形式相同而复制 cell-cell 数值。两个细胞、短程强非线性力且计算量很小，第一版将 material、particle velocity 以及两类 adhesion 的共用 timescale separation 都设为 `1`；等稳定性和收敛性确认后再考虑优化。

`RBC.xml` 保持现有高阶 RBC 本构参数和 `radius=3.91e-6 m` 不变；`config.xml` 中的 `<ibm><radius>` 必须与其一致。cell-wall LJ 分支已经承担近壁排斥，因此 `twoCellShear` 不再额外调用旧 `enableBoundaryParticles(...)`，避免对同一膜节点—壁面格点重复叠加两套排斥模型。旧接口只保留给其他算例回归使用。

## 17. `twoCellShear.cpp` 实施顺序

后续实现按以下顺序进行：

1. 从配置读取 `Lx/Ly/Lz/dx/dt/shearrate`，计算 `32×32×40` 网格并调用 `param::lbm_shear_parameters(*cfg,nz)`；日志同时打印名义尺寸、实际壁间距和目标上壁速度。
2. 校验 MPI size 为 2，显式创建 `1×1×2` 的规则 block 管理，再构建 HemoCell lattice；不依赖默认自动分解。
3. 在算例内建立 `x/y` 周期、`z` 非周期的上下 velocity boundary，将两壁初始速度设为零，并将整个流体初始化为静止平衡态。
4. 初始化 cell field，注册 `RbcHighOrderModel("RBC", RBC_FROM_SPHERE)`，从 `<ibm>` 设置 material/particle timescale。
5. 删除 `hemocell.setCellFixed(0)` 及其注释，不再按 cellId 改写粒子速度或位置推进。`cellId=0/1` 只保留为 lower/upper 的识别标签。
6. 从 `<cellCellAdhesion>` 读取五个参数并调用 `hemocell.setAdhesion(...)`；从 `<cellWallAdhesion>` 独立读取另五个参数并调用 `hemocell.setBoundaryAdhesion(...)`。两个调用都在 LBM 参数有效、上下壁 dynamics 已建立且 cell field 已初始化之后执行，adhesion timescale 设为 `1`。
7. `setBoundaryAdhesion(...)` 每次启动都重新构建本地 `boundaryParticles`，包括 checkpoint restart；该静态壁面列表不加入 checkpoint 粒子序列化。实现要保证重复初始化不会向同一 vector 重复追加壁面格点。
8. RBC 输出至少包含 `OUTPUT_POSITION`、`OUTPUT_TRIANGLES`、`OUTPUT_VELOCITY`、`OUTPUT_FORCE` 和 `OUTPUT_FORCE_REPULSION`。最后一项沿用旧名称，在本算例中保存 cell-cell 与 cell-wall LJ/Morse 力之和。
9. fresh run 时先完成可选的静止无粒子 fluid warmup，再读取两个 RBC；在粒子加载后、第一次 `iterate()` 前统计两个初始最近距离并立即写一次输出。checkpoint run 只加载 checkpoint，不能重新读 `RBC.pos`。
10. 主循环每一步先根据 `iter/tRelax/tRamp` 设置本步 top 速度，再调用一次 `hemocell.iterate()`；按 `tmeas` 输出、按 `tcheckpoint` 保存 checkpoint。
11. 日志遍历两个 cell ID，记录阶段、迭代数、两个中心、包围盒、面积、体积、中心相对位移、最小跨细胞节点距、下 RBC 到下壁的最小离散格点距和当前上壁速度。标签改为 `cellId=0 lower/wall-adhering`、`cellId=1 upper`，不能再写 fixed/free。
12. 日志统一使用 `TwoCellShear`，结束时释放 boundary-condition 对象并正常关闭自定义日志文件。

预松弛“稳定”第一版仍由固定 `tRelax` 控制，不增加复杂自动收敛判据。调参时同时观察两细胞中心距、跨细胞最近距离、下 RBC 到壁面最近距离、两个 RBC 的面积/体积/包围盒、速度以及合并 interaction force 的平台区间。

算例实现继续优先参考和复用 `src_hemocell/examples` 中已有写法：`oneCellShear` 的剪切参数和输出流程、`cube` 的显式上下壁边界、现有 `applyBoundaryRepulsionForce()` 的壁面粒子流程、`cellCollision_interior_viscosity` 的双细胞流程，以及其他算例的 timescale、checkpoint 和 CMake 接入方式。不新建 interaction manager、自动收敛控制器、配置封装类或大量兜底分支。

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
- 初始输出确认恰好有两个不同 `cellId`，盘面平行、`x/y` 同心，第一行 RBC 靠近下壁且两个 RBC 都未被边界加载逻辑删点。
- 初始最小跨细胞节点距接近 `cellCellAdhesion/r0`；下 RBC 到下壁离散格点的最小距离接近 `cellWallAdhesion/r0`。两者都不得明显落入各自 LJ 的强排斥区。
- 日志分别打印两套参数的六个 LBM 内部量和 `ceil(rc_lbm)` 搜索层数，并确认 `particleEnvelope=20` 大于两者最大搜索范围。
- 确认 `boundaryParticles` 在 fresh run 和 checkpoint restart 中数量一致，重复初始化没有重复追加；在本算例几何下，下壁格点被正确识别。
- 初始时两个 cellId 的节点均执行正常流体速度插值和粒子推进路径，不存在 fixed-cell 分支。

### 19.2 静止粘附预松弛

- relaxation 全程 top/bottom 速度均为零，fluid 无外加体力。
- 下 RBC 允许发生有限平移、转动和形变；其靠壁膜节点距离应向 `cellWallAdhesion/r0` 附近调整，中心不能持续远离下壁或穿入壁面。
- 上 RBC 在 cell-cell adhesion 下同步调整；跨细胞最近距离应向 `cellCellAdhesion/r0` 附近演化，不出现明显初始重叠导致的 LJ 爆发。
- cell-cell 节点对仍等量反向；cell-wall 是外部势，粒子系统净力不要求为零。壁面锚固反力最终由静止壁边界承担。
- `OUTPUT_FORCE_REPULSION` 中同时存在 cell-cell 和 cell-wall 贡献，且同一步两项相加而非相互覆盖或跨步累计。
- 通过下 RBC 中心高度/横向漂移、两套最近距离、两个 RBC 的面积/体积/包围盒和速度平台来选择足够的 `tRelax`，不能再用“下 RBC 完全不动”作为判据。
- 远离 RBC 与靠近 RBC 的流体速度/受力场应有差异，确认正常可变形下 RBC 通过 IBM 对流体产生阻碍，而不是流体无反馈地穿过一个被运动学锁死的网格。

### 19.3 上壁启动与剪切阶段

- bottom 始终为零，top 按 ramp 到达 `uTopTarget`，日志中的 LBM/SI 壁速与配置剪切率一致。
- 无细胞基准或远离细胞区域的速度剖面最终接近 `u_x(z)=gammaDot*z`。
- 观察上、下 RBC 间是否保持粘附、滑移、滚动或解离，同时观察下 RBC 对下壁是保持锚定、有限滑移/滚动还是脱附。
- 下 RBC 可变形且可有有限位移；“固定效果”改用其是否持续位于 wall cutoff 内、中心和接触区是否进入有界状态判断，不要求节点坐标严格不变。
- 检查 `OUTPUT_FORCE_REPULSION` 总 interaction force、两细胞相对位移、下 RBC 相对下壁位移以及两者形变，确认壁面粘附没有误施加成同细胞内部作用。
- checkpoint 分别在 relaxation、ramp 和恒定剪切阶段重启一次，确认壁速阶段不会重置。
- checkpoint 前后两类 adhesion 均保持启用、boundary particle 列表正确重建、轨迹连续，不恢复任何 fixed-cell 状态。

### 19.4 稳定性和参数扫描

- 首先以较短 `tmax`、高频 `tmeas` 做试运行，再增加总时长。
- 先固定 cell-cell 参数，只扫描 cell-wall 的 `epsilon/D0/alpha` 及必要时的 `r0/rc`；壁面锚定稳定后再单独扫描 cell-cell 参数。不要同时改变网格、timescale separation、两套粘附参数和 `dt`。
- 检查 `uTopLbm` 保持低 Mach 数；当前示例值约 `2.16e-4`，远低于常用稳定性上限。
- 检查膜节点是否进入壁面或被 `advanceParticles()` 删除；若发生，优先检查初始 wall 距离、wall LJ 参数和 `dt`，不叠加旧 boundary repulsion 掩盖问题。
- 上 RBC 若进入上壁 cutoff，也会按当前通用 wall 模型受力；记录这一事件。若后续物理目标明确要求“仅下壁粘附”，再单独扩展边界区域选择，不在第一版内核中硬编码 `z=0`。

## 20. `twoCellShear` 完成判据

只有同时满足以下条件，算例改造才算完成：

- `RBC.pos` 含两个平行 RBC，初始 cell-cell 最近距离接近 `cellCellAdhesion/r0`，下 RBC 到下壁离散格点的最近距离接近 `cellWallAdhesion/r0`，两处都没有明显初始重叠；
- `setCellFixed(...)` 及其源码字段、插值/推进分支和 main 调用已删除，`cellId=0` 与 `cellId=1` 都按正常 RBC 演化；
- 下 RBC 可以变形并对流体产生可观测反馈，在 relaxation 和剪切阶段由 cell-wall adhesion 保持在下壁附近；允许有限滑移、滚动和形变，不再要求逐节点位置与速度恒定；
- 域、周期性、单边壁速和 `1×1×2` MPI 分块与本方案一致；
- 五个 cell-cell 参数全部由 `<cellCellAdhesion>` 读取并传给 `setAdhesion(...)`，五个独立 cell-wall 参数全部由 `<cellWallAdhesion>` 读取并传给 `setBoundaryAdhesion(...)`，无算例内硬编码或跨段回退；
- cell-wall 内核的 LJ 段排斥、Morse 段吸引、cutoff 外为零，wall cutoff 大于 1 LU 时动态搜索不漏对；
- cell-cell 与 cell-wall 力在同一 `force_repulsion` 槽正确相加并铺展到流体，输出重算和正常迭代一致；
- relaxation 推进完整 HemoCell 耦合过程，两壁静止，两套粘附距离、两个 RBC 形变和速度进入可辨识的平台状态；
- 剪切阶段只有上壁运动，壁速严格由剪切率、`dt` 和离散壁间距计算；
- fresh run 和三个阶段的 checkpoint restart 均行为一致；
- 两个 cell ID、两套最近距离、下 RBC 相对壁面位移及总 interaction force 都进入日志或输出；
- 两核运行无跨 block 漏力、重复力、NaN/Inf 或粒子被意外删除；
- 旧 `setRepulsion()`、`enableBoundaryParticles()` 和 `applyBoundaryRepulsionForce()` 回归通过，未修改其他具体算例；
- `compile.sh` 能按“`src_hemocell/build/libhemocell.a` → 独立链接 twoCellShear”的两阶段流程工作，最终可执行文件位于 `twoCellShear/twoCellShear`，构建和结果文件均不进入 Git。

## 21. 实施前待用户确认项

1. 新公共 setter 暂定命名为 `setBoundaryAdhesion(...)`，底层函数按要求命名为 `applyBoundaryAdhesionForce()`。
2. 第一版 `<cellWallAdhesion>` 对 `populateBoundaryParticles()` 识别出的所有固体壁面生效；通过初始几何让下 RBC 只与下壁接触。若物理上必须永久禁止上壁粘附，需要在实施前改为只选择 bottom boundary。
3. cell-cell 与 cell-wall adhesion 共用 `adhesionTimescale=1` 和 `force_repulsion` 输出槽，`OUTPUT_FORCE_REPULSION` 是两类力之和；第一版不提供分项输出。
4. `<cellWallAdhesion>` 的五个数值尚待用户给定。最终 `RBC.pos` 高度必须在这些数值确定后，依据离散膜节点—壁面格点的真实最近距离校准。
5. 新策略下“固定”表示下 RBC 在观察时段内由壁面粘附保持在壁面附近，而不是节点严格零位移。允许的中心漂移、滑移和形变量需在参数试运行后确定验收阈值。
