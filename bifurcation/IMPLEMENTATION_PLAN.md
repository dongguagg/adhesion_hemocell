# 分叉管道中红细胞粘附流动算例改造与实现计划

## 1. 文档目的

本文档用于规划 `bifurcation` 算例从“包含 RBC、PLT、WBC 的预入口分叉流动算例”改造为“仅包含可相互粘附 RBC、并与管壁保持 HemoCell 原版排斥作用的分叉管道流动算例”。

本文档只描述待实施的代码、配置、构建和验证工作。本轮不修改：

- `bifurcation.cpp`；
- `config.xml`；
- `bifurcation3D.stl`；
- 用户已提供的 `RBC.pos`；
- HemoCell 的粘附力和原版壁面排斥力公式。

用户提供的 `RBC.pos` 作为既定输入文件使用，不再由本任务调用 `packCells` 重新生成，也不在构建或运行脚本中覆盖它。

本计划已经确认以下两项总体方案：

- 流体采用“周期 pre-inlet 内体积力驱动 → 向主分叉域提供速度入口 → 两个 daughter outlet 使用压力边界”的方式；主分叉域内部不施加统一体积力；
- HemoCell 源码增加 `enum class BoundaryParticleSelection` 和 `enableBoundaryParticles(...)` 四参数重载，由 `bifurcation` 显式选择 `SolidBounceBackOnly`，使开放入口/出口不参与 RBC-wall 排斥，同时保留旧三参数接口和旧行为。

## 2. 已确认的当前状态

### 2.1 算例文件

当前 `bifurcation` 目录包含：

| 文件 | 当前用途和状态 |
| --- | --- |
| `bifurcation.cpp` | 主程序，基于 HemoCell 的 STL pre-inlet 算例改写 |
| `bifurcation3D.stl` | 分叉管道几何，不计划修改 |
| `config.xml` | 流体、预入口、IBM 和模拟时间配置，尚无五参数粘附段 |
| `RBC.pos` | 用户提供的初始 RBC 位置，保持原样 |

目录中目前缺少 `RBC.xml`、独立构建文件和运行脚本。实施时直接复制 examples 中常用的标准高阶 RBC 配置：

```text
src_hemocell/examples/pipeflow_with_preinlet/RBC.xml
    -> bifurcation/RBC.xml
```

该文件与 `pipeflow`、`oneCellShear`、`cube`、`parallelplanes` 等常用算例中的 `RBC.xml` 内容一致。复制后保持其膜本构参数原样，不为本算例单独修改 `eta_m`、`kBend`、`kVolume`、`kArea`、`kLink`、网格分辨率、半径或参考体积。

### 2.2 `RBC.pos` 快照

已对文件做只读检查：

- 首行声明 RBC 数：`139`；
- 实际位置记录数：`139`；
- `x` 坐标范围约为 `0.0006–34.9750 um`；
- `y` 坐标范围约为 `0.0065–34.9992 um`；
- `z` 坐标范围约为 `0.0029–34.9966 um`；
- 每条记录均包含三个位置和三个欧拉角。

文件格式本身完整，位置覆盖约 `35 um × 35 um × 35 um`。本项目直接接受这份文件作为最终初始条件，不再依据 pre-inlet 的几何体积重新计算或校验压积，也不根据加载后的细胞数量调整代码、配置或位置文件。实施时不得自动缩放、平移、裁剪或重写 `RBC.pos`。

### 2.3 当前细胞和相互作用状态

当前主程序注册了：

- `RBC`；
- `PLT`；
- `WBC`。

原来的 `setRepulsion(...)` 及其时间尺度调用已经被注释，因此当前算例实际上没有启用细胞间排斥。HemoCell 源码中已经存在可直接调用的：

```cpp
setAdhesion(r0, rc, epsilon, D0, alpha)
```

现有粘附内核只处理不同 `cellId` 的膜节点对，并实现短程 LJ 排斥、中程 Morse 吸引和 cutoff 外零力。此次不重新实现或改写该公式。

原版 RBC-wall 排斥继续使用：

```cpp
enableBoundaryParticles(boundaryRepulsionConstant,
                        boundaryRepulsionCutoff,
                        timestep)
```

及现有 `applyBoundaryRepulsionForce()` 公式，不调用 `setBoundaryAdhesion(...)`。实施时新增带 `BoundaryParticleSelection` 的四参数重载；上述旧三参数形式继续调用原版“所有 `isBoundary()` dynamics”模式，保证其他算例源码兼容。

## 3. 目标物理模型

完成后的算例满足：

1. 只有一种细胞类型 `RBC`；
2. 不加载、不注册、不输出 PLT 或 WBC；
3. 不同 RBC 之间使用五参数 LJ/Morse 粘附模型；
4. 同一 RBC 内部节点不施加细胞间粘附；
5. RBC 与真正的固体管壁之间使用 HemoCell 原版边界排斥模型；
6. 入口 velocity boundary 和出口 pressure boundary 不参与壁面排斥；
7. RBC 通过周期 pre-inlet 持续进入主分叉域；
8. 两个 daughter outlet 都是开放压力出口，RBC 能正常离开；
9. 细胞粘附力、壁面排斥力和膜本构力共同铺展到流体；
10. 周期 pre-inlet 根据目标 `Re` 沿入口方向施加体积力，主分叉域只接受 pre-inlet 速度入口且不施加统一体积力；
11. 两个 daughter outlet 第一版采用相同参考压力，实际分流由几何、流动阻力和细胞分布自然决定；
12. checkpoint 重启后所有相互作用、流体驱动和边界选择保持一致。

## 4. 明确不做的改造

- 不生成或改写 `RBC.pos`；
- 不修改 STL 几何，除非体素化验证证明第二个出口无法正确打开，并经用户另行确认；
- 不启用 RBC-wall adhesion；
- 不删除或替换 HemoCell 原版 `setRepulsion(...)`、`enableBoundaryParticles(...)` 或 `applyBoundaryRepulsionForce()`；
- 不修改 Palabos 中 `Dynamics::isBoundary()`、`BoundaryCompositeDynamics::isBoundary()` 或 Zou/He 边界语义；
- 不在整个主分叉域施加固定方向的统一体积力，也不为两个 daughter branches 分别硬编码空间体积力；
- 不给 RBC 增加显式永久键、不可逆聚集或自定义 bond 数据；
- 不修改膜本构、IBM 插值/铺展、流体碰撞模型、粒子序列化或现有输出文件格式；
- 不修改其他算例；
- 不在本轮顺带进行通用 interaction manager 或大范围框架重构。

## 5. `bifurcation.cpp` 改造计划

### 5.1 删除 PLT/WBC

删除以下头文件：

```cpp
#include "pltSimpleModel.h"
#include "wbcHighOrderModel.h"
```

删除 PLT/WBC 的：

- `addCellType`；
- material timescale；
- `setInitialMinimumDistanceFromSolid`；
- `setOutputs`；
- 统计和日志标签。

只保留：

```cpp
hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
```

### 5.2 调整初始化顺序

初始化顺序改为：

1. 读取并体素化 STL；
2. 设置 LBM 参数；
3. 创建 pre-inlet 和主域 lattice；
4. 建立 STL bounce-back 管壁；
5. 建立主域入口 velocity boundary；
6. 识别两个真实 daughter outlet，并只在出口流体格点上建立 pressure boundary；
7. 设置初始流体平衡态和 driving force；
8. 调用 `lattice->initialize()`；
9. 初始化 cell field；
10. 只注册 RBC；
11. 配置 material、particle 和 interaction timescale；
12. 启用 RBC-RBC adhesion；
13. 以 `SolidBounceBackOnly` 模式启用原版 boundary repulsion，在生成 boundary-particle 列表时排除开放入口和出口；
14. fresh run 时完成无细胞流体 warmup，再直接读取既有 `RBC.pos`；
15. checkpoint run 只加载 checkpoint，不重新读取 `RBC.pos`；
16. 写入初始输出和初始化统计；
17. 进入主循环。

压力边界必须在 `lattice->initialize()` 和边界粒子构建之前完成，避免当前代码中“lattice 初始化后才设置出口”的顺序。

### 5.3 已确认的流体驱动方案

本算例保留原代码的总体 pre-inlet 架构，但修正出口、warmup 和时间推进顺序：

1. pre-inlet 是沿入口方向周期的直管，仅在该辅助管段中使用均匀体积力；
2. `calculateDrivingForce()` 根据 `<preInlet><parameters><Re>`、流体运动黏度和入口截面估算驱动力；
3. `Direction::Xneg` 对应 pre-inlet 沿全局 `+x` 方向驱动；
4. 主分叉域不设置全域 external force；
5. pre-inlet 截面的瞬时速度通过现有 `applyPreInletVelocityBoundary()` 写入主域 velocity inlet；
6. RBC 继续通过现有 `applyPreInletParticleBoundary()` 从周期 pre-inlet 传入主域；
7. 两个 daughter outlet 分别建立 pressure boundary，第一版共同使用 `outletDensity=1.0` 的 LBM 参考密度；
8. 不指定两个出口流量，满足 `Qin = Qout1 + Qout2` 的分流由两个分支的几何阻力和瞬时细胞分布自然形成。

固定驱动力必须在对应的 pre-inlet 流体推进使用它之前写入 lattice。fresh warmup 和正式主循环都要持续执行驱动力设置与 pre-inlet 速度交换，不能保留旧代码“warmup 只调用 `collideAndStream()`”的行为。对于现有一步耦合顺序，可以在推进前设置本步 pre-inlet force，在推进后把更新后的入口速度和粒子传给主域供下一步使用；启动时应先完成一次入口速度初始化，避免主域首步读取未设置的入口状态。

当前 `Re=0.1` 仍表示按无细胞圆管关系估算的名义 pre-inlet 驱动力。第一版只记录实际 `Qin` 和实际流速，不因 RBC 引起的表观黏度变化自动反馈调整驱动力；若后续要求严格恒定流量，再另立任务增加流量反馈控制。

### 5.4 读取并启用 RBC-RBC adhesion

从 `<cellCellAdhesion>` 独立读取：

```text
r0       [um]
rc       [um]
epsilon  [J]
D0       [J]
alpha    [um^-1]
```

随后调用：

```cpp
hemocell.setAdhesion(r0, rc, epsilon, D0, alpha);
```

不调用旧 `setRepulsion(...)`，也不为 adhesion 提供从旧 `kRep`/`RepCutoff` 自动回退的逻辑。

### 5.5 启用原版 RBC-wall 排斥

壁面排斥参数直接采用 `src_hemocell/examples` 中常见的基线值，并将当前配置中的同组数值迁移到语义明确的 `<boundaryRepulsion>` 段：

```xml
<boundaryRepulsion>
    <constant>2e-22</constant>
    <cutoff>0.7</cutoff> <!-- um -->
</boundaryRepulsion>
```

其中 `constant=2e-22`、`cutoff=0.7 um` 分别参考 examples 常见的 `kRep` 和 `RepCutoff`。本算例不再为这两个参数保留占位符，也不等待另行调参后才实施。

通过新增四参数重载调用原版排斥链：

```cpp
hemocell.enableBoundaryParticles(
    constant,
    cutoff,
    interactionEvery,
    BoundaryParticleSelection::SolidBounceBackOnly);
```

旧三参数调用仍然可用并保持 `AllBoundaryDynamics` 行为。壁面排斥的力公式保持不变；固壁候选格点必须在 STL 管壁、主域速度入口和两个压力出口全部建立完成后构建。

### 5.6 统一 interaction 更新时间

RBC-RBC adhesion 和 RBC-wall repulsion 共用 `sv.force_repulsion`。当前调度中：

- adhesion 更新会先清零 `force_repulsion`；
- boundary repulsion 只向其中累加，不自行清零。

因此两者必须使用同一个 `interactionEvery`。第一版建议：

```xml
<stepMaterialEvery>1</stepMaterialEvery>
<stepParticleEvery>1</stepParticleEvery>
<stepInteractionEvery>1</stepInteractionEvery>
```

待稳定性确认后才能提高 separation。若提高，则必须满足：

```text
adhesionTimescale == boundaryRepulsionTimescale
interactionEvery % stepParticleEvery == 0
stepMaterialEvery % stepParticleEvery == 0
```

`tmeas` 最好也是 `interactionEvery` 的整数倍，使输出重算与正常迭代的更新时间一致。

### 5.7 用 `BoundaryParticleSelection` 区分固壁和开放边界

当前 `populateBoundaryParticles()` 使用 `dynamics.isBoundary()` 判断壁面。以下三者都会返回 true：

- STL 管壁的 bounce-back dynamics；
- 主域入口的 velocity-boundary dynamics；
- 主域出口的 pressure-boundary dynamics。

若不筛选，入口和出口会产生非物理壁面排斥，阻止 RBC 进入或离开。

实施采用已经确认的源码方案：

```cpp
enum class BoundaryParticleSelection {
    AllBoundaryDynamics,
    SolidBounceBackOnly
};
```

1. 保留原三参数 `enableBoundaryParticles(...)` 的声明和实现入口，由它转发到 `AllBoundaryDynamics`，保证源码和静态库重编译后的旧调用行为不变；
2. 新增接收第四个 `BoundaryParticleSelection` 参数的重载；
3. `AllBoundaryDynamics` 完整复用当前 `dynamics.isBoundary()` 候选和邻域判断；
4. `SolidBounceBackOnly` 通过 `BounceBack<T,DESCRIPTOR>` 的 dynamics 类型或唯一 `getId()` 判断真实固壁，不把一般的 `BoundaryCompositeDynamics` 当作固壁；
5. solid-only 模式的当前格点和邻域格点都使用同一个 `isSolidWallDynamics(...)` 判据：当前格点必须是 solid wall，且至少有一个邻点不是 solid wall，才加入 `boundaryParticles`；
6. 不修改 Palabos 的 `isBoundary()` 返回值，因为 velocity/pressure boundary 仍必须保持合法的 LBM 边界语义；
7. 不修改 `applyBoundaryRepulsionForce()` 的公式、搜索范围或累加方式；
8. pre-inlet 挤出圆管侧壁和 STL 管壁当前均为 BounceBack，必须保留；主域 velocity inlet 和两个 pressure outlets 属于复合开放边界，必须排除；
9. fresh run 和 checkpoint restart 每次都用相同 selection 重新构建该静态列表；若以后启用 `doRestructure()`/动态负载均衡，重构后也必须重新构建；
10. 启动日志分别记录候选 boundary dynamics 数、最终 solid boundary-particle 数，以及被排除的 velocity/pressure boundary 数。

最终分类必须满足：

```text
solid bounce-back wall  -> 保留
velocity inlet          -> 排除
pressure outlet         -> 排除
```

不能只采用“排除所有 `BoundaryCompositeDynamics`”的黑名单，因为未来某些复合 dynamics 可能表示合法固壁。当前模式采用明确的 BounceBack 固壁白名单；若后续引入 moving wall、`VelocityBounceBack` 或其他固壁类型，再显式扩展 solid-wall 判据。长期若需要支持任意壁面模型，可另立任务增加由几何提供的 solid mask，不在本轮同时引入第二套选择机制。

### 5.8 修正双出口压力边界

当前代码对 `x=bb.x1` 的整个矩形 `y-z` 平面施加 pressure boundary，可能覆盖管腔外的 bounce-back 固体。

实施时改为：

1. 在体素化 flag 中寻找与下游端面相连的流体格点；
2. 对出口流体格点逐点设置 `addPressureBoundary0P`；
3. 不修改同一平面上的固体格点；
4. 将出口流体格点划分为两个连通分量，分别记录数量和中心；
5. 若两个出口不在同一个离散 `x` 层，分别在各自真实终止层设置压力边界；
6. 第一版两个出口分别设置相同的 `outletDensity=1.0`，不额外规定分流比；
7. 启动时打印两个 outlet 的格点数、中心、法向选择和参考密度，任一出口为空都立即报错退出。

STL 中两个出口的连续几何端面约位于 `x=152.19` 和 `x=152.65`，相差不足一个当前 LU，但仍需依据实际体素化结果判断，不能只依赖连续 STL 坐标。

### 5.9 处理 fresh run、warmup 和 checkpoint

fresh run：

1. 完成所有流体边界；
2. 在没有粒子的情况下执行可选 fluid warmup；
3. warmup 每一步在流体推进使用驱动力前设置 pre-inlet driving force，并同步推进、交换入口速度；
4. warmup 结束后加载现有 `RBC.pos`；
5. 依赖 HemoCell 现有加载流程读取细胞，不增加压积校验或调整逻辑；
6. 写一次 `iter=0` 输出；
7. 开始完整 HemoCell 迭代。

checkpoint run：

1. 从 checkpoint XML 读取同一套配置；
2. 不再次读取 `RBC.pos`；
3. 恢复 pre-inlet 和主域 lattice/particle fields；
4. 使用 `BoundaryParticleSelection::SolidBounceBackOnly` 重新构建静态 solid boundary-particle 列表；
5. 在恢复后的下一次流体推进前重新设置 driving force，并恢复 pre-inlet 速度交换；
6. 从已保存迭代数继续运行。

## 6. `config.xml` 改造计划

建议配置结构为：

```xml
<?xml version="1.0" ?>
<hemocell>

<preInlet>
  <parameters>
    <lengthN>60</lengthN>
    <Re>0.1</Re>
  </parameters>
</preInlet>

<flow>
    <outletDensity>1.0</outletDensity> <!-- LBM reference density for both outlets -->
</flow>

<parameters>
    <warmup>...</warmup>
    <outputDirectory>output</outputDirectory>
    <checkpointDirectory>checkpoint</checkpointDirectory>
    <logDirectory>log</logDirectory>
    <logFile>bifurcation.log</logFile>
</parameters>

<ibm>
    <stepMaterialEvery>1</stepMaterialEvery>
    <stepParticleEvery>1</stepParticleEvery>
    <stepInteractionEvery>1</stepInteractionEvery>
    <initialMinimumDistanceFromSolid>0.5</initialMinimumDistanceFromSolid>
</ibm>

<domain>
    <geometry>bifurcation3D.stl</geometry>
    <fluidEnvelope>2</fluidEnvelope>
    <rhoP>1025</rhoP>
    <nuP>1.1e-6</nuP>
    <dx>5e-7</dx>
    <dt>1e-7</dt>
    <refDir>2</refDir>
    <refDirN>40</refDirN>
    <blockSize>-1</blockSize>
    <particleEnvelope>35</particleEnvelope>
    <kBT>4.100531391e-21</kBT>
</domain>

<cellCellAdhesion>
    <r0>...</r0>
    <rc>...</rc>
    <epsilon>...</epsilon>
    <D0>...</D0>
    <alpha>...</alpha>
</cellCellAdhesion>

<boundaryRepulsion>
    <constant>2e-22</constant>
    <cutoff>0.7</cutoff> <!-- um -->
</boundaryRepulsion>

<sim>
    <tmax>...</tmax>
    <tmeas>...</tmeas>
    <tcheckpoint>...</tcheckpoint>
</sim>

</hemocell>
```

说明：

- `<preInlet><parameters><Re>` 只用于周期直管 pre-inlet 的名义驱动力计算，主分叉域不读取或施加统一 external force；
- `<flow><outletDensity>` 是两个 daughter outlets 共同使用的 LBM 参考密度，第一版固定为 `1.0`，两个出口仍分别建立和统计；
- `<cellCellAdhesion>` 五个参数必须由用户确认，不从 `twoCellShear` 自动复制；
- `<boundaryRepulsion>` 参数独立于 adhesion，并固定采用 examples 常见基线 `constant=2e-22`、`cutoff=0.7 um`；
- 当前 `<domain><kRep>` 和 `<domain><RepCutoff>` 删除或迁移，避免继续保留错误的 cell-cell 语义；
- `<parameters><maxPackIter>` 删除，因为本任务不再运行 `packCells`；
- `<maxFlin>`、`<tbalance>` 只有在确实启用 PARMETIS load balancing 时才保留；
- 配置中不增加 `targetHematocrit` 或其他有效压积校验参数。

## 7. `RBC.pos` 的直接使用

### 7.1 文件保护

构建和运行脚本必须满足：

- 不调用 `packCells`；
- 不删除 `RBC.pos`；
- 不覆盖 `RBC.pos`；
- 不在运行前自动裁剪或转换其坐标；
- 不因运行时实际加载数量而回写或重新生成位置文件。

### 7.2 使用原则

由于这是 pre-inlet 算例，`RBC.pos` 在周期预入口中加载，然后 RBC 被持续复制/传递到主分叉域。实现直接沿用这一加载机制，不新增下列行为：

- 不计算所谓“有效初始压积”；
- 不比较 `RBC.xml` 与 `packCells` 的单细胞体积口径；
- 不判断加载后的压积是否仍为 `0.2`；
- 不因压积偏离 `0.2` 而报警、停止或调整模拟；
- 不为了匹配目标压积修改 pre-inlet 尺寸或细胞数量。

正常运行日志仍可保留 HemoCell 原有的 RBC 数量信息，用于判断程序是否正常演化，但该数量不用于压积验收。

## 8. 输出和运行诊断

### 8.1 RBC 输出

至少启用：

```cpp
OUTPUT_POSITION
OUTPUT_TRIANGLES
OUTPUT_VELOCITY
OUTPUT_FORCE
OUTPUT_FORCE_REPULSION
OUTPUT_FORCE_VOLUME
OUTPUT_FORCE_BENDING
OUTPUT_FORCE_LINK
OUTPUT_FORCE_AREA
OUTPUT_FORCE_VISC
OUTPUT_CELL_ID
OUTPUT_VERTEX_ID
OUTPUT_RES_TIME
```

其中 `OUTPUT_FORCE_REPULSION` 保存 RBC-RBC adhesion 与 RBC-wall repulsion 的合力。第一版不增加两个独立粒子力字段；如后续需要分项输出，应另立任务。

### 8.2 流体输出

保留：

```cpp
OUTPUT_VELOCITY
OUTPUT_DENSITY
OUTPUT_FORCE
OUTPUT_BOUNDARY
```

### 8.3 日志

按 `tmeas` 输出：

- 当前迭代数和物理时间；
- pre-inlet RBC 数；
- 主域 RBC 数；
- 当前 pre-inlet driving force；
- 主入口流体体积流量 `Qin`；
- 主入口 RBC 通量；
- 两个 daughter outlet 各自的流体体积流量；
- 质量守恒残差 `Qin-(Qout1+Qout2)`；
- 两个出口各自的 RBC 数通量；
- 最大流体速度和低 Mach 检查；
- RBC 最大速度；
- 总 interaction force 最大值；
- 被删除的不完整 RBC 数；
- 可选的粘附聚集体数量和最大 cluster size。

聚集体定义第一版可在后处理中依据不同 `cellId` 间是否存在 `r < rc` 的膜节点对建立连通图，不必在主求解器中保存显式 bond。

## 9. HemoCell 源码最小改造范围

细胞间 adhesion 已存在，因此原则上不修改其 setter、单位换算或节点对公式。本算例需要的公共源码改动仅限于正确筛选固体壁面格点，预计涉及：

| 文件 | 计划职责 |
| --- | --- |
| `src_hemocell/core/hemoCellParticleField.h` | 在 `hemo` 命名空间声明 `BoundaryParticleSelection`，并增加接收 selection 的 `populateBoundaryParticles(...)` 重载 |
| `src_hemocell/hemocell.h` | 保留旧三参数 `enableBoundaryParticles(...)`，新增接收 `BoundaryParticleSelection` 的四参数重载 |
| `src_hemocell/core/hemoCell.cpp` | 旧三参数入口转发到 `AllBoundaryDynamics`；四参数入口保存排斥参数并按 selection 构建列表 |
| `src_hemocell/core/hemoCellFields.h/.cpp` | 让 `HemoPopulateBoundaryParticles` functional 携带 selection，并在各 MPI block 上按同一模式生成列表 |
| `src_hemocell/core/hemoCellParticleField.cpp` | 保留原模式；新增 BounceBack solid-wall 判据和 `SolidBounceBackOnly` 候选/表面判断；不改排斥力公式 |

计划接口为：

```cpp
enum class BoundaryParticleSelection {
    AllBoundaryDynamics,
    SolidBounceBackOnly
};

void enableBoundaryParticles(T constant, T cutoff,
                             unsigned int timestep = 1);

void enableBoundaryParticles(T constant, T cutoff,
                             unsigned int timestep,
                             BoundaryParticleSelection selection);
```

保留真正的三参数重载，而不是只把原函数直接改成带默认第四参数的新签名。这样旧调用路径和职责清楚，重新编译静态库后其他算例无需改源码；四参数重载只由本算例显式调用。

兼容要求：

- 旧三参数 `enableBoundaryParticles(...)` 继续可用；
- 旧三参数入口严格转发到 `BoundaryParticleSelection::AllBoundaryDynamics`；
- 旧 `applyBoundaryRepulsionForce()` 数学表达式不变；
- 旧 boundary adhesion 不受影响；
- 其他算例不自动切换到新的 solid-only 筛选行为；
- `bifurcation` 显式调用四参数入口和 `SolidBounceBackOnly`；
- `SolidBounceBackOnly` 的表面判断同时检查当前格点与邻域格点是否为 solid wall，不再用广义 `isBoundary()` 判断邻域是否为流体侧；
- 不更改 velocity/pressure boundary 的 dynamics，也不在列表生成后依赖算例端再次擦除坐标。

## 10. 原版壁面排斥的适用范围检查

原版 `applyBoundaryRepulsionForce()` 固定只访问壁面格点附近 `±1` 个 particle-grid bin。第一版需检查：

```text
boundaryRepulsionCutoffLbm = cutoff_um * 1e-6 / dx
```

本计划固定使用的 `cutoff=0.7 um`、`dx=0.5 um` 对应 `1.4 LU`，与原搜索范围相容。若后续任务要求把 cutoff 增大到可能跨越两个以上 bin：

- 不能假设原内核仍能找全作用对；
- 应先报告该限制；
- 若要扩展动态搜索范围，必须另行确认，因为这超出“保持 HemoCell 原版壁面排斥函数”的当前边界。

精确零距离会使原版公式发生除零。由于本轮要求保留公式，不增加任意力钳制；通过初始最小离壁距离、时间步和运行时有限值检查避免膜节点落到壁面格点上。

## 11. 构建与运行接入

`bifurcation` 位于仓库根目录，构建和运行脚本采用与 `twoCellShear` 相同的简洁形式：

1. `src_hemocell/build/libhemocell.a` 由 HemoCell 自身构建流程预先生成；
2. `compile.sh` 只在 `bifurcation/build` 配置并构建当前算例；
3. 最终可执行文件输出为 `bifurcation/bifurcation`；
4. `compile.sh` 和 `run.sh` 均从 `bifurcation` 目录调用，与 `twoCellShear` 的脚本使用方式一致；
5. `run.sh` 固定以 16 个 MPI rank 启动 `config.xml`；
6. 构建脚本不运行 `packCells`；
7. 构建目录、可执行文件、输出、checkpoint 和日志加入 `.gitignore`。
8. 固定使用 16 个 MPI rank：pre-inlet 按 `1×2×2` 使用 4 rank，主域按 `3×2×2` 使用 12 rank；程序启动时校验配置分块数之和与实际 MPI rank 数一致。

新增算例文件计划：

| 文件 | 用途 |
| --- | --- |
| `bifurcation/RBC.xml` | 原样复制 `src_hemocell/examples/pipeflow_with_preinlet/RBC.xml` |
| `bifurcation/CMakeLists.txt` | 独立链接 HemoCell 静态库 |
| `bifurcation/compile.sh` | 两阶段增量构建 |
| `bifurcation/run.sh` | 明确 MPI 进程数和配置路径 |

## 12. 分阶段实施顺序

### 阶段 A：建立可构建的纯 RBC 基线

1. 原样复制 examples 常用 `RBC.xml`，并补充独立构建文件；
2. 删除 PLT/WBC 注册和输出；
3. 将两个 daughter outlets 按真实流体连通分量分别建立为等参考压力边界；
4. 修正出口边界创建顺序、pre-inlet force 设置时序和 fluid warmup；
5. 验证只有周期 pre-inlet 受到体积力，主分叉域 external force 为零；
6. 只用 RBC、暂不启用 adhesion/boundary repulsion 做短程初始化；
7. 验证 pre-inlet 速度传递、主域流动、两个出口流量和 checkpoint 基线。

### 阶段 B：修正开放边界和固体壁面选择

1. 增加 `BoundaryParticleSelection`、旧三参数转发入口和四参数重载；
2. 在 `HemoCellFields` functional 与各 `HemoCellParticleField` block 中传递同一 selection；
3. 实现 `SolidBounceBackOnly` 候选和邻域表面判断；
4. 对主入口、两个出口、pre-inlet 侧壁和 STL 管壁分别计数；
5. 验证 `AllBoundaryDynamics` 保持旧行为，`SolidBounceBackOnly` 排除入口/出口；
6. 以 solid-only 模式启用原版 boundary repulsion；
7. 用单 RBC 通过入口、分叉和出口进行回归。

### 阶段 C：接入 RBC-RBC adhesion

1. 增加 `<cellCellAdhesion>`；
2. 读取五个参数并调用 `setAdhesion(...)`；
3. 统一 adhesion 和 boundary repulsion 的更新时间；
4. 增加 `OUTPUT_FORCE_REPULSION`；
5. 从极短、低输出间隔算例开始检查有限值和力符号。

### 阶段 D：接入现有 `RBC.pos`

1. 不改位置文件地完成加载；
2. 确认构建和运行脚本不会调用 `packCells`；
3. 确认 fresh run 读取该文件、checkpoint run 不重复读取；
4. 不增加有效压积计算、目标值比较或位置调整逻辑。

### 阶段 E：长程分叉流动验证

1. 检查流量充分发展；
2. 检查粘附聚集体形成、通过分叉和解离；
3. 统计两出口流体和 RBC 通量；
4. 检查 MPI block 边界处粘附力；
5. 检查 checkpoint 连续性；
6. 确认无 NaN/Inf、无入口/出口非物理反弹和异常细胞删除。

## 13. 验证矩阵

### 13.1 相互作用组合

| 模式 | RBC-RBC adhesion | RBC-wall repulsion | 目的 |
| --- | --- | --- | --- |
| A | 关 | 关 | 纯流固和边界基线 |
| B | 关 | 开 | 验证只有真实管壁排斥 |
| C | 开 | 关 | 验证粘附及 force slot 清零 |
| D | 开 | 开 | 最终目标及两类力正确累加 |

### 13.2 边界验证

- 管壁附近 RBC 受到指向管腔内部的排斥；
- pre-inlet 挤出侧壁和主域 STL 管壁在 `SolidBounceBackOnly` 中均被保留；
- 主入口附近不存在逆流向的伪壁面排斥；
- 两个出口附近不存在把 RBC 推回主域的伪壁面排斥；
- 旧 `AllBoundaryDynamics` 模式的列表与改造前一致；
- 两个 outlet 均有非零流量；
- 入口流量与两个出口流量之和在允许误差内一致；
- RBC 能从 pre-inlet 完整传入主域并从两个出口完整离开。

### 13.3 流体驱动验证

- `setDrivingForce()` 只对 `partOfpreInlet` 的周期辅助管段设置沿 `+x` 的力；
- 主分叉域所有普通流体格点的外加驱动力保持为零；
- fluid warmup 中 pre-inlet 流量从静止状态发展，主域入口同步收到非零速度；
- 两个出口采用相同参考密度，但分别计算 `Qout1` 和 `Qout2`；
- 稳态或统计稳定阶段满足 `Qin≈Qout1+Qout2`；
- 分流比不是由源码硬编码，而由分支阻力和瞬时细胞分布决定；
- checkpoint 后第一步恢复相同驱动力和入口交换，不出现重新从静止启动的瞬态。

### 13.4 adhesion 验证

- `r<r0` 时为排斥；
- `r=r0` 时为零；
- `r0<r<rc` 时为吸引；
- `r>=rc` 时为零；
- 不同 RBC 节点对满足等量反向；
- 同一 RBC 内部不产生该作用；
- interaction force 不跨时间步累计；
- MPI block/envelope 边界附近的节点对不漏算、不重复算。

### 13.5 回归和稳定性

- 默认 HemoCell 静态库可以构建；
- 原 `setRepulsion()` 接口仍可编译；
- 原三参数 `enableBoundaryParticles()` 的其他算例行为不变；
- 新四参数 `enableBoundaryParticles(..., SolidBounceBackOnly)` 在单进程和 MPI 下得到一致的固壁列表；
- `setBoundaryAdhesion()` 不被本算例启用但仍可编译；
- fresh run 和 checkpoint restart 结果连续；
- 所有粒子位置、速度和力保持有限；
- 没有因出口边界或壁面排斥错误造成的异常细胞删除。

## 14. 完成判据

只有以下条件全部满足，改造才算完成：

- 算例中只注册和输出 RBC；
- 使用用户提供的 `RBC.pos`，没有自动生成或改写；
- fresh run 直接读取包含 139 条记录的既有 `RBC.pos`；
- 五个 cell-cell adhesion 参数全部来自 XML；
- RBC-RBC 使用现有 LJ/Morse adhesion；
- RBC-wall 使用原版 boundary repulsion 公式；
- HemoCell 提供 `BoundaryParticleSelection::{AllBoundaryDynamics,SolidBounceBackOnly}`，旧三参数入口保持旧行为；
- `bifurcation` 通过四参数入口显式选择 `SolidBounceBackOnly`；
- inlet/outlet 不出现在 solid boundary-particle 列表中；
- 压力边界只作用于两个真实 outlet 的流体格点；
- 只有周期 pre-inlet 使用沿 `+x` 的体积力，主分叉域没有全域外力；
- 两个 daughter outlets 使用相同参考压力并分别统计流量，满足总体质量守恒；
- adhesion 与 boundary repulsion 使用相同刷新频率并正确累加；
- RBC 能从 pre-inlet 进入、通过分叉并从两个出口离开；
- 输出能够识别 cell ID、速度和合并 interaction force；
- 单进程基线和目标 MPI 配置下均无漏力、重复力、NaN/Inf 或异常删除；
- checkpoint 重启后轨迹、边界和相互作用连续；
- 构建和运行流程可重复，且不污染源码目录。

## 15. 实施前需要用户确认的项目

1. `<cellCellAdhesion>` 的 `r0`、`rc`、`epsilon`、`D0`、`alpha` 最终数值；
2. [已确认] 使用 16 个 MPI 进程，pre-inlet 为 `1×2×2`，主域为 `3×2×2`；
3. 是否需要在主程序中实时输出聚集体 cluster size，还是保留为基于 `cellId` 和距离的后处理任务。
