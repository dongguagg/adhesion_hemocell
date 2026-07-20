# 双红细胞光镊拉脱算例设计方案

## 1. 文档状态与实施边界

本文设计一个新的双红细胞光镊拉脱 benchmark，算例目录命名为
`twoCellPull`。第一版代码现已按本文方案实现；本文同时保留设计依据和实施后的接口说明。

新算例的 main 函数和 `config.xml` 位于：

```text
/home/jxh/adhesion_rbc/twoCellPull
```

可复用的表面节点选取与外力加载逻辑位于 HemoCell helper：

```text
src_hemocell/helper/hemoCellSurfaceForce.h
src_hemocell/helper/hemoCellSurfaceForce.cpp
```

本实现不修改膜本构、IBM 核函数、细胞间粘附、细胞-壁面粘附、粒子序列化或现有
粒子输出格式。外力复用现有 `sv.force`，不复用会被 adhesion 调度清零的
`sv.force_repulsion`。

## 2. 算例物理目标

算例包含两个相同的 RBC：

- `cellId=0` 是下方 RBC，通过 cell-wall adhesion 锚定在静止下壁附近；
- `cellId=1` 是上方 RBC，通过 cell-cell adhesion 与下方 RBC 粘附；
- 两个 RBC 都是正常可移动、可转动、可变形的 IBM 细胞；
- 上下壁均保持静止，本 benchmark 不施加 Couette 剪切；
- 先推进完整的细胞-流体耦合完成静止粘附预松弛，再施加光镊外力；
- 光镊外力可以施加到配置指定的任意 `cellId`；
- 光镊作用面可以从六个实验室坐标方向中选择；
- 实际拉力是用户配置的三维矢量，与作用面选择相互独立。

物理反力链为：

```text
光镊外力
    ↓
目标 RBC 的选定膜面补丁
    ↓
cell-cell adhesion
    ↓
另一 RBC
    ↓
cell-wall adhesion
    ↓
静止固体壁面
```

因此，双细胞粒子系统的外力合计不要求为零；最终反力由壁面边界条件承担。不要像
单细胞 `stretchCell` 那样在另一个细胞或相反表面人为施加一个等量反向力。

## 3. 计划中的目录和文件职责

```text
twoCellPull/
├── DESIGN.md                         # 本设计文档
├── twoCellPull.cpp                   # 算例 main
├── config.xml                        # 全部运行参数
├── RBC.pos                           # 两个 RBC 的初始位置
├── RBC.xml                           # RBC 本构参数
├── CMakeLists.txt                    # 独立链接 libhemocell.a
├── compile.sh                        # 两阶段增量构建
└── run.sh                            # MPI 启动脚本

src_hemocell/helper/
├── hemoCellSurfaceForce.h            # 公共接口和 functional 声明
└── hemoCellSurfaceForce.cpp          # 选点、MPI 汇总和施力
```

运行时额外生成：

```text
output/twoCellPull.csv                # 力、位移、补丁中心和细胞状态
output/grip_vertices.csv              # 选中 vertexId、初始位置、面积和权重
checkpoint/grip_state.dat             # checkpoint 对应的固定补丁状态
```

## 4. 配置接口

`config.xml` 中新增独立的 `<surfaceForce>` 段。建议固定使用以下接口：

```xml
<surfaceForce>
    <targetCellId>1</targetCellId>
    <face>zPositive</face>

    <!-- 实验室坐标系中的总拉力矢量，单位 pN。 -->
    <forceX>0.0</forceX>
    <forceY>0.0</forceY>
    <forceZ>50.0</forceZ>

    <!-- 第一版正式 benchmark 使用膜面测地半径，单位 um。 -->
    <patchMode>geodesicRadius</patchMode>
    <patchGeodesicRadius>0.70</patchGeodesicRadius>

    <!-- nodalArea 或 uniform；正式结果推荐 nodalArea。 -->
    <weighting>nodalArea</weighting>
</surfaceForce>

<sim>
    <tRelax>50000</tRelax>
    <tForceRamp>10000</tForceRamp>
    <tmax>300000</tmax>
    <tmeas>1000</tmeas>
    <tcheckpoint>50000</tcheckpoint>
</sim>
```

`<face>` 只接受以下六个区分大小写的枚举值：

| 配置值 | 表面搜索方向 |
| --- | --- |
| `xPositive` | `( 1, 0, 0)` |
| `xNegative` | `(-1, 0, 0)` |
| `yPositive` | `( 0, 1, 0)` |
| `yNegative` | `( 0,-1, 0)` |
| `zPositive` | `( 0, 0, 1)` |
| `zNegative` | `( 0, 0,-1)` |

表面方向只决定“在哪个膜面中心建立加载补丁”，`forceX/Y/Z` 决定“向哪里拉”。两者
不强制平行。例如可以在 `zPositive` 中央补丁上施加 `(+Fx,0,0)` 切向力。程序应打印
二者夹角；当夹角与常见的法向拉脱或切向加载不一致时给出提示，但不擅自改写用户矢量。

三个力分量同时为零是合法的零载荷回归模式。非零输入的单位换算只执行一次：

```text
force_lbm[d] = force_pN[d] * 1e-12 / param::df
```

## 5. 为什么不能用坐标极值选中央补丁

对于双凹 RBC，`zPositive` 表面的最大 `z` 位于环状隆起部，而不是中央凹陷。因此以下
规则不适用于本算例：

```text
选 z 最大的 N 个节点
选 z >= zmax-capDepth 的节点
```

`vertexId` 也不携带“上表面中心”之类的物理语义，其数值由网格生成和唯一编号顺序
决定。不能假定某个编号范围恒定对应某一表面。

本方案使用“从细胞中心看，节点方向与所选表面方向最一致”来定义表面中心。它能用
同一个算法覆盖六个方向，并能正确识别双凹 RBC 的中央凹陷。

## 6. 六个表面中心的统一选取算法

### 6.1 初始化几何

在 fresh run 中，加载两个 RBC 后、开始预松弛前执行一次选点。先同步粒子 envelope，
再收集目标细胞 owner 节点的：

```text
(baseCellId, vertexId, initialPosition)
```

MPI 下每个物理节点只由其 `localDomain` owner 提交一次。所有 rank 汇总后必须得到目标
cell type 的完整 `numVertex` 个不同 `vertexId`，否则立即报错。

目标细胞中心第一版取所有膜节点初始位置的算术平均：

```text
c = sum(x_i) / N
```

### 6.2 表面中心种子

把 `<face>` 转成单位向量 `n_face`。对每个目标细胞节点定义：

```text
r_i       = x_i - c
alignment = dot(r_i, n_face) / |r_i|
```

只考虑 `dot(r_i,n_face)>0` 的节点，在其中选择 `alignment` 最大者作为表面中心种子：

```text
seed = argmax alignment
```

相同 `alignment` 在浮点容差内并列时，用较小 `vertexId` 决定，保证不同 MPI 分解和不同
标准库排序实现得到相同结果。

这一规则对 `zPositive` 的意义是选择从细胞中心看正好沿 `+z` 的上层中央凹陷节点，
而不是选择离中心最远的隆起环。对 `xPositive`、`xNegative`、`yPositive`、`yNegative`
则自然落在相应的盘缘中心。

当前 `RBC_FROM_SPHERE`、`minNumTriangles=600`、`90 0 0` 姿态的启动回归期望为：

| 表面 | 期望中心 seed `vertexId` |
| --- | ---: |
| `xPositive` | `607` |
| `xNegative` | `600` |
| `yPositive` | `627` |
| `yNegative` | `620` |
| `zPositive` | `610` |
| `zNegative` | `615` |

这些编号只用于当前网格的回归检查，不能替代运行时几何选择。改变网格分辨率、构造方式
或初始姿态后，允许编号变化，但几何方向必须保持正确。

### 6.3 从 seed 扩展连续膜面补丁

不能简单选择 seed 周围三维欧氏距离小于某阈值的节点。双凹中心的上下两层膜空间距离
很小，欧氏球可能错误地同时选中上、下膜面。

本方案从 `HemoCellField::triangle_list` 建立无向膜网格图：

- 每个 `vertexId` 是一个图节点；
- 每个三角形的三条边是图边；
- 图边权重是两个端点的初始欧氏边长；
- 从 seed 运行 Dijkstra；
- 选择膜面测地距离不超过 `patchGeodesicRadius` 的节点。

```text
grip = { i | geodesicDistance(seed,i) <= patchRadiusLbm }
```

其中：

```text
patchRadiusLbm = patchGeodesicRadius_um * 1e-6 / param::dx
```

测地搜索沿膜拓扑扩展，因此小补丁不会从上层中央凹陷穿过空间直接跳到下层中央凹陷。

调试阶段可以额外支持 `patchMode=oneRing`：选择 seed 及其所有直接邻居。当前网格中
`zPositive` 的 one-ring 回归期望为：

```text
610, 264, 265, 274, 275, 402, 404
```

正式 benchmark 固定使用物理测地半径，避免固定节点数随网格加密改变实际加载面积。

## 7. 补丁节点的力分配

### 7.1 面积权重

对每个膜节点计算初始节点面积：

```text
A_i = (1/3) * sum(area_of_incident_triangle)
```

只在选中的补丁内归一化：

```text
w_i = A_i / sum(A_j for j in grip)
```

必须验证：

```text
sum(w_i) = 1
```

每个补丁节点在本步受到：

```text
f_i_lbm = loadScale(iter) * w_i * forceVector_lbm
```

所以目标细胞所受光镊合力严格为用户配置的矢量：

```text
sum(f_i_lbm) = loadScale(iter) * forceVector_lbm
```

`uniform` 模式仅用于与 `stretchCell` 类似的均匀节点分配回归：

```text
w_i = 1 / numberOfGripVertices
```

### 7.2 固定材料补丁

`vertexId` 和权重只在 fresh run 初始化时确定一次。预松弛、加载和形变过程中不根据
当前位置重新选择表面节点，也不重新计算权重。这样模拟的是光镊或粘附微球抓住一块
固定膜材料，而不是载荷区域在膜面上滑动。

如果未来要模拟持续作用在“当前空间最外侧”的压力或牵引，应作为另一种 follower-load
模型单独设计，不能混入本 benchmark。

## 8. helper 接口方案

建议新增以下公共类型：

```cpp
enum class CellSurfaceFace {
    XPositive,
    XNegative,
    YPositive,
    YNegative,
    ZPositive,
    ZNegative
};

enum class SurfacePatchMode {
    GeodesicRadius,
    OneRing
};

enum class SurfaceForceWeighting {
    NodalArea,
    Uniform
};

struct SurfaceForceVertex {
    plint vertexId;
    T weight;
    hemo::Array<T,3> initialPosition;
    T geodesicDistance;
    T nodalArea;
};
```

helper 类建议为：

```cpp
class HemoCellSurfaceForce {
public:
    HemoCellSurfaceForce(
        HemoCell &hemocell,
        const std::string &cellType,
        plint targetBaseCellId,
        CellSurfaceFace face,
        const hemo::Array<T,3> &forcePicoNewton,
        SurfacePatchMode patchMode,
        T patchGeodesicRadiusMicrometer,
        SurfaceForceWeighting weighting);

    void initializeFromCurrentGeometry();
    void restoreGripState(const std::string &filename);
    void saveGripState(const std::string &filename) const;
    void applyForce(T loadScale);

    const std::vector<SurfaceForceVertex> &selectedVertices() const;
    hemo::Array<T,3> forcePicoNewton() const;
    hemo::Array<T,3> forceLbm() const;
};
```

设计约束：

1. 不使用 `static lower_lsps/upper_lsps` 一类全局状态；每个 helper 实例拥有自己的目标
   cell、补丁和力矢量。
2. `targetBaseCellId` 通过 `HemoCellFields::base_cell_id(...)` 匹配，兼容周期复制 cell ID。
3. 构造函数只保存并校验配置；fresh run 的选点由显式
   `initializeFromCurrentGeometry()` 完成。
4. `applyForce()` 只向选中目标节点的 `sv.force` 做加法。
5. 不写入 `sv.force_repulsion`，否则会在 `iterate()` 开头被 adhesion 清零。
6. material timescale 和 particle velocity update timescale 第一版都必须为 `1`；否则直接
   拒绝运行，避免外力在 `sv.force` 中跨步累积。
7. `applyForce()` 只修改 owner 节点，随后调用 `syncEnvelopes()` 将新力同步给相邻 block
   的粒子副本。
8. `loadScale` 必须有限且位于 `[0,1]`；不在 helper 内维护迭代计数或加载阶段。

## 9. main 函数调度方案

新 main 复用 `twoCellShear` 已有的域、静止上下壁、RBC 注册、cell-cell adhesion、
cell-wall adhesion、checkpoint 和输出写法，但不设置上壁运动速度。

初始化顺序固定为：

```text
读取 config
→ 计算 LBM 参数和建 lattice
→ 建立静止上下壁
→ initializeCellfield()
→ 注册 RBC
→ 设置 material/particle/adhesion timescale=1
→ setAdhesion(...)
→ setBoundaryAdhesion(...)
→ 配置输出
→ fresh: loadParticles()
→ fresh: initializeFromCurrentGeometry()
→ fresh: 保存 grip_vertices.csv 和 grip_state.dat
→ restart: loadCheckPoint()
→ restart: restoreGripState(...)
→ 进入主循环
```

主循环为：

```cpp
while (hemocell.iter < tmax) {
    T scale = 0.0;

    if (hemocell.iter >= tRelax) {
        if (tForceRamp == 0) {
            scale = 1.0;
        } else {
            scale = std::min<T>(
                1.0,
                static_cast<T>(hemocell.iter - tRelax) /
                static_cast<T>(tForceRamp));
        }
    }

    if (scale > 0.0) {
        surfaceForce.applyForce(scale);
    }

    hemocell.iterate();

    // measurement/output/checkpoint
}
```

`applyForce()` 必须位于 `hemocell.iterate()` 之前。正常迭代开始时 adhesion 只清零并重算
`sv.force_repulsion`，随后 `spreadParticleForce()` 把：

```text
sv.force + sv.force_repulsion
```

共同铺展到流体。迭代末尾膜本构重算 `sv.force`，因此下一步必须重新施加一次光镊力。

`writeOutput()` 不应再次调用 `applyForce()`。现有 `OUTPUT_FORCE` 在迭代末尾和输出重算后
通常不包含刚刚已经铺展的瞬时光镊分量，因此光镊矢量和补丁合力由专用 CSV 明确记录。

## 10. checkpoint 策略

选中的材料补丁不是现有 HemoCell checkpoint 序列化字段。不能在变形后的 checkpoint
几何上重新寻找当前表面中心，否则 restart 前后可能抓住不同的膜节点。

第一版使用轻量 sidecar：

```text
checkpoint/grip_state.dat
```

其中保存：

```text
targetBaseCellId
face enum
force vector in pN
patch mode and radius
number of selected vertices
每个 vertexId、weight、initialPosition、geodesicDistance、nodalArea
```

保存 HemoCell checkpoint 时同步覆盖该文件；restart 时在第一次 `iterate()` 前读取并校验：

- target cell 存在；
- 所有 `vertexId` 都在目标 cell type 范围内；
- 所有权重有限且非负；
- 权重和在容差内等于 `1`；
- sidecar 的 face、force 和 patch 配置与当前 XML 一致；
- 不允许缺失 sidecar 时静默重新选点。

## 11. 输出与测量

RBC 输出至少包含：

```cpp
OUTPUT_POSITION
OUTPUT_TRIANGLES
OUTPUT_VELOCITY
OUTPUT_FORCE
OUTPUT_FORCE_REPULSION
OUTPUT_VERTEX_ID
OUTPUT_CELL_ID
```

`grip_vertices.csv` 在 fresh run 初始化时写一次：

```text
target_cell_id,face,vertex_id,x_um,y_um,z_um,
alignment,geodesic_distance_um,nodal_area_um2,weight
```

`twoCellPull.csv` 每次测量至少记录：

```text
iteration
stage                         # relaxation/ramp/constant-force
target_cell_id
face
load_scale
configured_force_x/y/z_pN
applied_force_x/y/z_pN
selected_vertex_count
sum_weights
weighted_grip_center_x/y/z_um
target_cell_center_x/y/z_um
other_cell_center_x/y/z_um
relative_center_x/y/z_um
minimum_cell_cell_distance_um
minimum_lower_cell_wall_distance_um
两个 RBC 的 area、volume 和 bbox
```

建议额外记录光镊补丁相对目标细胞中心的力矩：

```text
tau = sum((x_i-cellCenter) cross f_i)
```

法向中心拉脱 benchmark 中，非加载轴方向的合力和附加力矩应很小。面积权重不保证离散
网格下力矩严格为零，因此要记录而不是假定。

第一版不自动判定或提前终止“拉脱”。后处理可把以下条件作为候选判据：连续多个测量
周期内没有任何跨细胞节点对位于 `cellCellAdhesion/rc` 内，同时两细胞中心沿加载方向
持续分离。自动停止条件另立后续任务。

## 12. 参数校验与失败策略

启动时必须检查：

- `targetCellId>=0` 且加载粒子中确实存在该 base cell ID；
- `<face>` 严格属于六个允许值；
- 三个力分量均为有限数；
- `param::df>0`；
- `patchMode` 和 `weighting` 属于允许值；
- `geodesicRadius` 模式下半径严格大于零；
- `tRelax<=tmax`，`tRelax+tForceRamp<=tmax`；
- material、particle velocity 和 adhesion timescale 都为 `1`；
- 目标细胞所有顶点和三角形拓扑完整；
- 选中节点数大于零且小于目标细胞总节点数；
- 补丁没有包含当前回归已知的相反表面中心 seed；
- 力换算、节点权重、位置和面积均为有限数。

所有配置错误都应输出明确消息并退出，不采用默认 target cell、默认 face、跨字段回退或
硬编码拉力。

## 13. 源码级与算例级验证

### 13.1 选点验证

1. 当前 642 节点 RBC 对六个 face 分别返回本文件第 6.2 节列出的 seed ID。
2. `zPositive` 返回中央凹陷上层 seed `610`，不能返回最大 `z` 的隆起环节点。
3. `zNegative` 返回 `615`；`zPositive` 小补丁不得包含 `615`。
4. `oneRing + zPositive` 返回 seed 和六个已知邻居，共七个不同 vertexId。
5. 测地距离为零的节点只有 seed；所有被选节点都能由 triangle graph 从 seed 连通。
6. 改变细胞平移位置不改变选中的 vertexId。
7. 对受支持的初始整体旋转，使用相应实验室 face 后几何方向正确，不按旧硬编码 ID 选点。

### 13.2 力验证

1. `sum(weights)=1`。
2. owner 节点上的外力合计严格等于 `scale*forceVectorLbm`。
3. 未选节点和非目标 `cellId` 不得到光镊力。
4. 分别测试六个 face，并用与 face 平行的单位测试力确认符号。
5. 测试 face 与 force 不平行的切向加载，确认 helper 不擅自旋转力矢量。
6. 连续时间步外力不在 `sv.force` 中累计。
7. `writeOutput()` 不导致下一步外力重复或丢失。
8. 零拉力模式与无 helper 基准轨迹一致到数值容差。

### 13.3 MPI 与 checkpoint

1. 1 rank 和至少 2 ranks 得到完全相同的 seed、补丁 vertexId、权重和总外力。
2. 目标补丁跨 block/envelope 时没有漏力或重复力。
3. checkpoint 前后载荷阶段、vertexId、权重和轨迹连续。
4. 删除或篡改 `grip_state.dat` 时 restart 明确失败，不能从变形几何静默重选。

### 13.4 物理回归

1. relaxation 阶段光镊力严格为零，两类 adhesion 正常工作。
2. ramp 阶段合力线性增至配置值。
3. constant-force 阶段补丁合力保持配置值。
4. 下方 RBC 不是运动学固定，允许有限平移、转动和变形，但应由 wall adhesion 保持在壁面
   附近。
5. 上方 RBC 沿外力方向产生合理位移或拉脱，且没有因误选隆起环或相反膜面产生明显非物理
   折叠。
6. 上壁必须足够远，保证拉脱过程中目标 RBC 不进入上壁 cell-wall adhesion cutoff；若进入，
   该次结果不能作为纯 cell-cell 拉脱 benchmark。

## 14. 推荐实施顺序

1. 新建 `hemoCellSurfaceForce.h/.cpp`，先完成 face 枚举、MPI owner 节点汇总和 seed 选择。
2. 建立 triangle graph，实现 one-ring 和 Dijkstra 测地补丁。
3. 实现节点面积、权重归一化和 `grip_vertices.csv`，完成六方向选点测试。
4. 实现 owner 节点 `sv.force` 累加和 envelope 同步，完成总力、非目标节点和不累积测试。
5. 新建 `twoCellPull.cpp` 和 `config.xml`，复用双细胞与两类 adhesion 初始化。
6. 接入 relaxation、force ramp、专用 CSV 和 checkpoint sidecar。
7. 完成单进程、MPI、checkpoint、零载荷和六方向短程回归。
8. 最后进行力、补丁半径、`dt` 和粘附参数扫描，再确定正式 benchmark 参数。

## 15. 第一版明确不做的内容

- 不按每一步的当前空间坐标动态重新选择补丁；
- 不实现 follower force 或随膜法向旋转的力；
- 不把光镊力写入 `force_repulsion`；
- 不新增粒子序列化字段或修改现有 HDF5 数据布局；
- 不增加自动拉脱停止器、PID 力控制或位移控制；
- 不模拟有限尺寸光镊微球与膜的显式接触；
- 不修改现有 `HemoCellStretch` 行为；
- 不改变 cell-cell/cell-wall adhesion 公式或参数单位；
- 不在其他 HemoCell examples 中接入该 helper。

## 16. 当前配置基线与后续扫描参数

第一版 `config.xml` 采用当前 `twoCellShear/config.xml` 的盒子、流体、MPI 分块、RBC、
cell-cell adhesion 和 cell-wall adhesion 参数作为可运行基线。表面力默认目标为上方
`cellId=1`，作用面为 `zPositive`，总力为 `(0,0,50) pN`，测地半径为 `0.70 um`，
面积加权；`tForceRamp=0`，即预松弛结束后直接按配置总力加载。这些都不是 helper 中的
硬编码值，正式 benchmark 仍应在 XML 中确定或扫描：

1. 默认目标 `cellId`；当前双细胞布局建议为上方 `cellId=1`。
2. 第一组 benchmark 的 face；法向上拉建议为 `zPositive`。
3. `forceX/Y/Z` 的扫描范围和步长，单位 pN。
4. 光镊加载补丁的物理测地半径；`0.70 um` 只是与当前一圈邻居尺度接近的起始建议。
5. `tRelax`、`tForceRamp` 和恒力观察时间。
6. 域高度，必须覆盖最大预期拉脱位移并使目标 RBC 远离上壁 cutoff。
7. cell-cell 和 cell-wall adhesion 的正式参数及拉脱判据。

上述数值全部由配置读取，不在 helper 或 main 中设置物理调参默认值。

## 17. 第一版实施验证记录

第一版已完成以下构建和短程回归：

- 默认 HemoCell 静态库和独立 `twoCellPull` 可执行文件编译通过；
- 当前 642 节点 RBC 的六方向 seed 分别为 `607/600/627/620/610/615`，与
  第 6.2 节的 `x+/x-/y+/y-/z+/z-` 回归值一致；
- `geodesicRadius=0.70 um` 选中 seed 及其六个邻居，共 7 个节点；
- `nodalArea` 权重和为 1，专用 CSV 记录的加载步总力为 `(0,0,50) pN`；
- 1 rank 与 `1×1×2` 的 2 ranks 运行产生完全相同的 `grip_vertices.csv` 和
  `grip_state.dat`；
- checkpoint 能恢复 seed 610 和同一组 7 个材料节点，不会在变形后的膜面重新选点。
