# AMReX 通用三维局部加密后端修改说明

本文对应分支 `C++20_AMReX`，修改日期为 2026-11-09。

> 当前工作树已经清理了 `demo/` 根目录下的旧算例输入。本文中提到
> `demo/mpm-amr.json`、`demo/amr-layout.example.json` 等文件的内容属于开发过程记录；
> 当前可运行的标准算例是 `demo/amr_two_level_sine/`。

## 1. 本次实现范围

本次实现的是通用的静态块结构三维 AMR backend：

- 任意三维物理域；
- 任意整数粗网格数量；
- 任意整数 refinement ratio；
- 任意数量 refinement level；
- 每级可包含多个互不重叠的矩形 Box；
- 生成 AMReX `BoxArray`；
- 生成 AMReX `DistributionMapping`；
- 将物理坐标映射到指定 level；
- 查找某点所在的 finest valid level；
- 判断粗单元是否被细层覆盖；
- 支持 MPI patch 分区；
- 将 finest valid AMR cells 转换为 MPM Hex8 Cell；
- 按物理坐标去重生成 MPM NodeBase；
- 建立 `cell_id -> AmrCellKey` 映射；
- 接入 `MPMBase::initialise_mesh()` 的 `io_type="Amr"` 分支；
- 粒子通过 AMReX finest-level 查询定位；
- 提供独立布局预览工具；
- 提供自动测试。

本次**没有**完成：

- MPM 粗细界面节点力交换；
- 动态 `TagBox/regrid`；
- AMReX `ParticleContainer` 粒子迁移；
- 粗细界面守恒修正。

因此当前状态是：

```text
通用 AMR 网格拓扑 backend：已实现
AMR 布局和 finest-cell 查询：已实现
AMR 与 MPM Mesh/粒子定位接入：已实现
AMR 粗细界面物理耦合：尚未实现
```

这样可以先审查网格层级、Box 划分、MPI 分区和局部加密规则，再进入 MPM 界面耦合阶段。

## 2. 新增文件

```text
include/mesh_adaptivity/amr_config.h
include/mesh_adaptivity/amr_hierarchy.h
include/mesh_adaptivity/amrex_amr_backend.h

src/mesh_adaptivity/amr_config.cc
src/mesh_adaptivity/amrex_amr_backend.cc
src/mesh_adaptivity/amrex_amr_layout_tool.cc
src/mesh_adaptivity/amr_mesh_adapter.cc

tests/mesh_adaptivity/amrex_amr_backend_test.cc
tests/mesh_adaptivity/amrex_amr_mesh_adapter_test.cc

demo/amr-layout.example.json
demo/mpm-amr.json
demo/amr-particles-3d.txt
```

## 3. `amr_config.h`

### 第 8 行：引入 JSON

```cpp
#include "json.hpp"
```

该头文件使用项目已有的 nlohmann/json，不引入新的 JSON 依赖。

### 第 13 行：定义 AMR JSON 类型别名

```cpp
using AmrJson = nlohmann::json;
```

用于避免依赖 `io.h` 中的全局 `Json` 别名。

### 第 15-18 行：物理矩形 Box

```cpp
struct AmrPhysicalBox {
  std::array<double, 3> minimum{};
  std::array<double, 3> maximum{};
};
```

物理 Box 使用真实坐标，不直接使用 AMReX `IntVect`。

### 第 20-23 行：某一 refinement level 的 Box 集合

```cpp
struct AmrLevelRegion {
  int level{0};
  std::vector<AmrPhysicalBox> boxes;
};
```

同一级允许多个互不重叠的矩形加密区。

### 第 25-32 行：AMR 总配置

```cpp
struct AmrConfig {
  std::array<double, 3> domain_minimum{};
  std::array<double, 3> domain_maximum{};
  std::array<int, 3> coarse_cells{};
  std::array<int, 3> refinement_ratio{2, 2, 2};
  int max_grid_size{32};
  std::vector<AmrLevelRegion> levels;
};
```

逐项说明：

- `domain_minimum/maximum`：完整物理域；
- `coarse_cells`：Level 0 三个方向的整数单元数；
- `refinement_ratio`：相邻级三个方向的整数加密比；
- `max_grid_size`：AMReX patch 的最大单元数；
- `levels`：Level 1 以上的加密区域。

### 第 35-36 行：配置解析接口

```cpp
bool parse_amr_config(const AmrJson& json, AmrConfig* config,
                      std::string* error = nullptr);
```

## 4. `amr_hierarchy.h`

### 第 9-12 行：AMR cell 唯一键

```cpp
struct AmrCellKey {
  int level{0};
  std::array<int, 3> index{};
};
```

后续 MPM `Cell` 必须能映射到 `(level, i, j, k)`。

### 第 14-20 行：patch 描述

```cpp
struct AmrPatchInfo {
  int level{0};
  int patch{0};
  int rank{0};
  std::array<int, 3> lo{};
  std::array<int, 3> hi{};
};
```

用于布局预览、调试和后续生成 MPM Cell。

## 5. `amrex_amr_backend.h`

### 第 9-10 行：AMReX 类型

```cpp
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
```

本 backend 使用 AMReX patch 拓扑和 MPI patch 分配。

### 第 18-19 行：backend 类

```cpp
class AmrexAmrBackend {
```

它只负责网格层级，不保存 MPM 材料或粒子历史。

### 第 21 行：初始化

```cpp
bool initialise(const AmrConfig& config, std::string* error = nullptr);
```

### 第 23-27 行：读取 AMReX 层级数据

```cpp
int number_of_levels() const;
const amrex::BoxArray& box_array(int level) const;
const amrex::DistributionMapping& distribution_mapping(int level) const;
amrex::Box domain_box(int level) const;
std::array<double, 3> cell_size(int level) const;
```

这些接口供后续 MPM Mesh 适配器使用。

### 第 29-32 行：点定位和 finest-cell 判断

```cpp
std::optional<AmrCellKey> locate_finest(point) const;
bool is_finest_valid_cell(level, index) const;
```

`locate_finest` 从最高级向低搜索。

`is_finest_valid_cell` 用于生成粒子或 Cell 时排除被细层覆盖的粗单元。

### 第 34 行：patch 枚举

```cpp
std::vector<AmrPatchInfo> patches(int level) const;
```

### 第 40-45 行：内部层级存储

```cpp
AmrConfig config_;
std::vector<amrex::Box> level_domains_;
std::vector<amrex::BoxArray> level_boxes_;
std::vector<amrex::DistributionMapping> level_distributions_;
std::vector<std::array<double, 3>> level_cell_sizes_;
std::vector<std::array<int, 3>> level_cell_counts_;
```

## 6. `amr_config.cc`

### 第 16-31 行：浮点三元组解析

```cpp
bool parse_triplet(...)
```

要求：

- JSON 必须是三个数；
- 每个值必须是有限数。

### 第 33-48 行：整数三元组解析

```cpp
bool parse_int_triplet(...)
```

要求：

- JSON 必须是三个整数；
- 每个值必须大于零。

### 第 52-65 行：读取物理域和粗单元数

```cpp
json.at("domain_min")
json.at("domain_max")
json.at("coarse_cells")
```

### 第 67-77 行：读取 refinement ratio 和 patch 最大尺寸

```cpp
json.contains("refinement_ratio")
json.contains("max_grid_size")
```

### 第 79-86 行：检查物理域和 level 数组

确保 `domain_max > domain_min`，且 `levels` 非空。

### 第 88-117 行：逐级读取加密 Box

```cpp
region.level
region.boxes
```

同一 level 的 Box 允许有多个，但 backend 会检查它们是否重叠。

### 第 119-127 行：排序并检查 level 连续

要求 Level 1、2、3……连续出现。

## 7. `amrex_amr_backend.cc`

### 第 28-43 行：配置有效性检查

检查：

- 物理域边界；
- 粗单元数量；
- refinement ratio；
- 清空旧层级数据。

### 第 52-70 行：创建 Level 0

```cpp
coarse_boxes.maxSize(config_.max_grid_size);
```

Level 0 先用一个完整 `Box`，再由 AMReX 拆成多个 patch。

### 第 72-79 行：计算下一级网格尺寸

```cpp
child_cell_size[axis] /= refinement_ratio[axis];
child_cell_counts[axis] *= refinement_ratio[axis];
```

### 第 81-105 行：物理 Box 转 AMReX 索引 Box

```cpp
lower = (box_min - domain_min) / child_cell_size;
upper = (box_max - domain_min) / child_cell_size;
lo = floor(lower);
hi = ceil(upper) - 1;
```

这一步保证任意物理矩形区域都可以转成整数索引 Box。

### 第 107-112 行：检查同级 Box 重叠

同一级多个 patch 不允许相互覆盖。

### 第 114-120 行：检查 proper nesting

细层 Box 先 coarsen，再检查它是否被父层 `BoxArray` 覆盖。

### 第 122-128 行：创建 BoxArray 和 DistributionMapping

```cpp
amrex::BoxArray boxes(box_list);
boxes.maxSize(config_.max_grid_size);
level_distributions_.emplace_back(level_boxes_.back());
```

这一步真正创建 AMReX patch 和 MPI 分区。

### 第 130-136 行：保存该级 cell size 和 count

### 第 160-175 行：物理坐标转 level 索引

```cpp
relative = (point - domain_min) / level_dx;
index = floor(relative);
```

### 第 177-188 行：从最高级向下查找 finest cell

```cpp
for (int level = number_of_levels() - 1; level >= 0; --level)
```

找到第一个覆盖该点的 level 即返回。

### 第 190-204 行：排除被覆盖的粗单元

```cpp
center = cell_center;
finest = locate_finest(center);
return finest->level == level;
```

后续生成 MPM Cell 或粒子时应使用此判断，防止粗细层重复生成。

### 第 206-224 行：输出 patch 信息

返回每个 patch 的：

```text
level
patch id
MPI rank
lo index
hi index
```

## 8. `amrex_amr_layout_tool.cc`

### 第 7-18 行：读取配置文件

```text
amrex_amr_layout_tool AMR_CONFIG.json
```

### 第 20-27 行：解析 AMR JSON

### 第 29-33 行：初始化 backend

### 第 35-45 行：输出每级 patch 和 MPI rank

该工具用于在接入 MPM Cell 之前检查网格层级是否正确。

## 9. `amrex_amr_backend_test.cc`

### 第 10-18 行：构造三级测试配置

### 第 20-25 行：初始化 backend

### 第 27-36 行：检查 level 数和 BoxArray

### 第 38-52 行：检查 finest level 定位

### 第 54-57 行：检查计算域外点

### 第 59-66 行：检查粗单元被细层覆盖

### 第 68-80 行：检查每级 patch 和 MPI rank

## 10. `amr_mesh_adapter.h`

### 第 19-20 行：MPM Mesh 适配器

```cpp
class AmrMeshAdapter {
```

作用是把 AMReX 的 finest valid cells 转换为现有 MPM `Mesh<3>`。

### 第 22-25 行：填充 MPM Mesh

```cpp
bool populate(const AmrexAmrBackend& backend, Mesh<3>* mesh,
              const std::string& node_type,
              const std::shared_ptr<Element<3>>& element,
              std::string* error = nullptr);
```

### 第 27-29 行：查询 AMR cell 对应的 MPM cell ID

```cpp
const std::vector<AmrCellKey>& cell_keys() const;
std::optional<Index> cell_id(int level, const std::array<int, 3>& index) const;
```

## 11. `amr_mesh_adapter.cc`

### 第 23-26 行：检查目标 Mesh

只允许填充空的 MPM Mesh，避免与已有 mesh.txt 数据混用。

### 第 28-33 行：节点和单元临时容器

```cpp
node_coordinates
node_ids
cell_connectivity
cell_keys_
cell_ids_
```

### 第 35-43 行：按物理坐标去重节点

同一个物理位置的节点在不同 Cell 中只创建一个 MPM Node。

### 第 45-53 行：遍历 level、patch、i/j/k

```cpp
for (int level = 0; ...)
  for (const auto& patch : backend.patches(level))
    for k, j, i
```

只有 `is_finest_valid_cell()` 返回 true 的单元才会生成 MPM Cell。

### 第 55-66 行：计算 Hex8 八个角点

节点顺序必须是：

```text
0: x-,y-,z-
1: x+,y-,z-
2: x+,y+,z-
3: x-,y+,z-
4: x-,y-,z+
5: x+,y-,z+
6: x+,y+,z+
7: x-,y+,z+
```

### 第 68-71 行：保存 MPM cell ID 和 AMR key

```cpp
cell_keys_
cell_ids_
cell_connectivity
```

### 第 78-85 行：写入 MPM Mesh

```cpp
mesh->create_nodes(...)
mesh->create_cells(...)
```

## 12. `amrex_amr_mesh_adapter_test.cc`

### 第 10-18 行：构造三级 AMR 配置

### 第 20-25 行：初始化 AMReX backend

### 第 27-33 行：创建 Hex8 Element 并填充 MPM Mesh

### 第 35-40 行：检查 Cell 数量、Node 数量和 AMR key 映射

### 第 42-44 行：输出测试结果

## 13. `CMakeLists.txt`

### 第 27-46 行：GCC 13.3 ICE 兼容选项

新增 `MPM_GCC13_ICE_WORKAROUND`，只匹配 GCC 13.3，关闭会触发
`cfgcleanup` 内部错误的若干优化 pass。默认开启，可通过 CMake 关闭。

### 第 260-267 行：AMR backend 加入主 MPM target

```cmake
target_sources(mpm PRIVATE
  src/mesh_adaptivity/amr_config.cc
  src/mesh_adaptivity/amrex_amr_backend.cc
)
```

### 第 279-289 行：布局预览工具

### 第 291-301 行：backend 测试

### 第 303-316 行：AMR Mesh 适配器测试

所有 backend 代码只在 `MPM_WITH_AMREX=ON` 时构建。

## 14. `demo/amr-layout.example.json`

### 第 2-6 行：物理域和粗网格

```json
"domain_min": [0, 0, 0],
"domain_max": [8000, 8000, 8000],
"coarse_cells": [10, 10, 10],
"refinement_ratio": [2, 2, 2],
"max_grid_size": 4
```

### 第 7-26 行：两级局部加密

Level 1 和 Level 2 各自使用物理矩形 Box。

## 15. 构建和运行

```bash
cmake -S ~/mpm -B ~/mpm/build_amr_backend \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_PREFIX_PATH="$HOME/opt/amrex-26.11" \
  -DMPM_WITH_AMREX=ON \
  -DNO_KAHIP=True

cmake --build ~/mpm/build_amr_backend --parallel 8
```

运行测试：

```bash
~/mpm/build_amr_backend/amrex_amr_backend_test

mpirun --oversubscribe -np 2 \
  ~/mpm/build_amr_backend/amrex_amr_backend_test
```

预览布局：

```bash
~/mpm/build_amr_backend/amrex_amr_layout_tool \
  ~/mpm/demo/amr-layout.example.json
```

## 16. 已验证结果

- AMReX backend 编译成功。
- 主 MPM target 在 `MPM_WITH_AMREX=ON` 下编译成功。
- 单 MPI rank backend 测试返回 0。
- 双 MPI rank backend 测试返回 0。
- 示例 JSON 解析成功。
- 示例生成 3 个 level、每级 27 个 patch。
- finest-level 和 covered-coarse-cell 查询测试通过。
- AMR-to-MPM Mesh 适配器测试通过。
- 示例最终生成 2498 个 MPM Cell、3256 个 MPM Node 和 550 个悬点约束。
- 双 MPI rank 下每个 rank 生成相同全局网格：
  2498 个 Cell、3256 个 Node、550 个悬点约束。
- `io_type="Amr"` 已接入 `MPMBase::initialise_mesh()`。
- 5 个测试粒子均通过 AMReX finest-level 定位。
- 直径/边长不落在父单元边界上的细化盒会被拒绝。
- 跨两级的 2:1 非平衡界面会被拒绝。
- 悬点质量、动量、外力和内力投影的守恒测试通过。
- 悬点速度、加速度的 Hex8 形函数重建测试通过。
- AMR MPM 示例完成一步计算，返回 0。
- 双 MPI rank AMR MPM 示例完成一步计算，返回 0。

## 17. Mesh 接入

### `include/mesh/mesh.h` 第 213-241 行

新增：

```cpp
void assign_amr_backend(...);
void assign_amr_cell_key(...);
std::optional<Index> amr_cell_id(...);
```

### `include/mesh/mesh.h` 第 938-949 行

新增：

```cpp
std::shared_ptr<AmrexAmrBackend> amr_backend_;
std::map<std::array<int, 4>, Index> amr_cell_ids_;
```

### `include/mesh/mesh.tcc` 第 307-417 行

实现 AMR backend 指针、`(level,i,j,k) -> cell_id` 映射的保存与查询，以及悬点约束的保存、质量/动量凝聚、节点力投影和运动学重建。

### `include/mesh/mesh.tcc` 第 1219-1243 行

粒子定位增加 AMR 快速路径：

```text
particle coordinates
  -> AmrexAmrBackend::locate_finest
  -> find MPM cell id
  -> Cell::is_point_in_cell
  -> assign particle cell
```

只有在 `MPM_WITH_AMREX` 和 `Tdim == 3` 时启用；均匀网格路径保持不变。

## 18. MPMBase 接入

### `include/solvers/mpm_base.tcc` 第 254-289 行

新增 `io_type == "Amr"` 分支：

1. 读取 `node_type` 和 `cell_type`；
2. 解析 `mesh.amr`；
3. 创建 `AmrexAmrBackend`；
4. 绑定 backend 到 MPM Mesh；
5. 用 `AmrMeshAdapter` 生成 Cell 和 Node；
6. 计算邻居；
7. 跳过文件型 mesh 读取和实体集合读取；
8. 返回主 MPM 初始化流程。

### `include/solvers/mpm_base.tcc` 第 451-471 行

AMR 模式不创建不存在的 `IOMesh("Amr")`：

- `Ascii3D` 仍读取 particles-cells；
- `Amr` 直接调用 `Mesh::locate_particles_mesh()`；
- AMR 模式不写回 particles-cells 文件。

### `include/solvers/mpm_base.tcc` 第 881-886 行

加载条件中的文件 reader 改为惰性创建。AMR 模式无文件型载荷时不再查找
`IOMesh("Amr")`。

### `include/solvers/mpm_base.tcc` 第 982-990 行

数学函数只有存在 `file` 时才创建 reader；AMR 模式遇到文件型数学函数会给出
明确错误。

## 19. 最小 AMR MPM 示例

新增：

```text
demo/mpm-amr.json
demo/amr-particles-3d.txt
```

### `demo/mpm-amr.json` 第 4-39 行

定义 `io_type="Amr"`、Mesh 类型和嵌套的 `mesh.amr` 配置。

### `demo/mpm-amr.json` 第 41-53 行

使用现有的 `Ascii3D` 粒子文件生成器。

### `demo/mpm-amr.json` 第 54-65 行

定义线弹性岩石材料和零重力。

### `demo/mpm-amr.json` 第 67-94 行

定义单步 `MPMExplicit3D/USF` 分析。

### `demo/amr-particles-3d.txt` 第 1-6 行

包含 5 个测试粒子的数量和坐标。

运行时：

```bash
~/mpm/build_amr_integration_normal/mpm \
  -f ~/mpm/demo/ -p 2 -i mpm-amr.json
```

该示例使用：

```text
3 个 level
2498 个 MPM Cell
3256 个 MPM Node
550 个悬点约束
5 个粒子
1 个时间步
```

## 20. 当前边界

本次已经完成参考论文 3.1 节的直接悬点凝聚，但当前 AMR 模式仍是静态局部加密后端，不应被描述为完整的生产级动态 AMR-MPM。现阶段明确支持：

1. 三维、正交六面体 `Hex8`。
2. 静态层级的 finest valid cell 到 MPM mesh 转换。
3. 细化盒边界必须与父层单元边界对齐。
4. 相邻层差为 1，并检查 2:1 平衡。
5. 悬点质量、动量和节点力的直接凝聚。
6. 粗节点推进后重建悬点速度和加速度。

当前仍不支持：

1. 动态 `TagBox/regrid` 和每步变化拓扑。
2. AMReX `ParticleContainer` 粒子迁移。
3. 在悬点上直接施加独立的位移、速度或加速度约束；应把约束施加到粗支撑节点。
4. 文件型实体集合、边界条件、节点集中力到 AMR 层级的完整映射。
5. AMR 专用的高阶 VTK 输出和压力平滑界面修正。
6. 隐式/Newmark AMR 粗细分界面凝聚。

## 21. 既有 AMR backend 审查重点

建议按以下顺序审查：

1. `amr_config.h` 的输入数据是否足够。
2. `amr_config.cc` 的配置检查是否合理。
3. `amrex_amr_backend.cc` 第 72-134 行的物理坐标到 index Box 转换。
4. proper nesting 检查是否满足预期的 local refinement。
5. `locate_finest()` 和 `is_finest_valid_cell()` 是否覆盖边界点。
6. patch 是否由 MPI 正确分配。
7. 下一步应按第 22 节逐项复核悬点凝聚的代码位置。

## 22. 粗细界面悬点凝聚：逐行修改说明

### 22.1 论文 3.1 节采用的约束

设悬点为 \(h\)，它位于粗单元 \(K_h\) 的面或边上，但不是粗单元顶点。论文把悬点视为非独立自由度：

\[
u_h = \sum_{i \in N(K_h)} \phi_i^{K_h}(x_h)u_i
\]

因此：

1. 不增加拉格朗日乘子。
2. 不增加罚参数。
3. 不让悬点保留独立质量、动量和节点力。
4. 先把悬点贡献投影到 8 个粗单元支撑节点，再由粗节点推进。
5. 粗节点求速度、加速度后，按同一组形函数权重重建悬点值。

本次使用规则正交 Hex8 的三线性形函数：

\[
N_i(\xi,\eta,\zeta)=
\frac{1}{8}(1+\xi_i\xi)(1+\eta_i\eta)(1+\zeta_i\zeta)
\]

其中 \(\xi,\eta,\zeta\in[-1,1]\)。一个悬点最多得到 8 个非零权重；位于面心时通常只有 4 个非零权重，位于边中点时通常只有 2 个非零权重。

当前 MPM 显式求解器使用对角集中质量，因此质量凝聚采用一致质量矩阵
\(m_h\mathbf w\mathbf w^T\) 的行和：

\[
m_i \mathrel{+}= w_i m_h
\]

这不是完整的一致质量矩阵耦合，而是与现有显式 MPM 对角质量框架相容的集中质量近似。若以后改为隐式/一致质量矩阵，应把这里替换为矩阵形式的静力凝聚。

### 22.2 `include/mesh_adaptivity/amr_hierarchy.h` 第 24-30 行

新增 `HangingNodeConstraint`：

```cpp
struct HangingNodeConstraint {
  Index hanging_node{0};
  AmrCellKey parent_cell{};
  std::array<Index, 8> support_nodes{};
  std::array<double, 8> weights{};
};
```

逐行含义：

1. `hanging_node`：MPM 网格中的悬点 ID。
2. `parent_cell`：该悬点选定的粗父单元 `(level,i,j,k)`。
3. `support_nodes`：粗父单元 8 个顶点对应的 MPM Node ID。
4. `weights`：悬点在粗父单元中的 Hex8 形函数值。

### 22.3 `include/mesh_adaptivity/amr_mesh_adapter.h` 第 28-31、47-49 行

新增 `hanging_node_constraints()` 访问器和内部成员
`hanging_node_constraints_`。这样适配器既能生成约束，也能让测试在不进入求解器时检查约束总数。

### 22.4 `src/mesh_adaptivity/amr_mesh_adapter.cc` 第 18-25 行

新增坐标归一化与近似比较：

```cpp
double rounded_coordinate(double value);
bool near(double lhs, double rhs);
```

原因：细层节点和粗单元顶点可能通过不同的浮点运算路径生成，直接比较
`double` 会出现“数学上同一点、二进制上不相等”的问题。节点键统一四舍五入到
`1e-10`，实际写入 MPM mesh 的坐标仍保留原始精度。

### 22.5 `src/mesh_adaptivity/amr_mesh_adapter.cc` 第 43-49 行

每次执行 `populate()` 先清空旧的约束，并用归一化坐标建立“物理坐标到 MPM Node ID”的映射。

### 22.6 `src/mesh_adaptivity/amr_mesh_adapter.cc` 第 92-96 行

建立 `node -> incident cell levels` 表：

```cpp
std::vector<std::set<int>> incident_levels(node_coordinates.size());
```

对每个 finest valid cell 的 8 个角点，记录该点接触到的 AMR level。

### 22.7 `src/mesh_adaptivity/amr_mesh_adapter.cc` 第 99-108 行

对每个节点检查层级连续性。若同一节点同时接触 `L` 和 `L+2`，却没有任何
`L+1` 单元，则直接报错：

```text
AMR transition violates the 2:1 level-balance rule
```

这一步防止跨两级直接连接造成不唯一的悬点约束。

### 22.8 `src/mesh_adaptivity/amr_mesh_adapter.cc` 第 111-121 行

把当前节点接触的最高层记为 `finest_level`，父层取：

```cpp
const int parent_level = finest_level - 1;
```

然后在父层索引空间内检查 `base-1/base/base+1` 三个候选索引，处理节点正好落在粗单元边界上的情况。

### 22.9 `src/mesh_adaptivity/amr_mesh_adapter.cc` 第 123-165 行

遍历候选粗单元并检查：

1. 该粗单元必须是 finest valid cell，也就是没有被更细层完全覆盖。
2. 节点必须位于粗单元的闭区间内。
3. 节点必须落在 1 个或 2 个方向上，即粗单元的面或边上。
4. 若节点同时是粗单元的 3 方向顶点，则它不是悬点，跳过。
5. 若一个悬点同时接触多个粗单元，选 MPM `cell_id` 最大的父单元。

最后一条使约束唯一。按论文的直接约束思想，悬点只允许一个父单元形函数作为其从属关系；边和角附近相邻粗单元在共享边上的形函数限制相同，因此选择最大 `cell_id` 不会引入方向偏置。

### 22.10 `src/mesh_adaptivity/amr_mesh_adapter.cc` 第 169-191 行

按 MPM Hex8 节点顺序列出粗父单元 8 个顶点：

```text
0: (xmin,ymin,zmin)
1: (xmax,ymin,zmin)
2: (xmax,ymax,zmin)
3: (xmin,ymax,zmin)
4: (xmin,ymin,zmax)
5: (xmax,ymin,zmax)
6: (xmax,ymax,zmax)
7: (xmin,ymax,zmax)
```

通过归一化后的物理坐标查出每个顶点已经存在的 MPM Node ID。若粗父单元顶点在
finest-cell 节点表中不存在，说明层级或网格构造不一致，直接报错。

### 22.11 `src/mesh_adaptivity/amr_mesh_adapter.cc` 第 193-213 行

把悬点物理坐标转换到粗父单元的自然坐标：

```text
xi   = 2*(x-xmin)/(xmax-xmin) - 1
eta  = 2*(y-ymin)/(ymax-ymin) - 1
zeta = 2*(z-zmin)/(zmax-zmin) - 1
```

再计算 8 个三线性形函数值。该数组就是约束：

```text
u_h = sum(weights[i] * u_support_nodes[i])
```

### 22.12 `src/mesh_adaptivity/amr_mesh_adapter.cc` 第 214-221 行

检查权重和是否为 1，并保存约束。权重和检查保证常数速度场、常数加速度场能够被悬点精确重现。

### 22.13 `src/mesh_adaptivity/amr_mesh_adapter.cc` 第 224-234 行

禁止“悬点再作为另一个悬点支撑节点”的链式约束。当前实现只接受相邻层差 1 的单一约束层，避免两级以上插值链造成重复凝聚。

### 22.14 `src/mesh_adaptivity/amr_mesh_adapter.cc` 第 236-242 行

MPM Node 和 Cell 创建完成后，把约束表整体交给 `Mesh`：

```cpp
mesh->assign_hanging_node_constraints(hanging_node_constraints_);
```

### 22.15 `src/mesh_adaptivity/amrex_amr_backend.cc` 第 107-116 行

新增细化盒对齐检查：

```cpp
raw_box.smallEnd(axis) % ratio == 0
(raw_box.bigEnd(axis) + 1) % ratio == 0
```

原因：有限元悬点约束要求粗细分界面落在整粗单元边界上。若细层边界从粗单元内部穿过，悬点所在位置根本不是粗单元的面或边，无法构造论文 3.1 节的单一粗父单元约束。

### 22.16 `include/nodes/node_base.h` 第 193-201 行

新增纯虚接口：

```cpp
virtual void assign_velocity_acceleration(
    unsigned phase, const VectorDim& velocity,
    const VectorDim& acceleration) noexcept = 0;
```

它用于求解完成后一次性把插值得到的速度和加速度写回悬点。

### 22.17 `include/nodes/node.h` 第 199-203 行及 `include/nodes/node.tcc` 第 274-285 行

实现该接口。写入时加 `node_mutex_`，使行为与 Node 现有质量、动量、节点力写入保持一致。

### 22.18 `include/mesh/mesh.h` 第 226-241、938-949 行

新增 5 个 Mesh 接口和 1 个成员：

```text
assign_hanging_node_constraints
hanging_node_constraints
condense_hanging_node_mass_momentum
project_hanging_node_forces
reconstruct_hanging_node_kinematics
hanging_node_constraints_
```

### 22.19 `include/mesh/mesh.tcc` 第 334-347 行

保存并返回悬点约束表。非 AMR 构建不会编译这些成员，也不会改变均匀网格路径。

### 22.20 `include/mesh/mesh.tcc` 第 349-369 行

`condense_hanging_node_mass_momentum()` 逐约束执行：

```text
读 m_h、p_h
-> 悬点质量、动量清零
-> 对 8 个支撑节点累加 weight[i] * m_h
-> 对 8 个支撑节点累加 weight[i] * p_h
```

该操作保持：

\[
m_h=\sum_i \tilde m_i,\qquad
\mathbf p_h=\sum_i \tilde{\mathbf p}_i
\]

### 22.21 `include/mesh/mesh.tcc` 第 371-391 行

`project_hanging_node_forces()` 对外力与内力分别执行同一线性投影：

\[
\tilde f_i^{ext/grav/traction}\mathrel{+}=w_i f_h^{ext},\qquad
\tilde f_i^{int}\mathrel{+}=w_i f_h^{int}
\]

随后悬点自身外力和内力清零。由于投影发生在所有节点力装载完成之后，重力、体力、面力、接触力和集中节点力都能进入同一投影。

### 22.22 `include/mesh/mesh.tcc` 第 393-409 行

`reconstruct_hanging_node_kinematics()` 按：

\[
\mathbf v_h=\sum_i w_i\mathbf v_i,\qquad
\mathbf a_h=\sum_i w_i\mathbf a_i
\]

重建悬点速度和加速度，再通过
`Node::assign_velocity_acceleration()` 写回。

### 22.23 `include/solvers/mpm_scheme/mpm_scheme.tcc` 第 69-72 行

P2G 质量和动量映射、MPI halo exchange 完成后，在节点速度计算之前执行质量/动量凝聚。

顺序必须是：

```text
Particle -> hanging/support nodes
-> MPI halo exchange
-> condense hanging mass/momentum
-> compute support velocity
```

### 22.24 `include/solvers/mpm_scheme/mpm_scheme.tcc` 第 249-253、268-272 行

节点力全部装载后、计算节点加速度前，投影悬点力；粗支撑节点完成加速度/速度推进后，重建悬点速度和加速度。顺序是：

```text
assemble all nodal forces
-> project hanging forces
-> advance coarse support nodes
-> reconstruct hanging velocity/acceleration
-> Particle G2P update
```

### 22.25 `include/solvers/mpm_base.tcc` 第 279-283 行

启动日志现在同时报告悬点数和约束数。示例运行结果为：

```text
3 levels, 2498 cells, 3256 nodes, 550 hanging nodes
```

### 22.26 `tests/mesh_adaptivity/amrex_amr_mesh_adapter_test.cc` 第 39-145 行

新增测试覆盖：

1. 适配器确实生成悬点约束。
2. 每条约束权重和为 1。
3. 悬点质量凝聚前后总质量不变。
4. 悬点动量凝聚前后总动量不变。
5. 外力和内力投影前后分别守恒。
6. 重建速度等于支撑节点速度的形函数加权和。
7. 重建加速度等于支撑节点加速度的形函数加权和。
8. 跨两级直接连接的配置会被 2:1 检查拒绝。

### 22.27 `tests/mesh_adaptivity/amrex_amr_backend_test.cc` 第 82-89 行

新增反例：把 Level 1 细化盒下边界从 `1600` 改成 `2000`，该边界不再落在
800 大小的粗单元边界上，backend 必须拒绝初始化。

### 22.28 示例配置

`demo/amr-layout.example.json` 第 12-22 行和 `demo/mpm-amr.json` 第 24-34 行：

```text
Level 1: [0,0,1600] 到 [4000,4000,6400]
Level 2: [0,0,3200] 到 [1600,1600,4800]
Refinement ratio: 2
```

所有边界都能被父层单元尺寸整除，满足有限元粗细界面对齐要求。

## 23. 后续在程序里使用 AMReX 时必须遵守

1. 新建 AMR 算例时，先根据粗单元尺寸计算每个细化盒边界；边界必须是父单元边界坐标。
2. Level `L+1` 的细化区不能跳过 Level `L` 直接贴到 Level `L-1` 单元。
3. 需要加载荷或位移约束时，优先选粗支撑节点。悬点 ID 可以查看，但不应作为独立边界条件自由度。
4. AMR 模式目前只用于三维显式 MPM。切换为隐式/Newmark 前，需要为矩阵装配增加对应的约束凝聚和静力凝聚。
5. 每次修改 AMR layout 后，先运行 `amrex_amr_backend_test` 和 `amrex_amr_mesh_adapter_test`，再运行一两个时间步的 MPM 算例。

推荐验证命令：

```bash
cd ~/mpm/build_amr_integration_debug
./amrex_amr_backend_test
./amrex_amr_mesh_adapter_test

cd ~/mpm/demo
OMP_NUM_THREADS=1 ../build_amr_integration_debug/mpm \
  -f ~/mpm/demo/ -p 1 -i mpm-amr.json

OMP_NUM_THREADS=1 mpirun -np 2 \
  ~/mpm/build_amr_integration_debug/mpm \
  -f ~/mpm/demo/ -p 1 -i mpm-amr.json
```

## 24. 下一步建议

完成悬点凝聚后，剩余工作不再是“界面公式”本身，而是系统级完整性：

1. 做常数应力/常数速度 patch test，确认粗细分界面不会产生伪加速度。
2. 做粗到细、细到粗的一维应力波测试，比较界面反射率。
3. 把节点集、边界条件、集中载荷统一映射到 AMR `level/index/node` 坐标。
4. 为 AMR 输出单元级别、父单元 ID、悬点标记和约束权重。
5. 实现动态重划分时，需要同步重建 `HangingNodeConstraint`、MPM cell ID 映射和粒子所属单元。

## 25. 标准验证算例

两级立方体局部加密和正弦激励的标准算例位于：

```text
demo/amr_two_level_sine/
```

算例参数、运行方法、VTK 输出和应力曲线提取流程见：

```text
docs/AMReX两级立方体正弦激励标准算例说明.md
```
