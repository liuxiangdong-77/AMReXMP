# 渐变 Hex8 网格生成器使用教程

这个工具用于在 MPM 初始化阶段生成一张：

- 全正交六面体 Hex8 网格；
- 核心区域使用最小单元尺寸；
- 从核心区边界开始，沿 x、y、z 三个方向分别渐变；
- 整个计算域都参与尺寸过渡；
- 无悬挂节点；
- 可直接由 CB-Geo MPM 的 `IOMeshAscii` 读取。

工具只使用 Python 标准库，不需要 MATLAB、NumPy 或其他第三方包。

## 1. 目录内容

```text
tools/graded_hex_mesh/
  generate_graded_hex_mesh.py     主生成程序
  validate_graded_hex_output.py   生成结果校验程序
  graded_hex.example.json         小型示例配置
  test_generate_graded_hex_mesh.py 单元测试
  vscode-task.example.json        VS Code 任务示例
  README.md                       English quick reference
  USAGE.zh-CN.md                  本教程
```

## 2. 快速运行

进入工具目录：

```bash
cd ~/mpm/tools/graded_hex_mesh
```

先做 dry-run，只检查配置和估算网格规模，不写文件：

```bash
python3 generate_graded_hex_mesh.py graded_hex.example.json --dry-run
```

正式生成：

```bash
python3 generate_graded_hex_mesh.py graded_hex.example.json
```

校验生成结果：

```bash
python3 validate_graded_hex_output.py generated_graded_hex
```

如果输出目录中已有同名文件，使用 `--force` 覆盖生成文件：

```bash
python3 generate_graded_hex_mesh.py graded_hex.example.json --force
```

也可以覆盖输出目录：

```bash
python3 generate_graded_hex_mesh.py graded_hex.example.json \
  --output-dir /home/liuxiangdong/mpm/demo/graded_example
```

## 3. 配置结构

示例配置为 `graded_hex.example.json`。

### 3.1 计算域

```json
"domain": {
  "min": [0.0, 0.0, 0.0],
  "max": [30000.0, 30000.0, 35000.0]
}
```

这里定义整个背景网格范围。

### 3.2 核心加密区

```json
"core": {
  "min": [0.0, 0.0, 26000.0],
  "max": [500.0, 500.0, 27000.0]
}
```

核心区域内使用 `spacing.fine`。

核心区必须完全位于 `domain` 内，但可以与计算域边界重合。

例如：

```text
core.min.x = domain.min.x
```

表示 x 方向没有左侧过渡区，只在核心区右侧进行渐变。

### 3.3 粒子生成域与空单元缓冲层

网格域 `domain` 和粒子生成域 `particle_domain` 可以不同。例如：

```json
"domain": {
  "min": [0.0, 0.0, 0.0],
  "max": [30000.0, 30000.0, 35000.0]
},
"particle_domain": {
  "min": [1000.0, 1000.0, 1000.0],
  "max": [29000.0, 29000.0, 34000.0]
}
```

这样整张 Hex8 网格仍覆盖 0 到 30000，但只在内部区域生成粒子。
网格四周的单元会保留为空单元。

MPM 不要求每个单元必须包含粒子。空单元保持非活动状态，粒子进入空单元后
仍可通过 `locate_particles_mesh()` 正常找到该单元。

空单元只能提供有限缓冲层。如果粒子飞出整个 `domain`，MPM 仍然会报错，
因此缓冲层必须覆盖预计的大变形范围。

### 3.4 网格尺寸

```json
"spacing": {
  "fine": [25.0, 25.0, 25.0],
  "coarse": [800.0, 800.0, 800.0],
  "transition": "linear",
  "max_adjacent_ratio": 1.3
}
```

- `fine`：核心区 x、y、z 方向的最小单元尺寸。
- `coarse`：距离核心区最远处的目标最大尺寸。
- `transition`：`linear` 或 `smoothstep`。
- `max_adjacent_ratio`：允许的相邻单元尺寸比。

按你的要求，`linear` 使用：

```text
N = round(2L / (hf + hc))
h_i = hf + (hc - hf) * i / (N - 1)
```

其中：

- `L`：从核心区边界到计算域边界的距离；
- `hf`：核心区细尺寸；
- `hc`：远端粗尺寸；
- `N`：该侧过渡单元的个数。

由于浮点数和整数取整，最后一个尺寸会轻微偏离 `hc`。生成器会对所有尺寸做统一缩放，使总长度严格等于 `L`。

`smoothstep` 使用相同的离散编号，但把线性权重 `i/(N-1)` 替换为平滑权重：

```text
t = i / (N - 1)
w = 3t² - 2t³
h_i = hf + (hc - hf) * w
```

## 4. 核心区与过渡方式

三个方向分别生成节点坐标：

```text
x_nodes = graded_axis(domain, core, fine_x, coarse_x)
y_nodes = graded_axis(domain, core, fine_y, coarse_y)
z_nodes = graded_axis(domain, core, fine_z, coarse_z)
```

然后生成张量积网格：

```text
all_nodes = x_nodes × y_nodes × z_nodes
```

核心区和左右两侧分别处理：

```text
domain.min  <-- 左侧过渡 -->  core.min
core.min    <-- 核心细网格 --> core.max
core.max    <-- 右侧过渡 -->  domain.max
```

核心区不强制使用 2:1 跳变，整个距离都参与尺寸变化。

## 5. 单元类型与节点顺序

生成器只输出：

```text
cell_type = ED3H8
node_type = N3D
io_type = Ascii3D
isoparametric = false
```

每个单元的节点顺序为：

```text
0: (x-, y-, z-)
1: (x+, y-, z-)
2: (x+, y+, z-)
3: (x-, y+, z-)
4: (x-, y-, z+)
5: (x+, y-, z+)
6: (x+, y+, z+)
7: (x-, y+, z+)
```

节点编号从 0 开始。

## 6. 粒子生成

```json
"particles": {
  "per_cell": [2, 2, 2],
  "particle_type": "P3D",
  "material_id": 0,
  "pset_id": 0
}
```

`per_cell = [2, 2, 2]` 表示每个六面体内部生成 `2 × 2 × 2 = 8` 个粒子。

每个粒子的体积为：

```text
cell_volume / 8
```

因此细单元产生小体积粒子，粗单元产生大体积粒子，粒子间距会随局部单元尺寸变化。

## 7. Boolean 粒子删除

Boolean 不会切割背景六面体，只删除位于几何体内部的粒子。背景网格仍保持正交协调。

本功能默认关闭。关闭方式是把：

```json
"remove_regions": []
```

保留为空数组，或者直接省略 `remove_regions`。示例配置已经使用空数组。
需要启用时，在 `remove_regions` 中加入下面的几何体即可。

### 球体

```json
{
  "type": "sphere",
  "center": [0.0, 0.0, 26500.0],
  "radius": 179.0
}
```

### 圆柱体

```json
{
  "type": "cylinder",
  "center": [0.0, 0.0, 26500.0],
  "radius": 179.0,
  "axis": "z",
  "height": 245.0
}
```

### 长方体

```json
{
  "type": "box",
  "min": [0.0, 0.0, 0.0],
  "max": [100.0, 100.0, 100.0]
}
```

多个 `remove_regions` 按并集处理。

## 8. 核心中心节点和粒子

每次生成都会输出：

```text
core-center-ids.json
```

内容包括：

- 核心区精确中心坐标；
- 距离中心最近节点的 ID 和坐标；
- 距离中心最近的保留粒子 ID 和坐标；
- 中心节点所在的 `node_set` ID；
- 中心粒子所在的 `particle_set` ID。

默认 set ID 都是 100，可以在配置中修改：

```json
"entity_sets": {
  "center_node_set_id": 100,
  "center_particle_set_id": 100
}
```

生成器会自动把中心节点加入 `node_sets`，把中心粒子加入 `particle_sets`。
因此可以直接在 `mpm.json` 中使用这两个 set ID 施加节点力、速度约束或粒子载荷。

如果核心区内的所有粒子都被 Boolean 删除，`closest_particle` 会为 `null`。
中心节点仍然会正常输出。

## 9. 生成文件

生成目录中会出现：

```text
mesh-3d.txt
particles-3d.txt
particles-cells.txt
particles-volumes.txt
entity_sets.json
mpm-input-fragment.json
mesh-generation-report.json
axis-spacing.csv
core-center-ids.json
```

说明：

- `mesh-3d.txt`：MPM 网格文件。
- `particles-3d.txt`：粒子坐标。
- `particles-cells.txt`：粒子 ID 和所属单元 ID。
- `particles-volumes.txt`：粒子 ID 和体积。
- `entity_sets.json`：边界节点集合。
- `mpm-input-fragment.json`：可直接合并到 `mpm.json` 的 `mesh` 和 `particles` 配置。
- `mesh-generation-report.json`：节点数、单元数、粒子数、最小尺寸、最大相邻尺寸比和建议时间步。
- `axis-spacing.csv`：三个方向的节点坐标和单元尺寸。
- `core-center-ids.json`：核心中心附近节点和粒子的 ID、坐标及 set ID。

## 10. 接入 mpm.json

打开 `mpm-input-fragment.json`，将其中的 `mesh` 和 `particles` 内容合并到算例的 `mpm.json`。

主要字段应类似：

```json
"mesh": {
  "mesh": "mesh-3d.txt",
  "entity_sets": "entity_sets.json",
  "io_type": "Ascii3D",
  "check_duplicates": false,
  "isoparametric": false,
  "node_type": "N3D",
  "cell_type": "ED3H8",
  "particle_cells": "particles-cells.txt",
  "particles_volumes": "particles-volumes.txt"
}
```

如果生成文件与 `mpm.json` 不在同一目录，需要填写相对于 `-f` 工作目录的路径。

## 11. 时间步

MPM 当前不会根据最小单元尺寸自动计算稳定时间步。

生成报告会给出：

```json
"time_step": {
  "suggested_dt": 0.0123
}
```

建议把 `mpm.json` 中的：

```json
"analysis": {
  "dt": 0.0123
}
```

设置为不大于报告中的值。

报告中建议时间步的计算方式为：

```text
dt = CFL * min_cell_size / (wave_speed + max_velocity)
```

如果材料波速或当前速度更大，应进一步减小 `dt`。

## 12. 输出 radius_factor 到 VTK

MPM 已经支持把渐变网格的粒子半径因子输出到 VTK。在算例的
`post_processing.vtk` 中加入：

```json
"post_processing": {
  "path": "results/",
  "vtk": ["radius_factor"],
  "output_steps": 10
}
```

`radius_factor` 的定义为：

```text
radius_factor = 当前单元等效尺寸 / 全局最大单元等效尺寸
```

因此：

- 最粗网格区域的粒子约为 `1.0`；
- 最细网格区域的粒子小于 `1.0`；
- 在 ParaView 中可以用 `radius_factor` 控制粒子 Glyph 的缩放。

等效尺寸按单元体积计算，三维为：

```text
h_cell = cell_volume^(1/3)
```

MPI 计算时会使用所有 rank 的全局最大单元尺寸作为归一化基准。

### 12.1 合并的 ParaView 数据文件

MPM 以前会把每个变量写成单独的 VTP 文件，例如 `stresses1.vtp` 和
`radius_factor1.vtp`。因此打开 `stresses1.vtp` 时，Point Gaussian
看不到其他文件里的 `radius_factor`。

现在每次 VTK 输出还会写一个合并文件：

```text
particle_data0.vtp
particle_data1.vtp
...
```

该文件包含本次请求的所有标量、向量和张量，并额外包含：

- `radius_factor`：标量归一化粒子尺寸；
- `particle_spacing`：三分量的实际粒子间距。

如果各方向单元尺寸不同，仅使用一个标量无法正确填满各向异性的单元。
为了避免 Point Gaussian 在某些方向上重叠、另一些方向出现空隙，建议用
`3D Glyphs`：

1. 打开 `particle_data1.vtp`。
2. Representation 选择 `3D Glyphs`。
3. Glyph Type 选择 `Box`。
4. Scale Mode 选择 `Scale by Vector Components`。
5. Scale Array 选择 `particle_spacing`。
6. 用 `stresses` 或其他变量控制颜色。

`particle_spacing` 的三个分量由单元 x、y、z 方向尺寸和该单元粒子数量
共同计算。对于生成器的 `2×2×2` 粒子布局，Box 的三边分别等于局部
`dx/2`、`dy/2`、`dz/2`，可以无缝填满单元。

## 13. 使用 VS Code 任务

可以把 `vscode-task.example.json` 中的任务合并到仓库的 `.vscode/tasks.json`。

之后在 VS Code 中运行：

```text
Tasks: Run Task
Generate graded Hex8 mesh
```

任务会以 `tools/graded_hex_mesh` 为工作目录调用生成器。

## 14. 常见问题

### 14.1 报错“requested fine/coarse size range is too large”

这说明过渡距离不足，无法同时满足指定的最细尺寸、最粗尺寸和最大相邻尺寸比。

可以：

- 增大核心区到计算域边界的距离；
- 减小 `spacing.coarse`；
- 放宽 `max_adjacent_ratio`；
- 增大 `spacing.fine`。

### 14.2 网格数量过大

张量积网格的单元数是：

```text
nx * ny * nz
```

三个方向同时加密会快速增加总数。先运行：

```bash
python3 generate_graded_hex_mesh.py config.json --dry-run
```

查看 `nnodes`、`ncells` 和 `max_slots` 后再正式生成。

### 14.3 粒子数为零或过少

检查：

- Boolean 是否删除了整个区域；
- `domain` 与 `core` 是否写反；
- `particles.per_cell` 是否为正整数。

### 14.4 生成的节点集编号不符合算例

修改配置中的：

```json
"entity_sets": {
  "x_min_id": 1,
  "x_max_id": 2,
  ...
}
```

这些 ID 应与 `mpm.json` 中的边界条件 `nset_id` 对应。
