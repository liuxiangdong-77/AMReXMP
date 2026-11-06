# 组合粒子输出与 radius_factor

## 1. 输出目标

原先每个 VTK 字段分别写入 `.vtp`。后处理时应通过
`particle_data<step>.vtp` 查看组合数据，其中同一粒子同时具有：

- `stresses`
- `velocities`
- `displacements`
- `radius_factor`
- `particle_spacing`
- 配置中请求的其他标量、向量、张量和状态变量

单独字段文件仍会保留，用于兼容旧流程。

## 2. VtkWriter 的组合数据机制

`src/io/vtk_writer.cc` 构造函数创建固定的 `points_` 和
`combined_data_` 两个 `vtkPolyData` 对象。

新增接口：

```cpp
add_scalar_point_data(...);
add_vector_point_data(...);
add_tensor_point_data(...);
write_combined_point_data(...);
```

前三个函数把字段追加到 `combined_data_`，最后一个函数只把这份组合数据写出
一次。标量、向量和张量仍分别构造原来的单字段数据集，所以新旧输出并存。

## 3. MPMBase 中的修改

`include/solvers/mpm_base.tcc` 中的修改包括：

1. 在变量注册表加入 `radius_factor`（Scalar）和
   `particle_spacing`（Vector）。
2. 解析 `post_processing.vtk` 时按变量类型放入 `vtk_vars_`。
3. 写每个单独字段时同时调用对应的 `add_*_point_data`。
4. 遍历完所有请求字段后，检查 `radius_factor` 和 `particle_spacing`
   是否已经加入组合数据；若没有则自动追加。
5. 生成 `particle_data<step>.vtp` 并一次性写出组合数据。

自动追加使这两个诊断字段不必出现在 `post_processing.vtk` 中；但建议仍明确
写出，便于阅读输入文件。

## 4. Mesh 中的数据提取

`include/mesh/mesh.tcc` 增加两个特殊分支。

### radius_factor

```text
radius_factor = 当前单元等效尺寸 / 全局最大单元等效尺寸
```

单元等效尺寸使用单元体积的 `Tdim` 次方根。MPI 运行时用 `MPI_Allreduce`
取得全局最大值。找不到单元时写入 NaN。

该值来自旧多级网格可视化方案，适合各方向等比例缩放的粒子。对于三个方向
尺寸不同的网格，一个标量无法准确表达三个方向的粒子尺寸。

### particle_spacing

```text
particle_spacing[x] = particle_fraction * cell_extent_x
particle_spacing[y] = particle_fraction * cell_extent_y
particle_spacing[z] = particle_fraction * cell_extent_z
particle_fraction = (particle_volume / cell_volume)^(1/Tdim)
```

因此它输出三个方向的实际粒子间距。ParaView 中建议使用
`3D Glyphs -> Box -> Scale by Vector Components -> particle_spacing`，
而不是 Point Gaussian。

## 5. demo/mpm.json

示例输入在 `post_processing.vtk` 中加入了：

```json
"radius_factor",
"particle_spacing"
```

这用于展示字段注册方式。实际计算可继续选择需要的速度、位移、应力等字段。

## 6. 增加新输出字段的步骤

1. 在 `MPMBase` 的 `variables` 中登记字段类型。
2. 在 `mesh.tcc` 增加数据来源；普通粒子字段可能已经由
   `ParticleBase::scalar_data/vector_data` 自动支持。
3. 在 `write_vtk()` 中确认该类型已调用对应的 `add_*_point_data`。
4. 在 `mpm.json` 的 `post_processing.vtk` 中声明字段。
5. 打开 `particle_data<step>.vtp`，同时检查单字段文件，确认字段数量和
   粒子数量一致。

不要把派生量直接写入粒子历史状态；如果它只用于后处理，优先在
`Mesh::particles_*_data()` 中计算，避免影响求解器内存和状态。
