# 空单元网格与 ID 设计

## 1. 把 core 和 domain 设成相同尺寸会怎样

会生成均匀网格，但前提是把 `core.min/max` 和 `domain.min/max` 设置为相同值。
此时每条轴只有核心区，没有左右过渡区：

```text
core_intervals = round(domain_length / fine_spacing)
actual_spacing = domain_length / core_intervals
```

最后一步会重新均分长度，所以实际间距可能和 `fine` 有很小的浮点差异。
`coarse` 在该配置下不起作用，但建议仍设为与 `fine` 相同，避免配置歧义。

例如：

```json
"domain": {"min": [0, 0, 0], "max": [2, 2, 2]},
"core":   {"min": [0, 0, 0], "max": [2, 2, 2]},
"spacing": {
  "fine": [0.25, 0.25, 0.25],
  "coarse": [0.25, 0.25, 0.25],
  "transition": "linear",
  "max_adjacent_ratio": 1.0
}
```

每条轴会得到 8 个 0.25 单元。若长度不能被间距整除，生成器会按
`round(length / fine)` 取整数个单元，再略微缩放实际间距。

## 2. 空单元是怎样实现的

网格和粒子使用两个不同的盒子：

```text
domain          整个背景网格范围
particle_domain 允许生成粒子的范围
```

生成器始终建立完整的 `domain` 节点和 Hex8 单元。粒子中心只有落在
`particle_domain` 内才输出，所以 `domain` 比 `particle_domain` 大的部分就是
没有粒子的空单元。

当前 MPM 不需要每个单元都有粒子。单元是否有粒子可以通过 `nparticles()`
和 `status()` 判断。空单元不会从 `mesh.txt` 中删除，因为粒子以后进入这些
区域时，`locate_particles_mesh()` 仍然需要在这些单元中定位。

空单元只能缓冲网格域内部的大变形。如果粒子飞出完整 `domain`，程序仍会
把它视为无法定位并报错。

## 3. 节点编号

网格采用张量积顺序：

```text
node_id(i,j,k) = k*(nx+1)*(ny+1) + j*(nx+1) + i
```

其中 `i` 最快变化，然后 `j`，最后 `k`。节点文件按 `z -> y -> x` 嵌套输出，
因此文件行号减 1 就是节点 ID。

每个 Hex8 单元的节点顺序固定为：

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

这个顺序必须与 `ED3H8` 单元定义一致，否则 Jacobian 和积分点会错误。

## 4. 单元编号

```text
cell_id(i,j,k) = k*nx*ny + j*nx + i
```

脚本遍历顺序和 ID 公式完全一致。`particles-cells.txt` 中的 cell ID 是完整
背景网格的原始单元 ID，即使该单元后来成为空单元也不会重新编号。

## 5. 粒子编号

粒子 ID 不是原始张量积槽位编号，而是过滤后的连续编号：

1. 按 `z -> y -> x` 和单元内 `pz -> py -> px` 遍历；
2. 若粒子不在 `particle_domain` 或位于 Boolean 删除区，则不输出；
3. 只有实际写出的粒子获得下一个连续 ID；
4. `particles-cells.txt` 和 `particles-volumes.txt` 使用同一个新 ID。

因此删除粒子后，`particles-3d.txt` 不会保留空洞 ID。所有 entity set 和
`core-center-ids.json` 也都使用这个输出后的粒子 ID。

## 6. 相关文件

```text
mesh-3d.txt                 完整节点和所有背景单元
particles-3d.txt            过滤后的粒子坐标
particles-cells.txt         输出粒子 ID -> 原始单元 ID
particles-volumes.txt       输出粒子 ID -> 粒子体积
entity_sets.json            边界、核心区和中心点集合
core-center-ids.json        核心中心最近的节点和粒子
mesh-generation-report.json 空单元、删除粒子和数量统计
axis-spacing.csv            每条轴的实际节点坐标和单元尺寸
```

生成后建议运行：

```bash
python3 tools/graded_hex_mesh/validate_graded_hex_output.py generated_graded_hex
```
