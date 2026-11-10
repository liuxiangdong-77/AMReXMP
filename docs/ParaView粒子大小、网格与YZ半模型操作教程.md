# ParaView 粒子大小、网格与 YZ 半模型操作教程

如果只需要非常详细的 `Clip` 半模型步骤，请同时参考：

```text
docs/ParaView使用Clip查看YZ半边三维模型教程.md
```

## 1. 本教程解决什么问题

本教程对应标准算例：

```text
demo/amr_two_level_sine/
```

完成后可以在 ParaView 中：

1. 用 `radius_factor` 让不同 AMR 层的粒子显示成不同大小。
2. 查看 `grid_data0.vtu` 中的三级 Hex8 网格。
3. 用 X 法向的平面裁掉一半模型，从 YZ 方向观察应力波。
4. 用脚本自动完成上述视图设置。

## 2. 先重新运行算例

```bash
cd ~/mpm/demo/amr_two_level_sine

OMP_NUM_THREADS=4 \
~/mpm/build_amr_standard/mpm \
  -f ~/mpm/demo/amr_two_level_sine/ \
  -p 4 \
  -i mpm-two-level-sine.json
```

输出目录：

```text
results-amr-two-level-sine/amr-two-level-sine-000/
```

本版本会在其中生成：

```text
grid_data0.vtu
particle_data000.vtp
particle_data005.vtp
...
particle_data200.vtp
```

## 3. 粒子大小：radius_factor

### 3.1 必须打开组合文件

打开：

```text
particle_data000.vtp
```

不要用单独的 `stresses000.vtp` 设置粒子大小。原因是单字段文件只保存应力，不包含 `radius_factor`；组合文件 `particle_data000.vtp` 才同时包含：

```text
stresses
velocities
displacements
radius_factor
particle_spacing
```

### 3.2 设置 Point Gaussian

1. 在 Pipeline Browser 中选中 `particle_data000.vtp`。
2. 把 Representation 从 `Surface` 改为 `Point Gaussian`。
3. 勾选 `Use Scale Array`。
4. 在 `Gaussian Scale Array` 中选择：

```text
radius_factor
```

5. 将 `Gaussian Radius` 设置为约：

```text
0.12
```

6. 将 `Scale Factor` 设为：

```text
1.0
```

标准算例中的 `radius_factor` 为：

```text
Level 0: 1.00
Level 1: 0.50
Level 2: 0.25
```

因此细层粒子会显示为粗层粒子的 `1/2` 和 `1/4`。

如果粒子看起来仍然太小，只需增大 `Gaussian Radius`。如果某些粒子过大，降低 `Scale Factor`。

### 3.3 如果看不到 Use Scale Array

通常有以下原因：

1. 当前选中的是 `stresses000.vtp`，不是 `particle_data000.vtp`。
2. 当前 Representation 不是 `Point Gaussian`。
3. 当前打开的是旧结果文件。
4. 当前选中的是网格 `grid_data0.vtu`，它没有粒子点数据。

## 4. 网格显示：grid_data0.vtu

`grid_data0.vtu` 是真正的非结构六面体网格文件，不是粒子文件。

打开方式：

1. `File -> Open`。
2. 选择 `grid_data0.vtu`。
3. Property 面板中的 Representation 选择：

```text
Wireframe
```

或：

```text
Surface With Edges
```

4. Coloring 选择：

```text
amr_level
```

该文件包含的单元字段为：

```text
amr_level        AMR 层级 0/1/2
cell_size        单元尺寸 1.0/0.5/0.25
radius_factor    单元尺寸相对值 1.0/0.5/0.25
```

这样可以直接用颜色区分粗网格、Level 1 和 Level 2。

## 5. YZ 半模型：Clip 平面

“沿 YZ 平面切掉一半”应使用 `Clip`，不是简单的 `Slice`。

区别是：

```text
Slice：只保留二维切面
Clip：保留平面一侧的三维半模型
```

### 5.1 手动操作

1. 在 Pipeline Browser 中选中 `particle_data000.vtp`。
2. 点击 `Filters -> Common -> Clip`。
3. 点击 Apply。
4. Clip Type 选择：

```text
Plane
```

5. Origin 设置为：

```text
5.0  5.0  5.0
```

6. Normal 设置为：

```text
1.0  0.0  0.0
```

这个法向表示切割平面就是 YZ 平面，位置为：

```text
x = 5 m
```

7. `Invert` 不勾选时保留一侧，勾选后保留另一侧。
8. 把显示方式改为 `Point Gaussian`，Scale Array 选择 `radius_factor`。
9. 点击相机工具栏的 `-X` 或 `+X`，正对 YZ 截面观察。

### 5.2 网格也一起切半

对 `grid_data0.vtu` 重复同样的 Clip：

```text
Origin: 5 5 5
Normal: 1 0 0
```

然后把 Representation 设置为 `Wireframe`，Coloring 设置为 `amr_level`。

这样你会得到：

```text
半模型粒子 + 半模型网格线框
```

便于同时观察应力波和网格加密位置。

### 5.3 如果只想要 YZ 截面

选择 `Filters -> Common -> Slice`：

```text
Slice Type: Plane
Origin: 5 5 5
Normal: 1 0 0
```

Slice 会得到一个二维截面，不会保留另外半边三维体。

## 6. 自动配置 ParaView

仓库中新增了：

```text
tools/setup_paraview_amr_standard.py
```

它会自动：

1. 读取 `grid_data0.vtu`。
2. 读取 `particle_data*.vtp` 时间序列。
3. 将粒子设置为 `Point Gaussian`。
4. 用 `radius_factor` 控制粒子大小。
5. 对粒子和网格都建立 `x=5` 的 YZ 半模型 Clip。

### 6.1 在 ParaView Python Shell 中运行

先在 Windows 上准备一个可访问的数据目录，例如：

```text
D:/mpm_results/amr-two-level-sine-000/
```

在 ParaView 中打开：

```text
Tools -> Python Script Editor
```

然后执行：

```python
exec(open(r"D:/CODE/workspace/package/Ubuntu/mpm_live_snapshot/tools/setup_paraview_amr_standard.py", encoding="utf-8").read())

configure(
    r"D:/mpm_results/amr-two-level-sine-000",
    r"D:/mpm_results/amr-two-level-sine-000/amr_standard.pvsm",
    plane_x=5.0,
    particle_radius=0.12
)
```

### 6.2 使用 pvpython

```powershell
pvpython.exe setup_paraview_amr_standard.py `
  --data-dir "D:/mpm_results/amr-two-level-sine-000" `
  --state "D:/mpm_results/amr-two-level-sine-000/amr_standard.pvsm"
```

之后可以在 ParaView 中：

```text
File -> Load State
```

打开生成的 `.pvsm`。

## 7. 推荐观察顺序

1. 先打开 `grid_data0.vtu`，按 `amr_level` 看网格层级。
2. 对网格做 `x=5` Clip，观察 YZ 方向的加密块。
3. 打开 `particle_data*.vtp`，设置 Point Gaussian 和 `radius_factor`。
4. 对粒子做同样的 `x=5` Clip。
5. 时间播放，观察应力波从最细层向粗层传播。
6. 如果颜色动态范围太大，点击 `Rescale to Data Range Over All Timesteps`。

## 8. 常见问题

### 粒子全部一样大

确认 Scale Array 是 `radius_factor`，不是固定 Radius，也不是 `particle_spacing` 的某一个分量。

### 粒子消失或只剩少量点

检查是否设置了一个值为 0 的 Clip Origin/Normal，或打开了 `Invert` 后只保留了空的一侧。

### 网格看不见

把 `grid_data0.vtu` 的 Representation 改成 `Wireframe` 或 `Surface With Edges`，并暂时降低粒子大小。

### 截面上看不到波

`Slice` 只取严格位于平面上的粒子，点云可能很少。建议先用 `Clip` 保留一半三维模型，再从 `-X` 或 `+X` 方向观察。
