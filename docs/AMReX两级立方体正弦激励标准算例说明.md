# AMReX 两级立方体正弦激励标准算例说明

## 1. 算例目的

本算例用于快速检查以下功能是否协同工作：

1. 三维正交 Hex8 网格。
2. AMReX 两级 `2:1` 立方体局部加密。
3. 粗细界面悬点凝聚。
4. 单个粒子上的正弦激励荷载。
5. 应力场 VTK 时序输出。
6. 指定粒子的应力历史 CSV 提取。

算例目录：

```text
~/mpm/demo/amr_two_level_sine/
```

主要文件：

```text
mpm-two-level-sine.json   主输入文件
excitation-particle.txt   激励粒子坐标
```

## 2. 几何与网格

计算域为 `10 m x 10 m x 10 m` 的正方体：

```text
domain_min = [0, 0, 0]
domain_max = [10, 10, 10]
coarse_cells = [10, 10, 10]
```

因此粗网格尺寸为：

```text
h0 = 1.0 m
```

两级加密块分别为：

```text
Level 1: [2,2,2] 到 [8,8,8]
Level 2: [4,4,4] 到 [6,6,6]
```

对应的网格尺寸：

```text
Level 0: h0 = 1.00 m
Level 1: h1 = 0.50 m
Level 2: h2 = 0.25 m
```

所有加密块边界都落在父层单元边界上，满足当前有限元悬点约束的要求。

最终 finest valid cell 统计：

```text
3 个 AMR level
2960 个 MPM Cell
3789 个 MPM Node
936 个悬点
```

## 3. 粒子布置

### 3.1 背景粒子

背景区域使用输入文件里的 `gauss` 粒子生成器：

```json
{
  "type": "gauss",
  "nparticles_per_dir": 2,
  "particle_type": "P3D",
  "material_id": 0,
  "cset_id": -1,
  "pset_id": 0
}
```

每个 finest valid cell 内生成 `2 x 2 x 2 = 8` 个粒子：

```text
2960 x 8 = 23680 个背景粒子
```

### 3.2 激励粒子

激励粒子单独放在粒子集 `pset_id = 1`，坐标为：

```text
(5.125, 5.125, 5.125) m
```

该点位于最细加密区中心附近的小单元中心。

这里没有直接使用几何中心 `(5,5,5)`，原因是两级 `2:1` 加密块的对称中心正好落在有限元节点上。粒子如果放在节点上，会同时位于多个单元的公共边界，容易产生单元归属和形函数歧义。`(5.125,5.125,5.125)` 在几何中心附近，但明确位于一个最细单元内部。

## 4. 材料参数

算例采用线弹性岩石参数：

```text
密度            2650 kg/m^3
杨氏模量        20 GPa
泊松比          0.25
```

这些参数只用于标准算例验证。进行实际地下爆炸或岩石动力学分析时，应替换为相应材料模型和现场参数。

## 5. 正弦激励

激励作用在 `pset_id = 1` 的单个粒子上，方向为 `z`：

```json
{
  "pset_id": 1,
  "math_function_id": 0,
  "dir": 2,
  "traction": 5.0e6
}
```

时间函数为：

```json
{
  "id": 0,
  "type": "Sin",
  "x0": 0.0,
  "a": 3141.592653589793,
  "xrange": [0.0, 0.002]
}
```

对应：

```text
f = 500 Hz
omega = 2*pi*f = 3141.592653589793 rad/s
duration = 0.002 s
```

实际粒子面力为：

\[
t_z(t)=t_0\sin(\omega t),\qquad 0\leq t\leq 0.002
\]

其中：

```text
t0 = 5 MPa
```

## 6. 时间积分参数

```text
MPM scheme      USL
velocity update FLIP
dt              1.0e-5 s
nsteps          200
total time      0.002 s
output_steps    5
```

共输出 41 帧：

```text
initial frame + 40 time-history frames
```

如果只做快速调试，可以先把 `nsteps` 改成 `20`，把 `output_steps` 改成 `5`。

## 7. 运行方法

### 7.1 配置 AMReX 构建

```bash
cmake -S ~/mpm -B ~/mpm/build_amr_standard \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_PREFIX_PATH="$HOME/opt/amrex-26.11" \
  -DMPM_WITH_AMREX=ON \
  -DNO_KAHIP=True
```

### 7.2 编译

```bash
cmake --build ~/mpm/build_amr_standard --parallel 8
```

### 7.3 运行

```bash
cd ~/mpm/demo/amr_two_level_sine

OMP_NUM_THREADS=4 \
~/mpm/build_amr_standard/mpm \
  -f ~/mpm/demo/amr_two_level_sine/ \
  -p 4 \
  -i mpm-two-level-sine.json
```

本机实测使用 4 个 OpenMP 线程时，200 步约运行 170 秒。

### 7.4 在 VS Code 中运行

`.vscode/tasks.json` 已增加以下任务：

```text
AMR standard case: configure
AMR standard case: build
AMR standard case: run two-level sine
AMR standard case: extract stress history
AMR standard case: run then extract
```

打开 VS Code 的 `Terminal -> Run Task...`，选择：

```text
AMR standard case: run then extract
```

即可依次完成编译、计算和 `excitation_stress.csv` 提取。

## 8. 输出文件

输出目录：

```text
~/mpm/demo/amr_two_level_sine/
results-amr-two-level-sine/
amr-two-level-sine-000/
```

主要文件：

```text
particle_data000.vtp
particle_data005.vtp
particle_data010.vtp
...
particle_data200.vtp
```

每个组合文件包含：

```text
displacements
velocities
stresses
radius_factor
particle_spacing
```

其中 `stresses` 为粒子 Cauchy 应力张量。

完整 41 帧输出在本机约为：

```text
210 MB
```

如果不需要这么多帧，可以通过增大 `output_steps` 减少输出量。

## 9. 生成激励粒子应力曲线

构建会同时生成：

```text
~/mpm/build_amr_standard/particle_stress_history
```

该工具读取目录中的所有 `particle_data*.vtp`，选择每一帧中离目标坐标最近的粒子，并输出：

```text
step
time
sxx
syy
szz
sxy
syz
szx
von_mises
x
y
z
```

运行：

```bash
cd ~/mpm/demo/amr_two_level_sine

~/mpm/build_amr_standard/particle_stress_history \
  --dir results-amr-two-level-sine/amr-two-level-sine-000 \
  --x 5.125 \
  --y 5.125 \
  --z 5.125 \
  --dt 1e-5 \
  --output excitation_stress.csv
```

终端应显示：

```text
Wrote 41 frames to "excitation_stress.csv"
```

然后用 Excel、MATLAB、Python 或 ParaView 绘制：

```text
time -> von_mises
time -> szz
```

也可以直接把 CSV 的第 2 列作为横坐标，第 9 列作为纵坐标。

## 10. 本次运行结果

本次实际运行返回值为 0。

网格和悬点：

```text
3 levels
2960 cells
3789 nodes
936 hanging nodes
```

粒子：

```text
23680 个背景粒子
1 个激励粒子
```

应力历史共 41 帧。

激励粒子的 von Mises 应力：

```text
最大值约为 3082.66 Pa
最大值出现在 t = 0.0013 s
最终值约为 1505.08 Pa
```

这些数值不是材料强度判据，只用于确认正弦加载、AMR 界面和应力输出链路正常。

## 11. ParaView 查看方法

1. 打开一个 `particle_data000.vtp`。
2. 在时间序列中加载同目录下的全部 `particle_data*.vtp`。
3. 选择 `stresses` 作为着色变量。
4. 如果要看固定位置的时间曲线，选择激励粒子附近坐标
   `(5.125,5.125,5.125)`。
5. 使用 `Plot Selection Over Time` 查看该点随时间的应力变化。

需要注意：VTK 的 `stresses` 是 9 分量张量，不是单独的 `radius_factor` 标量。用于云图时可以查看 von Mises、最大主应力或指定分量。

## 12. 修改算例的入口

### 修改激励频率

修改 `math_functions[0].a`：

```text
a = 2*pi*f
```

例如 `f=1000 Hz`：

```text
a = 6283.185307179586
```

同时把 `analysis.nsteps` 或 `analysis.dt` 调整到至少覆盖一个周期。

### 修改激励幅值

修改 `particle_surface_traction[0].traction`，单位是 Pa。

### 修改激励方向

修改 `dir`：

```text
0: x 方向
1: y 方向
2: z 方向
```

### 修改加密区

修改 `mesh.amr.levels`。边界仍必须满足：

1. Level 1 边界对齐粗单元。
2. Level 2 边界对齐 Level 1 单元。
3. 相邻层差不得超过 1。

### 修改输出频率和长度

```json
"output_steps": 5
```

数值越小，输出帧越密。

```json
"nsteps": 200
```

用于控制总时间步数。

## 13. 验收标准

一次正常运行至少应满足：

1. 程序返回值为 0。
2. 日志报告 `3 levels` 和悬点数非零。
3. 输出目录中存在从 `particle_data000.vtp` 到
   `particle_data200.vtp` 的 41 个文件。
4. `excitation_stress.csv` 中存在 41 行数据。
5. `von_mises` 随时间出现明显非零变化。
6. 若使用双 MPI rank，应额外检查每个 rank 的 VTK 时序文件可以合并查看。

## 14. ParaView 高级显示

粒子 `radius_factor`、`grid_data0.vtu` 网格显示和 YZ 半模型 Clip 的完整操作见：

```text
docs/ParaView粒子大小、网格与YZ半模型操作教程.md
```
