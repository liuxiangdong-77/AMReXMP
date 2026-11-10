# AMReX 接入说明

本文对应分支 `C++20_AMReX` 的当前接入方式。逐行代码审查请同时阅读
[C++20 与 AMReX 直接接入修改说明](cxx20-amrex-migration.zh-CN.md)。

## 1. 当前结构

主工程已经统一使用 C++20，不再维护 C++14 主程序与 C++20 bridge 的混合构建。

```text
mpm C++20
  |
  +-- AmrexRuntime
        |
        +-- AMReX::amrex
```

AMReX 仍然是可选功能：

```text
-DMPM_WITH_AMREX=OFF  普通 C++20 MPM
-DMPM_WITH_AMREX=ON   C++20 MPM + AMReX
```

## 2. CMake 原理

`CMakeLists.txt` 第 5-7 行设置：

```cmake
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
```

`MPM_WITH_AMREX=ON` 时：

- 使用 `find_package(AMReX CONFIG REQUIRED)` 查找已安装 AMReX；
- 主程序直接链接 `AMReX::amrex`；
- 定义 `MPM_WITH_AMREX=1`；
- 构建独立的 `amrex_smoke_test`。

不再存在 `mpm_amrex_runtime` 静态库。

## 3. runtime 生命周期

`src/main.cc` 在 MPI 初始化之后创建：

```cpp
mpm::amr::AmrexRuntime amrex_runtime(argc, argv);
```

内部执行：

```cpp
amrex::Initialize(argc, argv, false, MPI_COMM_WORLD);
```

`false` 表示不解析 MPM 的 `-f/-p/-i` 参数。

对象析构时：

```cpp
amrex::Finalize(instance_);
```

生命周期顺序：

```text
MPI_Init
  -> amrex::Initialize
  -> MPM solve
  -> amrex::Finalize
  -> MPI_Finalize
```

## 4. 当前安装

```text
AMReX 版本：26.11
安装路径：/home/liuxiangdong/opt/amrex-26.11
静态库：  /home/liuxiangdong/opt/amrex-26.11/lib/libamrex_3d.a
MPI：     ON
OpenMP：  ON
Particles：ON
Linear solvers：OFF
```

## 5. 构建普通 C++20 MPM

```bash
cmake -S ~/mpm -B ~/mpm/build_cxx20 \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DMPM_WITH_AMREX=OFF \
  -DNO_KAHIP=True

cmake --build ~/mpm/build_cxx20 --parallel 8
```

## 6. 构建 C++20 MPM + AMReX

```bash
cmake -S ~/mpm -B ~/mpm/build_cxx20_amrex \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_PREFIX_PATH="$HOME/opt/amrex-26.11" \
  -DMPM_WITH_AMREX=ON \
  -DNO_KAHIP=True

cmake --build ~/mpm/build_cxx20_amrex --parallel 8
```

## 7. 验证

```bash
~/mpm/build_cxx20_amrex/amrex_smoke_test

mpirun --oversubscribe -np 2 \
  ~/mpm/build_cxx20_amrex/amrex_smoke_test

~/mpm/build_cxx20_amrex/mpm \
  -f ~/mpm/demo/ -p 8 -i mpm.json
```

运行 AMReX 版时，日志应出现：

```text
AMReX backend enabled: version 26.11
```

## 8. 后续如何使用 AMReX

下一步建议新增：

```text
include/mesh_adaptivity/amr_backend.h
src/mesh_adaptivity/amrex_backend.cc
```

AMReX 对象建议先映射为：

```text
amrex::Box                 一个规则 Hex8 逻辑块
amrex::BoxArray            某级网格的全部 patch
amrex::DistributionMapping MPI patch 归属
amrex::MultiFab            单元中心辅助数据
amrex::AmrMesh/AmrCore     层级网格和 regrid
amrex::TagBox              加密标记
amrex::ParticleContainer   粒子迁移和 MPI 通信
```

MPM 的应力、变形梯度和材料历史仍由 `ParticleBase` 保存。AMReX 不会自动完成
MPM 的粒子到节点映射、粗细界面动量传递或守恒修正。

第一版 AMR 建议实现：

```text
3D
两级网格
固定局部加密区
正交 Hex8
标准 MPM 单元
暂不动态 regrid
```

## 9. 注意事项

- 使用 `g++-13` 或更新版本。
- 不再使用旧的 `-fpermissive`。
- AMReX 编译需要 CMake 3.25 以上。
- 当前 `MPM_WITH_AMREX=OFF` 构建仍不依赖 AMReX。
