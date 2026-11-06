# AMReX 接入说明

本文对应 `/home/liuxiangdong/mpm` 当前代码，行号核对日期为 2026-11-05。

## 1. 这次接入完成了什么

这次完成的是“能让 MPM 调用 AMReX”的基础设施，不是已经完成 AMR 网格迁移。

具体包括：

- AMReX 是可选依赖，默认构建不需要安装 AMReX。
- `MPM_WITH_AMREX=ON` 时，CMake 查找并链接 AMReX 26.11。
- MPM 启动时初始化 AMReX，程序结束时关闭 AMReX。
- AMReX 复用 MPM 已经初始化的 MPI，并在 `MPI_Finalize` 之前关闭。
- 增加独立的 `amrex_smoke_test`，实际包含 AMReX 头文件、创建
  `amrex::Box` 并验证链接成功。
- 现有 `Mesh`、`Cell`、`MPMExplicit`、`MPMImplicit`、粒子映射和材料模型
  都没有改成 AMReX 数据结构。

因此现在的结果是：

```text
MPM 主程序可以创建 AMReX runtime，也可以继续调用任意 AMReX API；
但背景网格仍然是原来的 mesh.txt 单层网格。
```

## 2. 引入 AMReX 的整体原理

工程接入一个外部 C++ 库，通常包含四层：

```text
依赖发现：CMake find_package
构建链接：target_link_libraries
生命周期：Initialize / Finalize
代码隔离：MPM 自己的接口只暴露通用类型
```

本工程额外做了第五层处理：

```text
C++ 标准隔离：只有 AMReX bridge 使用 C++20，主 MPM 继续使用 C++14
```

原因是 AMReX 26.11 要求 C++20，而当前 MPM 和内置 TCLAP 的旧构造语法
不能直接整体迁移到 C++20。如果强行把 `mpm` 目标全部改成 C++20，会先产生
大量与 AMReX 无关的编译错误。

因此当前结构是：

```text
mpm                         C++14
  |
  +-- mpm_amrex_runtime     C++20
        |
        +-- AMReX::amrex
```

主程序只包含 `include/mesh_adaptivity/amrex_runtime.h`。该头文件不包含任何
AMReX 头文件，所以主 MPM 不需要以 C++20 编译。真正的 AMReX 类型只出现在
`src/mesh_adaptivity/amrex_runtime.cc` 和 smoke test 中。

## 3. 运行时生命周期

当前启动顺序是：

```text
main
  -> MPI_Init
  -> 创建 AmrexRuntime
       -> amrex::Initialize(..., false, MPI_COMM_WORLD)
  -> 创建 IO
  -> 创建并运行 MPM solver
  -> AmrexRuntime 析构
       -> amrex::Finalize
  -> MPI_Finalize
```

`amrex::Initialize` 的第三个参数为 `false`，表示 AMReX 不解析 MPM 的
命令行参数。这样 AMReX 不会消费 MPM 自己的：

```text
-f 工作目录
-p OpenMP 线程数
-i 输入 JSON
```

MPM 已经调用 `MPI_Init` 时，AMReX 会复制
`MPI_COMM_WORLD`，不会重复初始化 MPI。AMReX 析构发生在
`MPI_Finalize` 之前，这是当前生命周期设计中最重要的一点。

## 4. 逐文件、逐行修改说明

### 4.1 `CMakeLists.txt`

#### 第 60-62 行：定义可选功能开关

```cmake
option(MPM_WITH_AMREX "Build MPM with the optional AMReX backend" OFF)
```

作用：

- 默认值为 `OFF`，原构建命令保持不变。
- 只有显式传入 `-DMPM_WITH_AMREX=ON` 才查找 AMReX。
- 避免所有用户都被强制要求安装 AMReX。

#### 第 90-104 行：查找并检查 AMReX

```cmake
if(MPM_WITH_AMREX)
  ...
  find_package(AMReX CONFIG REQUIRED)
  ...
endif()
```

逐行作用：

- 第 91 行：只在开关打开时进入 AMReX 配置。
- 第 92-95 行：AMReX 26.11 要求 CMake 3.25 以上。
- 第 96 行：使用 AMReX 安装目录中的 `AMReXConfig.cmake`。
- 第 97-98 行：打印实际找到的 AMReX 目录和 MPI 配置。
- 第 99-103 行：如果 MPM 有 MPI，但 AMReX 是无 MPI 版本，给出警告。

`CONFIG` 很重要，它表示使用 AMReX 自己安装的 CMake package，而不是让
工程猜测头文件和库路径。

#### 第 247 行：统一 `target_link_libraries` 调用形式

```cmake
target_link_libraries(mpm PRIVATE git)
```

原代码使用的是不带 `PRIVATE/PUBLIC` 的旧式写法。后面又要给 `mpm` 链接
`mpm_amrex_runtime`；CMake 不允许同一个 target 混用两种签名，因此这一行
改为关键字形式。这个修改在 AMReX 关闭时同样有效。

#### 第 249-250 行：定义主程序的编译宏

```cmake
if(MPM_WITH_AMREX)
  target_compile_definitions(mpm PRIVATE MPM_WITH_AMREX=1)
```

作用：

- 只有打开 AMReX 时才编译 `main.cc` 中的 AMReX runtime 代码块。
- 不把 `MPM_WITH_AMREX` 暴露给其他 target。

#### 第 252-260 行：创建 C++20 AMReX bridge

```cmake
add_library(mpm_amrex_runtime STATIC ...)
target_include_directories(mpm_amrex_runtime PUBLIC include)
target_compile_features(mpm_amrex_runtime PRIVATE cxx_std_20)
target_compile_definitions(mpm_amrex_runtime PRIVATE MPM_WITH_AMREX=1)
target_link_libraries(mpm_amrex_runtime PRIVATE AMReX::amrex)
```

逐行作用：

- 第 252-254 行：创建一个独立静态库，只编译 `amrex_runtime.cc`。
- 第 255-257 行：让主工程能找到 `amrex_runtime.h`。
- 第 258 行：只有这个 bridge 使用 C++20。
- 第 259 行：bridge 源文件内部启用 AMReX 实现。
- 第 260 行：bridge 链接 AMReX 官方 imported target。

这里使用 `PRIVATE` 有两个目的：

- AMReX 头文件不会传给主 MPM target。
- AMReX 的 C++20 编译要求不会扩散到整个项目。

#### 第 262 行：把 bridge 链接给 MPM

```cmake
target_link_libraries(mpm PRIVATE mpm_amrex_runtime)
```

这样 `mpm` 可以调用 `AmrexRuntime`，但仍然不直接包含 AMReX 头文件。

#### 第 264-274 行：创建 AMReX smoke test

```cmake
add_executable(amrex_smoke_test ...)
target_compile_features(amrex_smoke_test PRIVATE cxx_std_20)
target_link_libraries(amrex_smoke_test PRIVATE
  mpm_amrex_runtime
  AMReX::amrex
)
```

作用：

- smoke test 自己使用 C++20。
- 它直接包含 AMReX 头文件，所以能验证 AMReX 的公开 API 是否可用。
- 它链接 runtime bridge 和 AMReX，验证库文件、MPI 和 OpenMP 链接。

### 4.2 `src/main.cc`

#### 第 14 行：包含 MPM 自己的 bridge 头文件

```cpp
#include "mesh_adaptivity/amrex_runtime.h"
```

这里没有包含 `AMReX.H`，因此主程序继续保持 C++14 兼容。

#### 第 41-45 行：创建 AMReX runtime

```cpp
#ifdef MPM_WITH_AMREX
    mpm::amr::AmrexRuntime amrex_runtime(argc, argv);
#endif
```

作用：

- 该对象的生命周期覆盖整个 `try` 块中的 MPM 求解过程。
- 离开 `try` 块时自动析构并关闭 AMReX。
- 对应宏未定义时，这一整段代码不参与编译。

位置选择在 `MPI_Init` 之后、MPM solver 之前。因此 AMReX 能看到已经存在的
MPI 通信器，同时能在 `MPI_Finalize` 之前完成关闭。

#### 第 56-59 行：输出 AMReX 版本

```cpp
console->info("AMReX backend enabled: version {}",
              mpm::amr::AmrexRuntime::version());
```

作用：

- 用于确认运行的确实是 AMReX 版本，而不是误用了普通构建。
- 版本查询通过 bridge 调用，不在 `main.cc` 中直接引用 AMReX API。

### 4.3 `include/mesh_adaptivity/amrex_runtime.h`

#### 第 4-5 行：只包含标准 C++ 头文件

```cpp
#include <memory>
#include <string>
```

该头文件故意不包含 AMReX 头文件，保证 C++14 主工程可以使用它。

#### 第 10-12 行：定义 bridge 类

```cpp
class AmrexRuntime {
```

它表示 MPM 对 AMReX runtime 的最小抽象，而不是 AMReX 网格类。

#### 第 14-15 行：构造和析构接口

```cpp
AmrexRuntime(int& argc, char**& argv);
~AmrexRuntime();
```

构造函数接收命令行参数，析构函数负责 AMReX 清理。

#### 第 17-18 行：禁止拷贝

```cpp
AmrexRuntime(const AmrexRuntime&) = delete;
AmrexRuntime& operator=(const AmrexRuntime&) = delete;
```

AMReX runtime 是全进程级资源，不应该被复制。

#### 第 20-24 行：状态和版本接口

```cpp
bool active() const noexcept;
static std::string version();
```

- `active()` 用于 smoke test 检查是否成功初始化。
- `version()` 用于日志和构建验证。

#### 第 26-28 行：Pimpl

```cpp
class Impl;
std::unique_ptr<Impl> impl_;
```

Pimpl 让 AMReX 的完整类型只存在于 `.cc` 文件中。这样头文件可以给 C++14
目标使用，而实现文件由 C++20 编译。

### 4.4 `src/mesh_adaptivity/amrex_runtime.cc`

#### 第 5-7 行：只在 bridge 实现中包含 AMReX

```cpp
#ifdef MPM_WITH_AMREX
#include <AMReX.H>
#endif
```

这是整个主 MPM 代码中唯一需要 AMReX 核心头文件的 runtime 文件。

#### 第 12-17 行：保存 AMReX 实例指针

```cpp
class AmrexRuntime::Impl {
 public:
  amrex::AMReX* instance{nullptr};
};
```

AMReX 26.11 的 `Initialize` 返回一个非拥有指针；该指针必须传回
`amrex::Finalize`。

#### 第 19-34 行：初始化 AMReX

```cpp
impl_->instance = amrex::Initialize(argc, argv, false, MPI_COMM_WORLD);
amrex::SetVerbose(0);
```

逐行作用：

- 第 22 行：判断 AMReX 本身是否启用了 MPI。
- 第 23 行：MPI 已由 MPM 初始化时，把 `MPI_COMM_WORLD` 交给 AMReX。
- 第 25 行：AMReX 无 MPI 时使用无通信器重载。
- 第 27 行：关闭 AMReX 自身的常规输出。
- 第 28-32 行：如果 bridge 被错误地在无 AMReX 构建中调用，抛出明确异常。

第三个参数 `false` 禁止 AMReX 解析命令行，避免破坏 MPM 参数。

#### 第 36-41 行：关闭 AMReX

```cpp
amrex::Finalize(impl_->instance);
```

该函数在 `AmrexRuntime` 析构时调用。由于对象位于 `main` 的 `try` 块内，
它一定先于第 99 行的 `MPI_Finalize` 执行。

#### 第 43-49 行：检查 runtime 状态

```cpp
return impl_ != nullptr && impl_->instance != nullptr;
```

只有 AMReX 已编译并且初始化成功时才返回 `true`。

#### 第 51-57 行：版本查询

```cpp
return amrex::Version();
```

用于 `main` 日志和 smoke test。

### 4.5 `tests/mesh_adaptivity/amrex_smoke_test.cc`

#### 第 3-4 行：直接包含 AMReX API

```cpp
#include <AMReX.H>
#include <AMReX_Box.H>
```

这是有意设计的：runtime bridge 隐藏 AMReX，但 smoke test 必须直接使用
AMReX，才能证明主工程真的能找到并链接 AMReX 头文件和库。

#### 第 9 行：初始化 runtime

```cpp
mpm::amr::AmrexRuntime runtime(argc, argv);
```

#### 第 11-12 行：创建 AMReX Box

```cpp
amrex::Box box(amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
               amrex::IntVect(AMREX_D_DECL(3, 3, 3)));
```

这是一个最小的 AMReX 数据结构，用来验证 `AMReX_Box.H` 和库符号可用。

#### 第 14-18 行：输出并返回测试结果

```cpp
return (runtime.active() && box.ok()) ? 0 : 1;
```

返回 0 表示 runtime、头文件和链接都正常。

### 4.6 `tools/amrex/setup_amrex.sh`

#### 第 9-17 行：配置路径和构建参数

```bash
install_prefix=${2:-"$HOME/opt/amrex-26.11"}
build_dir=${3:-"$HOME/opt/build-amrex-26.11"}
cc=${AMREX_CC:-gcc-13}
cxx=${AMREX_CXX:-g++-13}
```

作用：

- 默认使用 `g++-13`，因为系统默认 `g++` 是 9.5，不能编译 AMReX 26.11。
- 默认安装到 `~/opt/amrex-26.11`。
- 允许通过环境变量覆盖编译器和构建参数。

#### 第 24-34 行：配置 AMReX

```bash
-DAMReX_SPACEDIM=3
-DAMReX_MPI=ON
-DAMReX_OMP=ON
-DAMReX_PARTICLES=ON
-DAMReX_LINEAR_SOLVERS=OFF
```

作用：

- `SPACEDIM=3`：当前目标是 3D MPM。
- `MPI=ON`：为多 MPI rank 做准备。
- `OMP=ON`：与 MPM 的 OpenMP 并行一致。
- `PARTICLES=ON`：为后续 MPM 粒子容器迁移准备。
- `LINEAR_SOLVERS=OFF`：本阶段不需要 AMReX 线性求解器。

#### 第 36-37 行：编译并安装

```bash
cmake --build "$build_dir" --parallel "$build_jobs"
cmake --install "$build_dir"
```

当前脚本默认 `AMREX_BUILD_JOBS=4`。这台虚拟机的 `g++-13.3` 在 AMReX
`-O3` 下出现过内部编译器错误，因此默认 `AMREX_RELEASE_FLAGS="-O1 -DNDEBUG"`。

### 4.7 `README.md`

#### 第 17-21 行：增加本地文档入口

```markdown
Local extension notes:
- AMReX optional backend
- Graded Hex8 empty cells and ID design
- Combined particle output and radius_factor
```

作用只是让新增文档更容易找到，不影响构建。

## 5. 当前实际安装状态

当前虚拟机中已经安装：

```text
AMReX 版本：26.11
安装目录：/home/liuxiangdong/opt/amrex-26.11
静态库：  /home/liuxiangdong/opt/amrex-26.11/lib/libamrex_3d.a
包含：    AMR、Particles、MPI、OpenMP
未启用：  AMReX Linear Solvers
```

验证构建目录：

```text
/home/liuxiangdong/mpm/build_codex_plain
/home/liuxiangdong/mpm/build_codex_amrex2
```

## 6. 如何重新编译

### 6.1 原版构建

```bash
cmake -S ~/mpm -B ~/mpm/build_plain \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DMPM_WITH_AMREX=OFF \
  -DNO_KAHIP=True
cmake --build ~/mpm/build_plain --parallel 8
```

### 6.2 AMReX 构建

```bash
cmake -S ~/mpm -B ~/mpm/build_amrex \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_PREFIX_PATH="$HOME/opt/amrex-26.11" \
  -DMPM_WITH_AMREX=ON \
  -DNO_KAHIP=True
cmake --build ~/mpm/build_amrex --parallel 8
```

### 6.3 运行 smoke test

```bash
~/mpm/build_amrex/amrex_smoke_test
mpirun --oversubscribe -np 2 ~/mpm/build_amrex/amrex_smoke_test
```

### 6.4 运行 MPM

```bash
~/mpm/build_amrex/mpm -f ~/mpm/demo/ -p 8 -i mpm.json
```

打开 AMReX 后，日志应包含：

```text
AMReX backend enabled: version 26.11
```

## 7. 后续怎样在 MPM 中真正使用 AMReX

### 7.1 不要直接把 AMReX 类型放进现有 Mesh

不要直接在 `mesh.h`、`cell.h` 或 `MPMScheme` 中加入：

```cpp
amrex::BoxArray
amrex::MultiFab
amrex::ParticleContainer
```

这些类型会把整个 MPM 编译链拖入 C++20，并让原网格代码与 AMReX 强耦合。

推荐新增一层不依赖 AMReX 的后端接口：

```text
include/mesh_adaptivity/amr_backend.h
src/mesh_adaptivity/amrex_backend.cc
```

`amr_backend.h` 只暴露 MPM 自己的类型：

```text
LevelId
CellIndex
PatchId
RefineBox
ParticleTransfer
```

### 7.2 AMReX 对象与 MPM 数据的对应关系

建议的映射：

```text
amrex::Box                  -> 一个规则 Hex8 逻辑块
amrex::BoxArray             -> 某一级的所有 patch
amrex::DistributionMapping  -> patch 的 MPI 归属
amrex::MultiFab             -> 单元中心辅助场
amrex::AmrMesh / AmrCore    -> 多级网格和重网格
amrex::TagBox               -> 需要加密的单元标记
amrex::ParticleContainer    -> 粒子的存储和 MPI 迁移
amrex::FillPatch            -> 粗层给细层提供 ghost 数据
amrex::FluxRegister         -> 粗细层通量修正参考
```

注意：MPM 的应力、变形梯度和材料历史仍应由现有 `ParticleBase` 保存。
AMReX 不提供 MPM 材料模型、Piola 应力更新或粒子到节点的映射。

### 7.3 第一个可落地目标

建议先实现：

```text
3D
两级网格
固定局部加密区
正交 Hex8
标准 MPM 单元
先不动态重网格
```

第一阶段只验证：

- AMReX 的两级 patch 能映射成 MPM `Cell`。
- 粒子能被分配到正确层级。
- 节点质量和动量在粗细界面不会重复累加。
- 应力波穿过界面时没有明显非物理反射。

### 7.4 第二阶段再增加动态功能

第二阶段再处理：

- 误差标记和 `TagBox`
- 每隔固定步数 regrid
- 粒子跨层迁移
- 粒子分裂/合并
- 粗细层时间子循环
- FluxRegister 或 MPM 专用的粗细界面修正

### 7.5 新增 AMReX API 时的安全规则

如果需要在 bridge 中增加一个方法：

1. 在 `amrex_runtime.h` 中只声明 C++14 可编译的接口。
2. 接口参数优先使用 `double`、`int`、`std::array` 或 MPM 自己的数据结构。
3. 在 `amrex_runtime.cc` 或其他 C++20 backend `.cc` 中包含 AMReX 头文件。
4. 不在 `mesh.h`、`cell.h`、`particle_base.h` 中加入 `amrex::` 类型。
5. 每增加一种跨层操作，就增加对应 smoke test 和数值验证。

## 8. 当前验证结果

已经完成：

- `MPM_WITH_AMREX=OFF` 独立构建成功。
- `MPM_WITH_AMREX=ON` 独立构建成功。
- `amrex_smoke_test` 单 MPI rank 返回 0。
- `mpirun --oversubscribe -np 2 amrex_smoke_test` 返回 0。
- AMReX 版 MPM 单步 demo 返回 0。
- 现有均匀 `mesh.txt` 求解流程未修改。

## 9. 常见问题

### `find_package(AMReX) failed`

检查：

```bash
ls ~/opt/amrex-26.11/lib/cmake/AMReX/AMReXConfig.cmake
```

配置时加入：

```bash
-DCMAKE_PREFIX_PATH="$HOME/opt/amrex-26.11"
```

### 默认 `g++` 版本太低

使用：

```bash
-DCMAKE_CXX_COMPILER=g++-13
```

### AMReX 编译到 `-O3` 时崩溃

使用当前安装脚本默认的 `-O1`，或换用更稳定的 GCC/Clang 工具链。

### 运行时出现 `Authorization required`

这是虚拟机图形显示环境的 X11 提示，不是 AMReX 求解失败。程序返回 0 时
可以忽略。
