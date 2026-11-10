# C++20 与 AMReX 直接接入修改说明

本文对应分支 `C++20_AMReX`，修改日期为 2026-11-06。行号均为修改后行号。

## 1. 修改总览

本次修改完成两件事：

1. 主工程从 C++14 迁移到 C++20。
2. 删除上一版为 C++14 主程序保留的独立 C++20 AMReX bridge。

修改后的结构：

```text
mpm 主程序       C++20
AmrexRuntime     C++20
AMReX::amrex     C++20
```

修改前为：

```text
mpm 主程序 C++14 -> C++20 bridge 静态库 -> AMReX
```

## 2. CMakeLists.txt

### 2.1 第 5-7 行：统一切换到 C++20

修改前：

```cmake
set(CMAKE_CXX_STANDARD 14)
```

修改后：

```cmake
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
```

逐行说明：

- 第 5 行：所有 C++ target 使用 C++20。
- 第 6 行：编译器不支持 C++20 时配置直接失败，不允许静默降级。
- 第 7 行：关闭 GNU 扩展，实际参数为严格的 `-std=c++20`。

早期探测中严格 C++20 曾触发 GCC 内部优化错误；删除 `-fpermissive` 后，
严格 `-std=c++20` 已验证可以构建。

### 2.2 删除 Intel 分支中的 `-std=c++14`

原代码在 Intel 编译器分支中强制：

```cmake
-std=c++14
```

该分支已删除。语言标准现在只由 `CMAKE_CXX_STANDARD` 控制，不能再被
编译器专用 flags 覆盖。

### 2.3 删除 GNU 分支中的 `-fpermissive`

原工程对所有 GNU 编译无条件加入 `-fpermissive`。C++20 迁移后该选项会继续
放宽本应报错的旧语法，并可能干扰优化器。该设置已删除，代码必须在正常
编译规则下通过。

### 2.4 第 221 行：AMReX runtime 加入主程序源文件

```cmake
${mpm_SOURCE_DIR}/src/mesh_adaptivity/amrex_runtime.cc
```

原实现将该文件编入独立静态库。现在它直接属于 `mpm`，不再需要跨 C++ 标准
使用静态库连接。

### 2.5 第 239-252 行：简化 AMReX 链接

修改后：

```cmake
if(MPM_WITH_AMREX)
  target_compile_definitions(mpm PRIVATE MPM_WITH_AMREX=1)
  target_link_libraries(mpm PRIVATE AMReX::amrex)

  add_executable(amrex_smoke_test ...)
  target_compile_definitions(amrex_smoke_test PRIVATE MPM_WITH_AMREX=1)
  target_link_libraries(amrex_smoke_test PRIVATE AMReX::amrex)
endif()
```

逐行说明：

- 第 239 行：只有启用 AMReX 时进入配置。
- 第 240 行：主程序启用 AMReX runtime 调用。
- 第 241 行：主程序直接链接 `AMReX::amrex`。
- 第 243-246 行：创建 smoke test，并直接编译 runtime 源文件。
- 第 247-249 行：smoke test 使用 MPM 的 include 目录。
- 第 250 行：smoke test 启用 AMReX 编译宏。
- 第 251 行：smoke test 直接链接 AMReX。

删除的内容：

- `mpm_amrex_runtime` 静态库。
- bridge 专用的 `cxx_std_20` 设置。
- 主程序通过中间静态库链接 AMReX 的步骤。

## 3. C++20 语法兼容修改

### 3.1 `include/containers/map.h` 第 24 行

修改前：

```cpp
Map<T>() = default;
```

修改后：

```cpp
Map() = default;
```

C++20 不允许在类内构造函数声明中重复写模板实参。

### 3.2 `include/containers/vector.h` 第 18 行

修改前：

```cpp
Vector<T>() = default;
```

修改后：

```cpp
Vector() = default;
```

原因与 `Map` 相同。

### 3.3 `external/tclap/ValueArg.h` 第 243 行

修改前：

```cpp
ValueArg<T>(const ValueArg<T>& rhs);
```

修改后：

```cpp
ValueArg(const ValueArg<T>& rhs);
```

这是禁用拷贝的构造函数声明，C++20 不再接受旧式模板构造名。

### 3.4 `external/tclap/MultiArg.h` 第 228 行

修改前：

```cpp
MultiArg<T>(const MultiArg<T>& rhs);
```

修改后：

```cpp
MultiArg(const MultiArg<T>& rhs);
```

原因与 `ValueArg` 相同。

### 3.5 `include/mesh/mesh.tcc` 第 1103 行

修改前：

```cpp
[=, &particles](...)
```

修改后：

```cpp
[this, &particles](...)
```

C++20 不再推荐通过 `[=]` 隐式捕获 `this`。新写法明确捕获当前 `Mesh` 对象。

## 4. AMReX runtime 简化

### 4.1 `include/mesh_adaptivity/amrex_runtime.h` 第 4-30 行

修改前使用 Pimpl：

```cpp
#include <memory>

class Impl;
std::unique_ptr<Impl> impl_;
```

修改后只前置声明 AMReX 类型并直接保存指针：

```cpp
namespace amrex {
class AMReX;
}

class AmrexRuntime {
 private:
  amrex::AMReX* instance_{nullptr};
};
```

原因：

- 主工程已经是 C++20，不再需要 Pimpl 隔离编译标准。
- 仍然不把完整 `AMReX.H` 放入公共头文件。
- 删除动态分配的 `Impl`，runtime 对象更简单。

### 4.2 `src/mesh_adaptivity/amrex_runtime.cc` 第 12-49 行

核心修改：

```cpp
AmrexRuntime::AmrexRuntime(...) : instance_(nullptr) {
  instance_ = amrex::Initialize(...);
}

AmrexRuntime::~AmrexRuntime() {
  if (instance_ != nullptr) amrex::Finalize(instance_);
}
```

逐行说明：

- 第 12-13 行：直接初始化 AMReX 实例指针。
- 第 14-20 行：按 AMReX 是否启用 MPI 选择初始化重载。
- 第 21-26 行：未启用 AMReX 时保留明确异常。
- 第 29-33 行：析构时关闭 AMReX。
- 第 35-41 行：直接根据 `instance_` 判断状态。
- 第 43-49 行：版本查询逻辑不变。

## 5. src/main.cc

本次没有改变 `main.cc` 的逻辑。它仍然：

- 在 AMReX 构建时创建 `AmrexRuntime`；
- runtime 生命周期覆盖 MPM solve；
- AMReX 在 `MPI_Finalize` 前析构；
- 日志输出 AMReX 版本。

区别只是 runtime 实现现在直接属于主程序，不再经过独立静态库。

## 6. tests/mesh_adaptivity/amrex_smoke_test.cc

测试逻辑不变，仍然：

- 包含 `AMReX.H` 和 `AMReX_Box.H`；
- 创建 `amrex::Box`；
- 检查 runtime 状态。

CMake 现在把 runtime 源文件和测试源文件直接编译到 smoke test。

## 7. .vscode/tasks.json

### 第 7 行和第 40 行

Release/Debug 构建原来使用：

```text
-DCMAKE_CXX_COMPILER=g++
```

修改为：

```text
-DCMAKE_CXX_COMPILER=g++-13
```

原因：Ubuntu 默认 `g++` 仍是 9.5，不能可靠支持当前 C++20 和
AMReX 26.11。

## 8. 接入原理

修改后的依赖关系：

```text
MPM C++20
  |
  +-- AmrexRuntime
        |
        +-- AMReX::amrex
```

后续 AMR 后端可以直接使用：

```cpp
amrex::BoxArray
amrex::DistributionMapping
amrex::MultiFab
amrex::ParticleContainer
```

仍建议将复杂 AMReX 逻辑放在 `src/mesh_adaptivity/` 的 backend 中，不要直接
写进 `Mesh` 或 `MPMScheme`。这样网格后端仍可替换和独立测试。

## 9. 编译方式

默认 C++20：

```bash
cmake -S ~/mpm -B ~/mpm/build_cxx20 \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DMPM_WITH_AMREX=OFF \
  -DNO_KAHIP=True
cmake --build ~/mpm/build_cxx20 --parallel 8
```

AMReX C++20：

```bash
cmake -S ~/mpm -B ~/mpm/build_cxx20_amrex \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_PREFIX_PATH="$HOME/opt/amrex-26.11" \
  -DMPM_WITH_AMREX=ON \
  -DNO_KAHIP=True
cmake --build ~/mpm/build_cxx20_amrex --parallel 8
```

## 10. 验证结果

已完成：

- `MPM_WITH_AMREX=OFF` 的 C++20 Debug 构建成功。
- `MPM_WITH_AMREX=OFF` 的 C++20 标准 Release 构建成功。
- `MPM_WITH_AMREX=ON` 的 C++20 Debug 和 Release 构建成功。
- `amrex_smoke_test` 单 MPI rank 返回 0。
- `mpirun --oversubscribe -np 2 amrex_smoke_test` 返回 0。
- C++20 + AMReX 单步 demo 返回 0。
- 完整测试目标可以编译。

测试运行情况：

- `mpmtest_unit`：110 个测试，109 通过，1 失败。
- 失败项为 `tests/mesh/mesh_test_3d.cc:1043` 的 `radius_factor` 断言。
- 该测试前面故意写入了非法粒子单元 ID，后面却要求所有
  `radius_factor` 都是有限正数，属于既有测试逻辑问题。
- `mpmtest_small`：53 个测试，33 通过，20 失败；失败均集中在
  `checkpoint_resume()`，需要后续单独核查测试数据/工作目录。

这些失败不是本次 C++20 语法迁移或 AMReX 链接造成的编译错误。

## 11. 审查清单

1. `CMakeLists.txt` 第 5-7 行是否全局使用 C++20。
2. `mpm_amrex_runtime` 中间库是否已经完全删除。
3. `Map`、`Vector`、TCLAP 的构造语法是否已修正。
4. `AmrexRuntime` 是否保证 `Finalize` 早于 `MPI_Finalize`。
5. `MPM_WITH_AMREX=OFF` 是否可独立构建。
6. `MPM_WITH_AMREX=ON` 是否可构建且 smoke test 返回 0。
