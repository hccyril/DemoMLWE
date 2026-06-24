#==================================================================================
# CMakeLists.txt
#
# 基于 OpenFHE 的 MLWE 模块与演示程序的构建脚本。
#
# 依赖：
#   - OpenFHE (v1.x) 已安装，并能通过 find_package(OpenFHE) 找到。
#     安装方式参见 https://github.com/openfheorg/openfhe-development
#
# 使用：
#   mkdir build && cd build
#   cmake .. -DOpenFHE_DIR=<OpenFHE 安装路径>/lib/cmake/OpenFHE
#   cmake --build . --config Release
#==================================================================================

cmake_minimum_required(VERSION 3.16)

project(OpenFHE_MLWE_Demo
    VERSION 1.0.0
    DESCRIPTION "Module-LWE (MLWE) module and demo based on OpenFHE"
    LANGUAGES CXX)

# ---------- C++ 标准 ----------
# OpenFHE 要求 C++17。
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 默认构建类型为 Release（性能更优，OpenFHE 的 NTT 等会更快）
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------- 查找 OpenFHE ----------
# find_package 会自动设置 OpenFHE::openfhe 这个 imported target。
find_package(OpenFHE REQUIRED)

# ---------- 包含目录 ----------
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${OpenFHE_INCLUDE}
)

# ---------- 编译选项 ----------
# 与 OpenFHE 推荐一致：开启优化、O3、并行（如 OpenMP 可用）。
if(MSVC)
    add_compile_options(/W3 /bigobj /EHsc)
    add_compile_definitions(_USE_MATH_DEFINES NOMINMAX)
else()
    add_compile_options(-Wall -Wextra -Wno-unused-parameter -O3)
endif()

#==================================================================================
# 目标 1：MLWE 模块静态库
#==================================================================================
add_library(mlwe STATIC
    src/mlwe.cpp
)
target_include_directories(mlwe PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)
# 链接 OpenFHE
target_link_libraries(mlwe PUBLIC OpenFHE::openfhe)

#==================================================================================
# 目标 2：演示可执行文件
#==================================================================================
add_executable(mlwe_demo
    src/mlwe_demo.cpp
)
target_link_libraries(mlwe_demo PRIVATE mlwe)

#==================================================================================
# 安装规则（可选）
#==================================================================================
include(GNUInstallDirs)
install(TARGETS mlwe mlwe_demo
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
