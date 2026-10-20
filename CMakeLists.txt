cmake_minimum_required(VERSION 3.20)
project(fields LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# header only library
add_library(fields INTERFACE)
target_include_directories(fields INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(fields INTERFACE cxx_std_20)

add_executable(demo examples/main.cpp)
target_link_libraries(demo PRIVATE fields)
target_compile_options(demo PRIVATE -Wall -Wextra -Wpedantic)

enable_testing()
add_executable(tests tests/tests.cpp)
target_link_libraries(tests PRIVATE fields)
target_compile_options(tests PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME tests COMMAND tests)

# optional: Eigen adapter (include/fields/eigen.hpp), Eigen is downloaded with FetchContent
option(FIELDS_WITH_EIGEN "Fetch Eigen and build the Eigen adapter example and tests" ON)

if(FIELDS_WITH_EIGEN)
    include(FetchContent)

    # only the headers are needed, keep Eigen's own tests, docs and BLAS/LAPACK out of the build
    set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_BLAS OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_LAPACK OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_CMAKE_PACKAGE OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        Eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        SYSTEM
    )
    FetchContent_MakeAvailable(Eigen)

    add_executable(linear_algebra examples/linear_algebra.cpp)
    target_link_libraries(linear_algebra PRIVATE fields Eigen3::Eigen)
    target_compile_options(linear_algebra PRIVATE -Wall -Wextra -Wpedantic)

    add_executable(eigen_tests tests/eigen_tests.cpp)
    target_link_libraries(eigen_tests PRIVATE fields Eigen3::Eigen)
    target_compile_options(eigen_tests PRIVATE -Wall -Wextra -Wpedantic)
    add_test(NAME eigen_tests COMMAND eigen_tests)
endif()
