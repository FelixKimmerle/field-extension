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
