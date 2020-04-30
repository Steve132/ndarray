cmake_minimum_required(VERSION 3.10)


add_library(ndarray INTERFACE)
target_compile_features(ndarray INTERFACE cxx_std_11)
target_include_directories(ndarray INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(nd_ext)
add_executable(ndarray_test
	ndarray_test.cpp
	nd/ndarray.hpp
)
target_link_libraries(ndarray_test ndarray)
