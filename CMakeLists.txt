cmake_minimum_required(VERSION 3.10)

add_executable(ndarray_test
	ndarray_test.cpp
	nd/ndarray.hpp
)
target_compile_options(ndarray_test PRIVATE "-march=native")
set_target_properties(ndarray_test PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

add_subdirectory(nd_ext)
