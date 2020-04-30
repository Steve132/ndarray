cmake_minimum_required(VERSION 3.10)

add_library(nd_image
	image.hpp
	image.cpp
)

target_compile_options(nd_image PRIVATE "-march=native")
set_target_properties(nd_image PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

target_include_directories(nd_image PUBLIC $(CMAKE_CURRENT_SOURCE_DIR)/../../)
