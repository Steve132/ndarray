cmake_minimum_required(VERSION 3.10)

add_library(nd_image
	image.hpp
	image.cpp
)

target_compile_options(nd_image PRIVATE "-march=native")
target_link_libraries(nd_image ndarray)
