cmake_minimum_required(VERSION 3.10)

add_library(nd_image
	image.hpp
	image.cpp
)
target_compile_features(nd_image INTERFACE cxx_std_11)

target_link_libraries(nd_image ndarray)
