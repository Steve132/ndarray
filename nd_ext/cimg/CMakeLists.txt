cmake_minimum_required(VERSION 3.10)

add_library(nd_cimg
	display.hpp
	display.cpp
)

target_compile_options(nd_cimg PRIVATE "-march=native")
target_link_libraries(nd_cimg PUBLIC ndarray)
target_link_libraries(nd_cimg PRIVATE X11)
target_compile_features(nd_cimg INTERFACE cxx_std_11)
target_include_directories(nd_cimg PRIVATE deps)
