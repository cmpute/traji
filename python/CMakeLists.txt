cmake_minimum_required(VERSION 3.12)

project(traji_cython)

find_package(PythonExtensions REQUIRED)
find_package(Cython REQUIRED)
find_package(NumPy REQUIRED)
include_directories(${NumPy_INCLUDE_DIRS})

get_filename_component(main_dir ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
find_package(traji REQUIRED HINTS ${main_dir}/build)

add_subdirectory(traji)
