cmake_minimum_required(VERSION 2.6)

project(MQM)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules/")

add_subdirectory(src)

enable_testing()

