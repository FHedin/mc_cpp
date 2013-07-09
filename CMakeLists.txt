cmake_minimum_required (VERSION 2.6)

set(TGT "mc_cpp" CACHE STRING "The executable name.")

project(${TGT} CXX)

# The version number.
set(${TGT}_VERSION_MAJOR 0)
set(${TGT}_VERSION_MINOR 2)

# list of include folders with .h / .hpp files
include_directories(
  include
)

# list all source files
# WARNING : for each new source file added, it is necessary to run cmake again !
# TODO : For release versions, this has to be replaced by a hand-written list.
aux_source_directory(src SRCS)

# add custom options for compiler
set(CUSTOM_CXX_FLAGS "-std=c++11" CACHE STRING "Custom c++ compiler options.")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CUSTOM_CXX_FLAGS}")

add_executable(${TGT} ${SRCS})