cmake_minimum_required (VERSION 2.6)

set(TGT "mc_cpp" CACHE STRING "The executable name.")

project(${TGT} CXX)

# The version numbers, they are added to a dedicated header "Global_include.h" present everywhere .
set(VERSION_MAJOR 0 CACHE STRING "Major version number of this program.")
set(VERSION_MINOR 2 CACHE STRING "Minor version number of this program.")

configure_file (
  "include/Global_include.h.in"
  "${PROJECT_BINARY_DIR}/Global_include.h"
  )
  
# list of include folders with .h / .hpp files
include_directories(
  include 
  ${PROJECT_BINARY_DIR}
)

# list all source files
# WARNING : for each new source file added, it is necessary to run cmake again !
# TODO : For release versions, this has to be replaced by a hand-written list.
aux_source_directory(src SRCS)

# add custom options for compiler
set(CUSTOM_CXX_FLAGS "-std=c++11 -Wall -Wextra -Wno-unused-but-set-variable" CACHE STRING "Custom c++ compiler options.")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CUSTOM_CXX_FLAGS}")

add_executable(${TGT} ${SRCS})
