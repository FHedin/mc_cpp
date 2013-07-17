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
  rapidxml-1.13 
  ${PROJECT_BINARY_DIR}
)

# list all source files
set(SRCS
src/Angle.cpp     
src/Dihedral_improper.cpp  
src/FField_MDBAS.cpp    
src/main.cpp           
src/Atom.cpp      
src/Ensemble.cpp           
src/IO.cpp              
src/MC.cpp             
src/Parser.cpp
src/Bond.cpp      
src/Ens_NVT.cpp            
src/IO_MDBAS.cpp        
src/MC_metropolis.cpp  
src/PerConditions.cpp
src/Bond_UB.cpp   
src/List_Moves.cpp      
src/Selection.cpp
src/Dihedral.cpp  
src/FField.cpp             
src/List_nonBonded.cpp  
)

#Choose optimisation level
#set(CUSTOM_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose here the optimisation level.")
#set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE} ${CUSTOM_BUILD_TYPE}")

# add custom options for compiler
set(CUSTOM_CXX_FLAGS "-std=c++11 -Wall -Wextra -Wno-unused-but-set-variable" CACHE STRING "Custom c++ compiler options.")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CUSTOM_CXX_FLAGS}")

# common libraries
link_libraries(m)

add_executable(${TGT} ${SRCS})


