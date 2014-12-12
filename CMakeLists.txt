cmake_minimum_required (VERSION 2.6)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(TGT "mc_cpp" CACHE STRING "The executable name.")

project(${TGT} CXX)

# The version numbers, they are added to a dedicated header "Global_include.h" present everywhere .
set(VERSION_MAJOR 0 CACHE STRING "Major version number of this program.")
set(VERSION_MINOR 2 CACHE STRING "Minor version number of this program.")

configure_file (
  "include/Global_include.hpp.in"
  "${PROJECT_BINARY_DIR}/Global_include.hpp"
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

# add_definitions(-D...)

set(CMAKE_CXX_FLAGS_DEBUG   "-std=c++11 -Wall -Wextra -Wno-unused-but-set-variable -O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -Wall -Wextra -Wno-unused-but-set-variable -O3 -march=native")

add_executable(${TGT} ${SRCS})
target_link_libraries(${TGT} m)




