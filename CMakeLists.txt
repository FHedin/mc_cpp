cmake_minimum_required (VERSION 2.6)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug RelWithDebInfo MinSizeRel Release."
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
include/experimental
rapidxml-1.13
vectorclass-1.16
${PROJECT_BINARY_DIR}
)

# list all source files
set(SRCS
src/Angle.cpp
src/AtomList.cpp
src/Bond.cpp      
src/Bond_UB.cpp   
src/Dihedral.cpp  
src/Dihedral_improper.cpp  
src/Ensemble.cpp           
src/Ens_NVT.cpp            
src/FField.cpp             
src/FField_MDBAS.cpp    
src/IO.cpp              
src/IO_MDBAS.cpp        
src/List_Moves.cpp      
src/List_nonBonded.cpp  
src/mc_cpp.cpp           
src/MC.cpp             
src/MC_metropolis.cpp  
src/Parser.cpp
src/PerConditions.cpp
src/Selection.cpp
src/experimental/Ens_NPT.cpp
src/experimental/Ens_NPT.cpp
src/experimental/IO_CHARMM.cpp
src/experimental/MC_spav.cpp
)

#uncomment each line for enabling experimental features (not well supported for the moment)
add_definitions(
# -D RANGED_E_EXPERIMENTAL
# -D BALDRICH_EXPERIMENTAL
# -D NPT_EXPERIMENTAL
# -D SPAV_EXPERIMENTAL
# -D CHARMM_EXPERIMENTAL
# -D VECTORIZED_ENER_EXPERIMENTAL
-D VECTORCLASS_EXPERIMENTAL
)

#default
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG           "-mavx2 -O0 -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-mavx2 -O2 -g -DNDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL      "-mavx2 -Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE         "-mavx2 -O3 -DNDEBUG")

# Compiler-specific C++11 and OpenMP activation.

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if ( (GCC_VERSION VERSION_GREATER 4.8) OR (GCC_VERSION VERSION_EQUAL 4.8) )
#         set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -fopenmp")
        set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra")
    else ()
#         set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -Wextra -fopenmp")
        set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -Wextra")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    #openmp not available with clang !!	
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
# 	set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -openmp")
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra")
endif ()

add_executable(${TGT} ${SRCS})
	
target_link_libraries(${TGT} m)



