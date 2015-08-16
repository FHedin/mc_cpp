cmake_minimum_required (VERSION 2.8)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug RelWithDebInfo MinSizeRel Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(TGT "mc_cpp" CACHE STRING "The executable name.")

project(${TGT} CXX)

# set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
# find_package(OpenCL)
# message(STATUS "OpenCL found: ${OPENCL_FOUND}")
# message(STATUS "OpenCL includes: ${OPENCL_INCLUDE_DIRS}")
# message(STATUS "OpenCL CXX includes: ${OPENCL_HAS_CPP_BINDINGS}")
# message(STATUS "OpenCL libraries: ${OPENCL_LIBRARIES}")

# The version numbers, they are added to a dedicated header "Global_include.h" present everywhere .
set(VERSION_MAJOR 0 CACHE STRING "Major version number of this program.")
set(VERSION_MINOR 3 CACHE STRING "Minor version number of this program.")

configure_file (
"include/Global_include.hpp.in"
"${PROJECT_BINARY_DIR}/Global_include.hpp"
)
  
# list of include folders with .h / .hpp files
include_directories(
${PROJECT_BINARY_DIR}
include
include/experimental
rapidxml-1.13
vectorclass-1.17
# uncomment the 2 following lines for enabling OpenCL
# CL
# ${OPENCL_INCLUDE_DIRS}
)

# list all source files for executable
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
src/FField_MDBAS_VECT.cpp
src/IO.cpp
src/IO_MDBAS.cpp
src/List_Moves.cpp
src/List_nonBonded.cpp
src/MC.cpp
src/mc_cpp.cpp
src/MC_metropolis.cpp
src/Parser.cpp
src/PerConditions.cpp
src/Selection.cpp
src/experimental/Ens_NPT.cpp
src/experimental/Ens_NPT.cpp
src/experimental/IO_CHARMM.cpp
src/experimental/MC_spav.cpp
src/experimental/FField_MDBAS_CL.cpp
)

# list source files for library
set(LIBSRCS
src/FORTRAN_interface
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
src/FField_MDBAS_VECT.cpp
src/IO.cpp
src/IO_MDBAS.cpp
src/List_Moves.cpp
src/List_nonBonded.cpp
src/PerConditions.cpp
src/Selection.cpp
)

# uncomment each line for enabling experimental features
# only VECTORCLASS_EXPERIMENTAL recommended, will provide a good speed up on SSE/SSE2/SSE3/AVX capable cpus (almost all cpus compatible nowadays 
#     so only disable if you get errors while compiling)
add_definitions(
# -D RANGED_E_EXPERIMENTAL
# -D BALDRICH_EXPERIMENTAL
# -D NPT_EXPERIMENTAL
# -D SPAV_EXPERIMENTAL
# -D CHARMM_EXPERIMENTAL
# -D VECTORCLASS_EXPERIMENTAL
# -D OPENCL_EXPERIMENTAL
)

#default
set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG           "-O0 -g ")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-O2 -g -DNDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL      "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE         "-O3 -DNDEBUG")

# Compiler-specific C++11 and OpenMP activation.

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if ( (GCC_VERSION VERSION_GREATER 4.8) OR (GCC_VERSION VERSION_EQUAL 4.8) )
        set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -fabi-version=6 -Wabi -Wall -Wextra -fopenmp")
    else ()
        set(CMAKE_CXX_FLAGS "-std=c++0x -march=native -fabi-version=6 -Wabi -Wall -Wextra -fopenmp")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    #openmp not available with clang !!	
    set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -Wall -Wextra")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "-std=c++11 -march=native -Wall -Wextra -openmp")
endif ()

add_executable(${TGT} ${SRCS})
target_link_libraries(${TGT} m)

add_library(charmmene STATIC ${LIBSRCS})

# uncomment the following line for enabling OpenCL
# target_link_libraries(${TGT} m ${OPENCL_LIBRARIES})


