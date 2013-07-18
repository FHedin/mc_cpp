#!/bin/bash

#SCOREP_DIR=/opt/scorep/bin
#scorep --nompi --noopenmp --nocompiler --user g++ -O3 -DNDEBUG -std=c++11 -I ../include -I ../rapidxml-1.13 -I . ../src/*.cpp -o mc_cpp_scorep

export PATH=$PATH:$SCOREP_DIR

#export SCOREP_EXPERIMENT_DIRECTORY=scorep_serial_tracing
export SCOREP_EXPERIMENT_DIRECTORY=scorep_serial_profiling
#export SCOREP_FILTERING_FILE=scorep.filt
#export SCOREP_METRIC_PAPI=PAPI_L1_DCM,PAPI_L2_DCM
export SCOREP_ENABLE_TRACING=true
export SCOREP_ENABLE_PROFILING=true
export SCOREP_TOTAL_MEMORY=128M

../build/mc_cpp_scorep -i input.xml > out.log


