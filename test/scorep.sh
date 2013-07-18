#!/bin/bash

SCOREP_DIR=/opt/scorep/bin

export PATH=$PATH:$SCOREP_DIR

export SCOREP_EXPERIMENT_DIRECTORY=scorep_serial_tracing
export SCOREP_FILTERING_FILE=scorep.filt
#export SCOREP_METRIC_PAPI=PAPI_FP_OPS
export SCOREP_ENABLE_TRACING=true
export SCOREP_ENABLE_PROFILING=true
export SCOREP_TOTAL_MEMORY=2G

../build/mc_cpp_scorep -i input.xml > out.log


