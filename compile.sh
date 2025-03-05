#!/bin/bash

FC=ifort

FOPT="-O0 -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn all -warn interfaces -fp-stack-check"

rm -f *.o

# MODULES: compile to an object (.o) file only

# MAIN: compile to an object (.o) file only
$FC -c $FOPT co2_adsorption_generator.f90

# executable file
$FC -o co2_adsorption_generator.x $FOPT *.o

exit 0
