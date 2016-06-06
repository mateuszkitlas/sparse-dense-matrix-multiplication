#!/bin/bash

NP=4
C=3
E=3
SEED=10
FILE=exported_tests/sparse05_00010_000

cat $FILE

make clean && make ARGS="-DDEBUG" && mpirun -np $NP ./matrixmul -f $FILE -s $SEED -c $C -e $E
