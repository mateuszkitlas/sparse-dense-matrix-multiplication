#!/bin/bash

C=2
let NP="$C*2"
E=3
SEED=10
FILE=exported_tests/sparse05_00010_000

cat $FILE
echo "NP=$NP C=$C E=$E SEED=$SEED FILE=$FILE"

make clean && make ARGS="-DDEBUG" && mpirun -np $NP ./matrixmul -f $FILE -s $SEED -c $C -e $E
