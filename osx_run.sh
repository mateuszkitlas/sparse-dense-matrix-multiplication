#!/bin/bash

NP=4
FILE=exported_tests/sparse05_00010_000

cat $FILE

make ARGS="-DDEBUG" && mpirun -np $NP ./matrixmul -f $FILE -s seed_for_dense_matrix -c repl_group_size -e exponent
