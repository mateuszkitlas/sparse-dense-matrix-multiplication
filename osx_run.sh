#!/bin/bash

NP=4
FILE=exported_tests/sparse05_00010_000

cat $FILE

f(){ ./matrixmul -f $FILE -s seed_for_dense_matrix -c repl_group_size -e exponent ; }
#make \
#  && mpirun -np $NP ./matrixmul -f exorted_tests/matrix01_00010_00042 -s seed_for_dense_matrix -c repl_group_size -e exponent

make clean && make ARGS="-DDONT_USE_MPI -DDEBUG" && f
