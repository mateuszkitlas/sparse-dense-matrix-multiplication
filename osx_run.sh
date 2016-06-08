#!/bin/bash

#make clean && make ARGS="-DDEBUG -DIDENTITY_MATRIX"
#make clean && make ARGS="-DDEBUG"
make clean && make
if [[ $? == 0 ]] ; then echo ""; else exit ; fi


C=2
let NP="$C*2"


X=1
Z=000
Y=10


A=$(ls -1 exported_tests/matrix* | sed -E "s/(.*)_(.*)_(.*)/\3/g" | head -n $(let Z1="$Z+1" ; echo $Z1) | tail -n 1)


INFILE=exported_tests/sparse05_000${Y}_$Z
OUTFILE=exported_tests/result_${X}_000${Y}_$Z
BFILE=exported_tests/matrix01_000${Y}_$A

#cat $INFILE

diff <(mpirun -np $NP ./matrixmul -f $INFILE -s $A -c $C -e $X -v) $OUTFILE
#cat $BFILE
