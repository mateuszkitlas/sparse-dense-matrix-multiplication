#!/bin/bash

PASS=0;
FAIL=0;
ERROR=0;

C=$1
let NP="$C*$2"

#make clean && 
make
if [[ $? == 0 ]] ; then echo ""; else exit ; fi

rm -f fails/*.txt

for OUTFILE in `ls exported_tests/result*` ; do
  X=$(basename $OUTFILE | sed -E "s/result_(.*)_(.*)_(.*)_(.*)/\1/g")
  Y=$(basename $OUTFILE | sed -E "s/result_(.*)_(.*)_(.*)_(.*)/\2/g")
  Z=$(basename $OUTFILE | sed -E "s/result_(.*)_(.*)_(.*)_(.*)/\3/g")
  A=$(basename $OUTFILE | sed -E "s/result_(.*)_(.*)_(.*)_(.*)/\4/g")

  INFILE=exported_tests/sparse05_${Y}_$Z
  #OUTFILE=exported_tests/result_${X}_${Y}_${Z}_$A

  echo -n "test $INFILE $X... "

  mpirun -np $NP ./matrixmul -f $INFILE -s $A -c $C -e $X -ge 1
  mpirun -np $NP ./matrixmul -f $INFILE -s $A -c $C -e $X -v 2>/dev/null > out.txt
  if [[ $? == 0 ]] ; then
    python diff_numbers.py out.txt $OUTFILE

    if [[ $? == 0 ]] ; then
      let PASS="$PASS+1"
      echo OK
    else
      let FAIL="$FAIL+1"
      echo FAIL
      diff -w out.txt $OUTFILE >/dev/null >fails/Y${Y}_X${X}_Z${Z}_C${C}.txt
    fi
  else
    let ERROR="$ERROR+1"
    echo ERROR
  fi
done

echo "PASS $PASS, FAIL $FAIL, ERROR $ERROR"
