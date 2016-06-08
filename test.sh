#!/bin/bash

PASS=0;
FAIL=0;
ERROR=0;

C=$1
let NP="$C*2"

make clean && make
if [[ $? == 0 ]] ; then echo ""; else exit ; fi

rm -f fails/*.txt

for Y in 10 64 ; do
  for X in `seq 3` ; do
    for I in `seq 0 9` ; do
      Z=00$I
      A=$(ls -1 exported_tests/matrix* | sed -E "s/(.*)_(.*)_(.*)/\3/g" | head -n $(let Z1="$I+1" ; echo $Z1) | tail -n 1)

      INFILE=exported_tests/sparse05_000${Y}_$Z
      OUTFILE=exported_tests/result_${X}_000${Y}_$Z
      #BFILE=exported_tests/matrix01_000${Y}_$A

      echo -n "test $INFILE $X... "

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
  done
done

echo "PASS $PASS, FAIL $FAIL, ERROR $ERROR"
