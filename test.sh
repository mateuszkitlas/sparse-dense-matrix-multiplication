#!/bin/bash

PASS=0;
FAIL=0;
ERROR=0;

C=$1
let NP="$C*$2"
let USE_IDENTITY_ID="$3+0"
ARGS=""
if [[ $USE_IDENTITY == 1 ]] ; then
  ARGS="-DIDENTITY_MATRIX"
fi
let ONLY_FIRST="$4+0"

diff make_args <(echo $ARGS)
if [[ $? != 0 ]] ; then make clean ; fi
echo $ARGS > make_args

make ARGS=$ARGS
if [[ $? == 0 ]] ; then echo ""; else exit ; fi

rm -f fails/*.txt

for OUTFILE_PATH in `ls exported_tests/result*` ; do
  for OPTIONS in "" "-i" ; do
    OUTFILE=$(basename $OUTFILE_PATH)
    X=$(echo $OUTFILE | sed -E "s/result_(.*)_(.*)_(.*)_(.*)/\1/g")
    Y=$(echo $OUTFILE | sed -E "s/result_(.*)_(.*)_(.*)_(.*)/\2/g")
    Z=$(echo $OUTFILE | sed -E "s/result_(.*)_(.*)_(.*)_(.*)/\3/g")
    A=$(echo $OUTFILE | sed -E "s/result_(.*)_(.*)_(.*)_(.*)/\4/g")

    INFILE=exported_tests/sparse05_${Y}_$Z
    CMD="mpirun -np $NP ./matrixmul -f $INFILE -s $A -c $C -e $X"

    echo -n "test $INFILE $X... "

    if [[ $ONLY_FORST == 1 ]] ; then
      $CMD -v $OPTIONS
    else
      $CMD -v $OPTIONS 2>/dev/null > out.txt
    fi

    if [[ $? == 0 ]] ; then
      python diff_numbers.py out.txt $OUTFILE_PATH

      if [[ $? == 0 ]] ; then
        let PASS="$PASS+1"
        echo OK
        $CMD -ge 1 >/dev/null
      else
        let FAIL="$FAIL+1"
        echo FAIL
        diff -w out.txt $OUTFILE >/dev/null >fails/$OUTFILE
      fi
    else
      let ERROR="$ERROR+1"
      echo ERROR
    fi
  done
  if [[ $ONLY_FIRST == 1 ]] ; then
    let RESULT="$ERROR + $FAIL"
    exit $RESULT
  fi
done

echo "PASS $PASS, FAIL $FAIL, ERROR $ERROR"
