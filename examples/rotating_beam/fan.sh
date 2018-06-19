#!/bin/bash
STEP=0.1
EXE=./rotating_beam
echo "Executing for normalized omega:" $STEP
${EXE} Omega=0.1 | grep -v "number"| grep -v "Norm" | grep -v "\[0\]" > freqdata.dat
for i in {2..12}
do
    ANSWER=$(expr $STEP*$i | bc)
    echo "Executing for normalized omega:" $ANSWER
    ${EXE} Omega=${ANSWER} | grep -v "number"| grep -v "Norm" | grep -v "\[0\]" >> freqdata.dat
done
