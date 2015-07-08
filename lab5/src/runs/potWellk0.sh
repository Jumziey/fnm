#!/bin/bash

cd ..
for k in {150..300}
do
	runfile="parameters/potWellk0_$k"
	tail -n +2 parameters/potWellk0 > $runfile
	echo "wf_output_text =wavefunc_potWellk0_$k.dat" >> $runfile
	command="./potWellWavek0 $runfile -k=$k"
	echo $command
	eval $command
done
