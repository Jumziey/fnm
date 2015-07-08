#!/bin/bash

cd ..
for k in {150..300}
do
	runfile="potWellk0_$k"
	tail -n +2 parameters/potWellk0 > parameters/$runfile
	echo "wf_output_text =wavefunc_potWellk0_$k.dat" >> parameters/$runfile
	command="./potWellWavek0 parameters/$runfile -k=$k"
	echo $command
	eval $command
done
