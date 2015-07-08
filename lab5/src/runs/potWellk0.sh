#!/bin/bash

cd ..
for k in {150..300}
do
	runfile="parameters/potWellk0_$k"
	cat parameters/potWellk0 > $runfile
	echo "wf_output_text =wavefunc_potWellk0_$k.dat" >> $runfile
	
	kdiff=$(($k-150))
	nt=$((16000 - $kdiff*100))
	echo "nt =$nt" >> $runfile
	command="./potWellWavek0 $runfile -k=$k"
	echo $command
	#eval $command
done
