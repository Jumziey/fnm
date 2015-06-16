#!/bin/bash

cd ..
for k in {150..200}
do
	command="potWellWavek0 parameters/potWellk0 -k=$k"
	echo $command
	eval $command
done
