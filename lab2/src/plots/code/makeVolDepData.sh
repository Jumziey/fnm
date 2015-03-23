#!/bin/bash
declare -a gases=("Ne" "Ar" "Kr" "Xe")
declare -a directions=("2 0 0" "2 2 0" "2 2 2")
res=200
Prop="VolDep"

for gas in "${gases[@]}" 
do
	for dir in "${directions[@]}" 
	do
		#Just removing whitespace
		dirnw="$(echo -e "${dir}" | tr -d '[[:space:]]')"
		
		commando="../../phonons $gas gamma 0 0 0 $dir $res"
		echo $commando
		eval $commando > ../data/"$gas$Prop$dirnw"
	done
done

