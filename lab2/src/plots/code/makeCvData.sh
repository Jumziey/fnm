#!/bin/bash
declare -a gases=("Ne" "Ar" "Kr" "Xe")
declare -a temp=("1" "200")
res=400
prop="Cv"

for gas in "${gases[@]}" 
do
		#Just removing whitespace
		dirnw="$(echo -e "${temp}" | tr -d '[[:space:]]')"
		
		commando="../../phonons $gas cv ${temp[0]} ${temp[1]} $res"
		echo $commando
		eval $commando > ../data/"$gas$prop"
done

