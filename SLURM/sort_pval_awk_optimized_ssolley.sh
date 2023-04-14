#!/bin/bash

tmp=${4##*/}
pheno=`echo $tmp | cut -f2 -d'_'`

for ((i=$1; i<=$2; i++));
do
	genepair=`sed -n ${i}p $3 | cut -f1`
	score=`sed -n ${i}p $3 | cut -f3`
	if [[ "$score" == -* ]]
	then	
		oppscore=${score#-}
	else
		oppscore=-${score}
	fi
	n=`awk -v score="$score" '$1 > score {print NR;exit}' $4`
	nopp=`awk -v score="$oppscore" '$1 > score {print NR;exit}' $4`
	nless=$((n-1))
	ntotal=2898960
	ngreat=$((ntotal - nless))
	noppless=$((nopp - 1))
	noppgreat=$((ntotal - noppless))
	if [[ "$score" == -* ]]
        then
		noppless=$((nopp))
        	noppgreat=$((ntotal - noppless))
        else
		noppless=$((nopp - 1))
        	noppgreat=$((ntotal - noppless))
        fi
	echo $genepair $score $nless $ngreat $noppless $noppgreat >> /scratch/tmp/ssolley/pvals/${pheno}/${tmp%%_*}_genecombination_nulldist_counts_${1}_${2}.txt
done
