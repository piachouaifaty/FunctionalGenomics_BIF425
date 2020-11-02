#!/bin/bash

#$1: fastq file
#$2: specified length

let n=$(zcat $1 | wc -l) #number of lines
let i=2 #index of first sequence
let count=0 #count of sequences with specified length
declare -a avgs=() #array with all averages

while $i -lt n;
do
  lngth=$(sed -n "${i}p" $1 | wc) #get sequence length
  if [[$lngth=$2]]; #if sequence length is equal to specified length
  then count=$count+1 #increment count
        let sc=$i+2 #the index for the sequence's QS
        declare a- seq_v=(sed "${sc}q;d" $1 | od -An -t d1) #QS values for each nucleotide
        sum=$(IFS=+; echo "$((${seq_v[*]}))") #sum of QS values for this sequence
        avg=$((${sum}/${lngth}))
        avgs+=($avg) #add this average to the list of averages
    fi
let i=$i+4; #indeces of sequences
done

echo "there are $count sequences with length = $2\n"
echo "their average quality scores are:\n"
printf '%s\n' "${avgs[@]}" #prints each average score for each sequence
