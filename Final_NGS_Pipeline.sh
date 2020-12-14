#!/bin/bash

#Copying the files to my directory
cp -R /mnt/gkhazen/NGS-Fall2020/FinalProject/* .
zcat 392_1.fastq.gz | more

#Number of lines in fastq file

zcat 392_1.fastq.gz | wc -l
#output: 123580348
#reads=30895087

zcat 392_2.fastq.gz | wc -l
#output: 123580348
#reads=30895087

#Reference chromosome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr13.fa.gz

#Running FastQC
time fastqc -o FastQCResults 392_1.fastq.gz  392_2.fastq.gz

#real    12m33.402s
#user    12m24.276s
#sys     0m14.858s

#Copy html report files into my machine
scp -r pia.chouaifaty@linuxdev.accbyblos.lau.edu.lb:FunctionalFinalProject/FastQCResults /Users/piachouaifaty

#Checking adapter files
cat TruSeq2-PE.fa #yes
cat TruSeq3-PE.fa #no
cat TruSeq3-PE-2.fa #yes

#trimming step to trim adapters
java -jar NGS/trim/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-threads 8 \
-trimlog ./FunctionalFinalProject/392.log \
./FunctionalFinalProject/392_1.fastq.gz \
./FunctionalFinalProject/392_2.fastq.gz \
./FunctionalFinalProject/392_1_trimmed_R1_paired.fastq.gz \
./FunctionalFinalProject/392_1_trimmed_R1_unpaired.fastq.gz \
./FunctionalFinalProject/392_2_trimmed_R2_paired.fastq.gz \
./FunctionalFinalProject/392_2_trimmed_R2_unpaired.fastq.gz \
ILLUMINACLIP:NGS/trim/Trimmomatic-0.39/adapters/TruSeq2-PE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:20 \
MINLEN:36

##specify full path to trimmomatic: mine is located
#in home/pia.chouaifaty/NGS/trim
##threads 8 on server
##trimlog: the log (errors saved there)
## r1raw r2raw r1paired r1unpaired r2paired r2unpaired
##Illuminaclip: location of the file with the adapters
#TruSeq2-PE.fa:2:30:10
#we allow 2 mismatches to the adapter sequence
#30 is for the palidrome clip threshold used only in PE
#require a score of at least 10 for alignment between adapter and read
##Trailing/leading 3: minimum quality required to keep a leading 5’ or trailing 3’ base
##SLIDINGWINDOW: Window size of 4 (BASES) minimum mean quality in window=20 since the reads are of good quality
##MINLEN: discard sequences that are smaller than 36 base pairs after the other trimming operations.
