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
