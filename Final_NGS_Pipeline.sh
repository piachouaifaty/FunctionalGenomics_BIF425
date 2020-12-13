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
