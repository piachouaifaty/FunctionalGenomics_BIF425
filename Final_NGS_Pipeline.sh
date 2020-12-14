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
ILLUMINACLIP:NGS/trim/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 \
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
#TruSeq3-PE-2.fa:2:30:10
#we allow 2 mismatches to the adapter sequence
#30 is for the palidrome clip threshold used only in PE
#require a score of at least 10 for alignment between adapter and read
##Trailing/leading 3: minimum quality required to keep a leading 5’ or trailing 3’ base
##SLIDINGWINDOW: Window size of 4 (BASES) minimum mean quality in window=20 since the reads are of good quality
##MINLEN: discard sequences that are smaller than 36 base pairs after the other trimming operations.

#Trimmomatic Output
Input Read Pairs: 30895087
Both Surviving: 30320928 (98.14%)
Forward Only Surviving: 0 (0.00%)
Reverse Only Surviving: 0 (0.00%)
Dropped: 574159 (1.86%)
TrimmomaticPE: Completed successfully

#FastQC Post-Trim
fastqc -o PostTrim_FastQCResults 392_1_trimmed_R1_paired.fastq.gz 392_2_trimmed_R2_paired.fastq.gz

scp -r pia.chouaifaty@linuxdev.accbyblos.lau.edu.lb:FunctionalFinalProject/PostTrim_FastQCResults /Users/piachouaifaty

#BWA

#Indexing
bwa index -p chr13bwaidx -a bwtsw chr13.fa
#-p filename, by convention genome|algo|idx
#-a index algo (bwtsw for long genomes and is for short ones)

#Output:
[bwt_gen] Finished constructing BWT in 72 iterations.
[bwa_index] 110.82 seconds elapse.
[bwa_index] Update BWT... 0.71 sec
[bwa_index] Pack forward-only FASTA... 0.81 sec
[bwa_index] Construct SA from BWT and Occ... 32.72 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index -p chr13bwaidx -a bwtsw chr13.fa
[main] Real time: 146.446 sec; CPU: 146.302 sec


#READ GROUP
@RG\tID:rg1\tSM:392\tPL:ILLUMINA\tLB:lib1\t:PU:HNLHYDSXX:1:GCCGGACA+TGTAAGAG

#ID = Read group identifier

#PU = Platform Unit
#The PU holds three types of information, the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}.
#The {FLOWCELL_BARCODE} refers to the unique identifier for a particular flow cell.
#The {LANE} indicates the lane of the flow cell and the {SAMPLE_BARCODE} is a sample/library-specific identifier.
#Although the PU is not required by GATK but takes precedence over ID for base recalibration if it is present.
#In the example shown earlier, two read group fields, ID and PU, appropriately differentiate flow cell lane, marked by .2,
#a factor that contributes to batch effects.

#SM = Sample
#The name of the sample sequenced in this read group. GATK tools treat all read groups with the same SM value as containing
#sequencing data for the same sample, and this is also the name that will be used for the sample column in the VCF file.
#Therefore it is critical that the SM field be specified correctly.
#When sequencing pools of samples, use a pool name instead of an individual sample name.

#PL = Platform/technology used to produce the read
#This constitutes the only way to know what sequencing technology was used to generate the sequencing data.
#Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO.

#LB = DNA preparation library identifier
#MarkDuplicates uses the LB field to determine which read groups might contain molecular duplicates,
#in case the same DNA library was sequenced on multiple lanes.

#MAPPING
bwa mem -t 16 \
-R '@RG\tID:rg1\tSM:392\tPL:ILLUMINA\tLB:lib1\t:PU:HNLHYDSXX:1:GCCGGACA+TGTAAGAG' \
/ref_chrom/chr13bwaidx \ #index file path #full path better
392_1_trimmed_R1_paired.fastq.gz  \ #file 1
392_2_trimmed_R2_paired.fastq.gz \ #file 2
> 392_aln.sam #redirect output to sam file

bwa mem -t 16 -R '@RG\tID:rg1\tSM:392\tPL:ILLUMINA\tLB:lib1\t:PU:HNLHYDSXX:1:GCCGGACA+TGTAAGAG' ref_chrom/chr13bwaidx 392_1_trimmed_R1_paired.fastq.gz 392_2_trimmed_R2_paired.fastq.gz > 392_aln.sam

Processed 429766 reads in 120.953 CPU sec, 7.628 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 16 -R @RG\tID:rg1\tSM:392\tPL:ILLUMINA\tLB:lib1\t:PU:HNLHYDSXX:1:GCCGGACA+TGTAAGAG ref_chrom/chr13bwaidx 392_1_trimmed_R1_paired.fastq.gz 392_2_trimmed_R2_paired.fastq.gz
[main] Real time: 1084.125 sec; CPU: 17303.121 sec

#Cleaning up and converting sam to bam
samtools fixmate -O bam 392_aln.sam 392_aln.bam

#GATK
#Validating
/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk --java-options "-Xmx16g" ValidateSamFile INPUT=392_aln.bam MODE=SUMMARY

No errors found
[Mon Dec 14 20:32:42 EET 2020] picard.sam.ValidateSamFile done. Elapsed time: 5.01 minutes.
Runtime.totalMemory()=2985820160
Tool returned:
0

#Sorting the SAM file
/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk --java-options "-Xmx16g" SortSam INPUT=392_aln.bam OUTPUT=392_sorted.bam SORT_ORDER=coordinate

[Mon Dec 14 21:41:00 EET 2020] picard.sam.SortSam done. Elapsed time: 8.33 minutes.
Runtime.totalMemory()=6858735616
Tool returned:
0

#Marking duplicates
/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk --java-options "-Xmx16g" MarkDuplicates INPUT=392_sorted.bam OUTPUT=392_dedup.bam METRICS_FILE=392.metrics

[Mon Dec 14 21:52:24 EET 2020] picard.sam.markduplicates.MarkDuplicates done. Elapsed time: 6.75 minutes.
Runtime.totalMemory()=14134280192
Tool returned:
0
