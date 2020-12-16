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

#TILES!! all good

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

###SAMTOOLS STATS

#-f only output reads with that bit
#-F only output reads WITHOUT that bit

#checking duplicates (possible only after MarkDuplicates)
samtools view -f 0x400 392_dedup.bam

#to count them without printing them
samtools view -c -f 0x400 392_dedup.bam
934746

# get the total number of reads of a BAM file (may include unmapped and duplicated multi-aligned reads)
#(unmapped reads + each aligned location per mapped read
samtools view -c 392_dedup.bam
61204902

# counting only mapped (primary aligned) reads
samtools view -c -F 260 392_dedup.bam
#read unmapped (0x4)
#not primary alignment (0x100)
#Flag: 260
#-F to exclude them
5805578

#unique (without multimapping)
#read unmapped (0x4)
#excludes unmapped reads, sorts and keeps only unique
samtools view -F 0x4 392_dedup.bam | cut -f 1 | sort | uniq | wc -l
2619706

(mapped/total)*100

#Number of reads without a pair complement
#An alignment with an unmapped mate is marked with a ‘*’ in column 7.
samtools view 392_dedup.bam | cut -f 7 | grep -c '*'
55399324

#Use the CIGAR string, to compute the number of reads without any Insertion or Deletion
#WITH INSERTIONS/DELETIONS
#column 6 has insertions and deletions
samtools view 392_dedup.bam | cut -f 6 | grep -c -E 'I|D'
997161
#WITHOUT = TOTAL - WITH

#- Number of supplementary and secondary reads
#supplementary alignment (0x800)
samtools view -c -f 0x800 392_dedup.bam
566166

#Average Mapping score/quality for the mapped reads
#-F 0X4 excludes unmapped reads
#awk gets the MapQ
samtools view -F 0x4 392_dedup.bam |  awk '{sum+=$5} END {print "Mean MAPQ =",sum/NR}'
Mean MAPQ = 19.1783

###END

#Realignment, need an index file and a dictionary file
/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk CreateSequenceDictionary \
R=ref_chrom/chr13.fa \
O=ref_chrom/chr13.dict

/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk CreateSequenceDictionary R=ref_chrom/chr13.fa O=ref_chrom/chr13.dict

[Mon Dec 14 22:14:21 EET 2020] picard.sam.CreateSequenceDictionary done. Elapsed time: 0.03 minutes.
Runtime.totalMemory()=2084569088
Tool returned:
0

samtools faidx ref_chrom/chr13.fa

#BaseRecalibrator
/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk BaseRecalibrator \
-I 392_dedup.bam \
-R ref_chrom/chr13.fa \
--known-sites /mnt/NGSdata/snpdb151_All_20180418.vcf \
-O recal_data.table

/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk BaseRecalibrator -I 392_dedup.bam -R ref_chrom/chr13.fa --known-sites /mnt/NGSdata/snpdb151_All_20180418.vcf -O recal_data.table

[December 15, 2020 4:19:37 PM EET] org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator done. Elapsed time: 6.55 minutes.
Runtime.totalMemory()=13362003968
Tool returned:
SUCCESS

#ApplyBQSR

/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk ApplyBQSR \
-R ref_chrom/chr13.fa \
-I 392_dedup.bam \
--bqsr-recal-file recal_data.table \
-O output_recal.bam

/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk ApplyBQSR -R ref_chrom/chr13.fa -I 392_dedup.bam --bqsr-recal-file recal_data.table -O output_recal.bam

[December 15, 2020 4:52:01 PM EET] org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator done. Elapsed time: 6.26 minutes.
Runtime.totalMemory()=14309392384
Tool returned:
SUCCESS

#HaplotypeCaller
#HaplotypeCaller runs per-sample to generate an intermediate GVCF (not to be used in final analysis), which can then be used in GenotypeGVCFs for joint genotyping of multiple samples in a very efficient way.

/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk --java-options "-Xmx16g" HaplotypeCaller \
-R ref_chrom/chr13.fa \
-I output_recal.bam \
-O output.g.vcf.gz \
-ERC GVCF

/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk --java-options "-Xmx16g" HaplotypeCaller -R ref_chrom/chr13.fa -I output_recal.bam -O output.g.vcf.gz -ERC GVCF

#gvcf file instead of vcf, for grouping samples according to genotype / case/controls
[December 15, 2020 7:48:11 PM EET] org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 72.30 minutes.
Runtime.totalMemory()=3538944000

#GenotypeGVCF
#joint genotyping on a single input, which may contain one or many samples

/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk --java-options "-Xmx16g" GenotypeGVCFs \
-R ref_chrom/chr13.fa \
-V output.g.vcf.gz \ #input of this command is output of Haplotype Caller
-O output.vcf.gz

/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk --java-options "-Xmx16g" GenotypeGVCFs -R ref_chrom/chr13.fa -V output.g.vcf.gz -O output.vcf.gz

[December 15, 2020 8:06:46 PM EET] org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs done. Elapsed time: 0.92 minutes.
Runtime.totalMemory()=5145886720

#The last two steps didn't work
#Read GATK Documentation

#counts all variants
/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk CountVariants -V output.vcf.gz
#OUTPUT
[December 16, 2020 1:09:30 AM EET] org.broadinstitute.hellbender.tools.walkers.CountVariants done. Elapsed time: 0.02 minutes.
Runtime.totalMemory()=2227699712
Tool returned:
22711

#selects SNPS
/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk SelectVariants -R ref_chrom/chr13.fa -V output.vcf.gz  --select-type-to-include SNP -O SNP_392.vcf
[December 16, 2020 1:12:58 AM EET] org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants done. Elapsed time: 0.04 minutes.
Runtime.totalMemory()=2216689664
#counts SNPS
/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk CountVariants -V SNP_392.vcf
[December 16, 2020 1:16:23 AM EET] org.broadinstitute.hellbender.tools.walkers.CountVariants done. Elapsed time: 0.01 minutes.
Runtime.totalMemory()=2205155328
Tool returned:
21953

#selects INDELs
/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk SelectVariants -R ref_chrom/chr13.fa -V output.vcf.gz  --select-type-to-include INDEL -O INDEL_392.vcf
[December 16, 2020 1:17:10 AM EET] org.broadinstitute.hellbender.tools.walkers.variantutils.SelectVariants done. Elapsed time: 0.02 minutes.
Runtime.totalMemory()=2123890688
#counts INDELs
/mnt/gkhazen/NGS-Fall2020/gatk-4.1.9.0/gatk CountVariants -V INDEL_392.vcf
[December 16, 2020 1:17:51 AM EET] org.broadinstitute.hellbender.tools.walkers.CountVariants done. Elapsed time: 0.01 minutes.
Runtime.totalMemory()=2286944256
Tool returned:
755

#Scripts

scp -r pia.chouaifaty@linuxdev.accbyblos.lau.edu.lb:FunctionalFinalProject/output.vcf.gz /Users/piachouaifaty

#Subset Test Files
zcat 392_1_trimmed_R1_paired.fastq.gz | head -n 1000 > test1.fastq
zcat 392_2_trimmed_R2_paired.fastq.gz | head -n 1000 > test2.fastq
scp -r pia.chouaifaty@linuxdev.accbyblos.lau.edu.lb:FunctionalFinalProject/test* /Users/piachouaifaty

#do on home dir
Rscript Find_Read_Lengths.R test1.fastq
Rscript Min_10_Length_IDs.R test1.fastq
Rscript Max_Score_IDs.R test1.fastq

time Rscript Count_Hom_Hetero.R output.vcf

#[1] "Homozygous Wild Type: 0"
#[1] "Heterozygous: 4253"
#[1] "Homozygous Mutant: 18408"
#real	0m2.585s
#user	0m2.296s
#sys	0m0.224s
