#!/usr/bin/env Rscript

#Rscript Count_Hom_Hetero.R path

if (!requireNamespace("vcfR", quietly = TRUE))
  install.packages("vcfR")
library(vcfR)

main=function(vcfpath)
{
  vcf = read.vcfR(vcfpath, verbose = FALSE )
  opt=vcf@gt
  #hom_wild = c("0/0", "0|0") #homozygous wild type
  #hetero = c("0/1", "0|1") #heterozygous
  #hom_mut= c("1/1", "1|1") #homozygous mutant
  
  count_hom_wild=0
  count_hetero=0
  count_hom_mut=0
  
  for (i in 1:nrow(opt))
  {
    ind=strsplit(opt[i,"392"],":")[[1]][1] #split by first : 
    
    if (ind=="0/0"|ind=="0|0")
    {count_hom_wild=count_hom_wild+1}
    else if (ind=="0/1"|ind=="0|1")
    {count_hetero=count_hetero+1}
    else if (ind=="1/1"|ind=="1|1")
    {count_hom_mut=count_hom_mut+1}
    
  }
  
  print(paste0("Homozygous Wild Type: ", count_hom_wild))
  print(paste0("Heterozygous: ", count_hetero))
  print(paste0("Homozygous Mutant: ", count_hom_mut))
  
}

args = commandArgs(trailingOnly = TRUE)
vcfpath = args[1]

main(vcfpath)