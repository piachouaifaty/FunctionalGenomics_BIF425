#!/usr/bin/env Rscript

#Rscript Find_Read_Lengths.R path

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("ShortRead", quietly = TRUE))
  BiocManager::install("ShortRead")

library(ShortRead)

#Find_Read_Lengths
main = function(fstq_file_path)
{
  fstq_file=readFastq(fstq_file_path)
  l=length(sread(fstq_file)) #number of reads
  vectlen=c()
  
  for (i in 1:l)
  {
    seq_vect=as.vector(sread(fstq_file)[i]) #the sequence as a vector
    slen=nchar(seq_vect) #length of the vector
    vectlen=append(vectlen, slen)
  }
  
  print(paste0("The lengths after trimming (unique): "))
  return(sort(unique(vectlen)))
}

args = commandArgs(trailingOnly = TRUE)
fstq_file_path = args[1]

main(fstq_file_path)