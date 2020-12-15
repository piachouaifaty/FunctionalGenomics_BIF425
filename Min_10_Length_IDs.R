#!/usr/bin/env Rscript

#Rscript Min_10_Length_IDs.R path

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("ShortRead", quietly = TRUE))
  BiocManager::install("ShortRead")

library(ShortRead)

#Min_10_Length_IDs
main = function(fstq_file_path)
{
  vectlen=c()
  fstq_file=readFastq(fstq_file_path)
  n=length(sread(fstq_file)) #number of reads
  
  for (i in 1:n) #for each read
  {
    seq_vect=as.vector(sread(fstq_file)[i]) #the sequence as a vector
    slen=nchar(seq_vect) #length of the vector
    vectlen=append(vectlen, slen) #add length to vector of lengths
  }
  
  vectlen=cbind(vectlen, 1:n) #matrix with length and index of reads
  colnames(vectlen)=c("length", "index")
  
  srtd=vectlen[order(vectlen[,"length"]),] #sort matrix by length
  min10=srtd[1:10,] #10 shortest
  min10idx=min10[,"index"] #indeces of 10 shortes
  min10ids=as.vector(id(fstq_file)[c(min10idx)]) #ids of 10 shortest
  
  min10=cbind(min10, min10ids)
  colnames(min10)=c("length", "index", "id")
  print("10 Shortest Reads")
  print(min10)
}

args = commandArgs(trailingOnly = TRUE)
fstq_file_path = args[1]

main(fstq_file_path)