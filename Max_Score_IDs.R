#!/usr/bin/env Rscript

#Rscript Max_Score_IDs.R path

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("ShortRead", quietly = TRUE))
  BiocManager::install("ShortRead")

library(ShortRead)

#Max_Score_IDs

main = function(fstq_file_path)
  
{
  vectmeans=c() #empty vector that will contain average score for each read
  fstq_file=readFastq(fstq_file_path)
  n=length(sread(fstq_file)) #number of reads
  
  for (i in 1:n) #for each read
  {
    seq_vect=as.vector(sread(fstq_file)[i]) #the sequence as a vector
    slen=nchar(seq_vect) #length of the vector
    qual=as(quality(fstq_file)[i], "matrix")[,1:slen] #qualities of bases as a vector
    vectmeans=append(vectmeans, mean(qual)) #mean qs of read, add to mean quality vector
  }
  
  #which.max(vectmeans) #get index of read with highest average quality
  
  idxmax=which.max(vectmeans)
  
  max_score=vectmeans[idxmax] #maximum score
  
  print(paste0("Max Average quality: ", max_score)) 
  
  idxmaxes=which(vectmeans==max_score) #indeces of reads with the maximum score
  
  print("Indeces of reads with max average quality: ")
  print(idxmaxes)
  print(paste0("Read ID: ", as.vector(id(fstq_file)[c(idxmaxes)]))) #headers of reads with max score
}

args = commandArgs(trailingOnly = TRUE)
fstq_file_path = args[1]

main(fstq_file_path)