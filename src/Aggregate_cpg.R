#!/usr/bin/env Rscript

#Script to take in directory of txt files of CpG WGBS methylation and annotated regions to output percent methylation at annotations
#Txt file format: chr, position, strand (arbitrary), frequency of C, coverage
#Usage: Rscript Aggregate_cpg.R <directory name with all cytosine files i.e chromosome> <name of anno file, eg. AR.loci.amplification> 

#Get command line arguments
args = commandArgs(trailingOnly=TRUE)

#Load libraries
library(tidyverse)
library(methylKit)
library(GenomicRanges)


############################## Set path and sampleid ####################################

#Set paths
wd <-  "/storage1/fs1/ha/Active/maherlab/saha.d/projects/Seq_pipelines/BS_seq/"
# bs_path <- paste0(wd, "data/Methylkit_files/")
bs_path <- paste0(wd, "data/PRAD_CX_reports/")



#List full path for files in directory
bsfiles <- as.list(list.files(path = paste0(bs_path, args[1]), full.names = T))

#List of sample ids
sampleid <- as.list(list.files(path = paste0(bs_path, args[1])))
#sampleid <- lapply(str_split(sampleid, ".merged."), FUN = function(i){i[1]})
sampleid <- lapply(str_split(sampleid, "_1_val_1_"), FUN = function(i){i[1]})

#Treatment vector
tx <- rep(1, length(sampleid))




############################# Read in data and do median normalization 3###############################

#Read in files in methylkit format - min coverage = 0
# obj <- methRead(bsfiles, sample.id=sampleid, treatment = tx,
#                 assembly="hg38", context="CpG",
#                 pipeline =list(fraction=FALSE,chr.col=1,start.col=2,end.col=2, coverage.col=5,strand.col=3,freqC.col =4), mincov =1)

obj <- methRead(bsfiles, sample.id=sampleid, treatment = tx,
assembly="hg38", context="CpG", pipeline='bismarkCytosineReport', mincov = 1)


#Normalize coverage
obj <- normalizeCoverage(obj)


############################ Read in annotations ########################################################

#Get annotations from GRanges object
#annotation <- read_tsv(file = paste0(wd, "data/hmr.bed"), col_names = c("chrom", "start", "end")) %>% GRanges()

##################################### Annotate CpG in annotations ###########################################


#Sum coverage that resides in annotations for all samples, no strandedness
#obj <- regionCounts(object = obj, regions = annotation)
 
#Inner join regions across samples 
obj <- unite(object = obj) 


#Calculate percent methylation 
res <- data.frame(Coord = paste0(obj$chr, ".", obj$start, ".", obj$end),  percMethylation(obj)) 



######################### Save results ################################################################################



#Remove duplicate regions and write to disk
#res[!duplicated(res), ]  %>% write.table(paste0(wd, "results/", args[1], ".prad.hmr.txt"), quote = F, sep = "\t", row.names = F)

# res %>% write.table(paste0(wd, "results/", args[1], ".mcrpc.cpg.txt" ),  quote = F, sep = "\t", row.names = F)
res %>% write.table(paste0(wd, "results/", args[1], ".prad.cpg.txt" ),  quote = F, sep = "\t", row.names = F)











