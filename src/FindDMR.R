#!/usr/bin/env Rscript

#Script to take in directory of txt files of CpG WGBS methylation and output DMR, ie. reads in data using Methylkit and calculates DMR using DSS. 
#Txt file format: chr, position, strand (arbitrary), frequency of C, coverage
#Usage: Rscript FindDMR.R <directory name with all cytosine files i.e chromosome> <name of anno file, eg. AR.loci.amplification> 

#docker(dsaha0295/methylkit:latest)

#Get command line arguments
args = commandArgs(trailingOnly=TRUE)

#Load libraries
library(tidyverse)
library(methylKit)
library(DSS)



################################# Read in file names #######################################################


#Path for Methylkit format CX files from Bismark
bs_path <- "/storage1/fs1/ha/Active/maherlab/saha.d/projects/Seq_pipelines/BS_seq/data/Methylkit_files/"

#List full path for files in directory
bsfiles <- as.list(list.files(path = paste0(bs_path, args[1]), full.names = T))

#Get sample ids
sampleid <- as.list(list.files(path = paste0(bs_path, args[1])))
sampleid <- lapply(str_split(sampleid, ".merged."), FUN = function(i){i[1]})





################################ Read in annotations and subset  ####################################




#Read in annoical annotations df for WGS/WGBS samples and subset sampleids/files/anno list to those in common
#anno format: rownames = sampleid, columns = annotations for each sample
anno <- readRDS(paste0("/storage1/fs1/ha/Active/maherlab/saha.d/projects/lncRNA_analysis/data/Histo_mCRPC_WGBS_anno.rds"))

ids <- intersect(row.names(anno), sampleid)


bsfiles <- bsfiles[sampleid %in% ids] #Uses partial string matching
sampleid <- sampleid[sampleid %in% ids]
anno <- anno[row.names(anno) %in% ids,, drop = F]

group1 <- row.names(anno)[anno$Histology == "small_cell"]
group2 <- row.names(anno)[anno$Histology == "adeno"]


################################ Read in CpG files and format ########################################################


#Treatment vector (it's arbitrary)
tx <- rep(1, length(sampleid))


#Read in files in methylkit format -keep min coverage at 0 to get more CpGs
obj <- methRead(bsfiles, sample.id=sampleid, treatment = tx,
                assembly="hg38", context="CpG",
                pipeline =list(fraction=FALSE,chr.col=1,start.col=2,end.col=2, coverage.col=5,strand.col=3,freqC.col=4), 
                mincov = 1)

#Format each df per sample for DSS - <chr> <position> <Coverage> <# Methylation> - return list of df
meth <- lapply(1:length(sampleid), FUN = function(i){
  
  df <- obj[[i]] %>% GRanges() %>% as.data.frame() %>% dplyr::select(c(seqnames, start, coverage, numCs))
  colnames(df) <- c("chr", "pos", "N", "X")
  return(df)  
  
})
names(meth) <- sampleid




########################### Call DSS to detect DML/DMR ##############################################################



#Create BS object
bsobj <- makeBSseqData(dat = meth,sampleNames = names(meth))

#Use annoical anno df as design for LM formula - adjust design ***
#DMLfit = DMLfit.multiFactor(bsobj, design=anno, formula= ~ AR , smoothing = T)#Need to adjust formula per comparison ***
dmlTest.sm = DMLtest(bsobj, group1=group1, group2=group2, smoothing=TRUE)

#Test if coef of interest is different from 0 (first term is intercept) for each DML  - return significant DML
#DMLtest.cell = DMLtest.multiFactor(DMLfit, coef=2)
dmr = callDMR(dmlTest.sm, p.threshold=0.01)

#Aggregate significant DML into DMR
#dmr <- callDMR(DMLresult = DMLtest.cell, p.threshold = 1e-2) #Adjust pvalue threshold ***

#Write DMR bed to disk - need to use Aggregate_cpg.R to find methylation values 
write.table(x = dmr, file = paste0("/storage1/fs1/ha/Active/maherlab/saha.d/projects/lncRNA_analysis/results/", args[1], ".dmr.Histo.bed"), sep = "\t", quote = F, row.names = F)

