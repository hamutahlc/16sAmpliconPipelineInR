##clear objects from environment
rm(list=ls())

##set working directory
setwd('~/Dropbox/sunflower/analysis/microbiome')

##We are going to use an open-source software packages from the
##Bioconductor Project (Huber et al. 2015), which includes CRAN
##packages but also DADA2, a processing workflow for recovering
##Amplicon Sequence Variants (ASVs) from amplicon data. It has been
##suggested that ASVs replace OTUs in marker-gene data analysis and
##reporting (Callahan et al. 2017 Nature). While OTUs are produced
##when you cluster reads that differ less than some fixed threshhold,
##ASVs will resolve amplicon reads to the level of single-nucleotide
##differences over the sequenced region. We will be modifying the
##workflow published by Callahan et al. 2016
##doi.org/10.12688/f1000research.8986.2 I modified this pipeline from
##https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
##For notes on where to find packages for different analyses available
##both in and outside of this pipeline, check out
##https://github.com/mdozmorov/Microbiome_notes

##load in Bioconductor (works with R 3.6.0 or higher)

## move to a seperate script
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

##install packages

install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.14") # change the ref argument to get other versions

## not sure if these work yet
library("knitr")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}

## Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

##set seed: We may need a random number sequence later, but we want to
##use the same sequence throughout this pipeline. So we will set a
##seed at 100.  later, we want to remove this and see if the code
##still works
set.seed(100)


##unzip fastq.gz files
##read in data. CHANGE to the directory containing the fastq files


miseq_path <- "~/Documents/sunflowerdata"

##Quality filter and trim your files.  These steps remove low quality
##sequences, and remove adapters, primers, PhiX libraries, etc. You
##also will convert your data into different formats

##First, examine the quality of your reads: you'll need to read in the
## names of your fastq files, and perform string manipulation to get
## lists of forward and reverse fastq files in matched order.  The
## sort function ensures forward/reverse reads are in the same order

#the $ sign indicates that we want the string to end with fastq so that we dont pull up other files in our directory like zips and things

fnFs <- sort(list.files(miseq_path, pattern="_pair1.fastq$"))
fnRs <- sort(list.files(miseq_path, pattern="_pair2.fastq$"))

##Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "flowcell1179lane1_"), `[`, 1)

##I don't know if this is applicable, our sample names are found in
##sample file maps, and we might need to give RBCL and 16s Samples
##that have the same ID new names so that we can compare

## Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)


##Most Illumina sequencing data shows reduced quality towards the end
##of the reads. To determine where to truncate reads, plot quality
##profiles First two forward reads
fwplotF<-dada2:::plotQualityProfile(fnFs[1:2])
fwplotF

##We want to truncate reads that fall below a quality score of 30, so we are going to add lines to figure out what cycle that would be

fwplotF + ggplot2:::geom_vline(xintercept=267) +
  ggplot2:::geom_hline(yintercept=30)


##First two reverse reads
fwplotR<-dada2:::plotQualityProfile(fnRs[1:2])
fwplotR
fwplotR + ggplot2:::geom_vline(xintercept=267) +
  ggplot2:::geom_hline(yintercept=30)


##Based on these plots, determine where to truncate. We found that
##quality drops at 267 for forward reads and also 267 for reversee reads.
##Often, the first 10 nucleotides of each
##read are likely to contain errors as well. There are other filtering
##parameters: commonly, we should allow a maximum of 2 expected errors
##per read (Edgar and Flyvbjerg 2015). We will perform these filter
##steps on paired ends jointly, meaning both reads must pass the
##filter for the whole pair to pass. We also want to remove PhiX libraries

##Define a filename for the filtered fastq.gz. Place filtered files in
##filtered/ subdirectory

filt_path <- file.path(miseq_path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0("sampleNames", "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0("sampleNames", "_R_filt.fastq.gz"))

##filter the reverse and forward reads truncLen: truncates reads after
##specified number of bases. Will discard reads shorter than these
##specified lengths.
##maxEE: reads with more than 2 errors will be discarded
##maxN: ???  
##truncQ: truncates reads at the first instance of a quality score less than what it is set to.  Default is 2.
##rm.phix: remove reads that match the PhiX genome (a reference
##library used in illumina prep) 
##compress: if TRUE, output is gzipped
##multithread: if TRUE, input files are filtered in parallel by
##forking. Windows cannot perform this function. If using Windows, set
##to FALSE

out <- dada2:::filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(267,267),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)



##Currently, your sequence run consists of hundreds of samples pooled
##together and in "fastq" file format.  Demultiplex Step: You must
##sort your sequence reads into separate files for each barcoded
##sample, so that you have 1 fastq file for each sample Dereplication
##Step: combines all identical sequencing reads into into “unique
##sequences” with a corresponding “abundance”: the number of reads
##with that unique sequence. Dereplication substantially reduces
##computation time by eliminating redundant comparisons. We name the
##resulting "derep-class" objects by their sample name.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
## Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

##DADA2 can distingish sequencing errors from real, biological
## variation by using a parameterized model of substitution errors the
## parameters for this model come from your data itself by using
## unsupervising learning in which sample inference is alternated with
## parameter estimation until both are jointly consistent

##Parameter learning is computationally intensive, as it requires
##multiple iterations of the sequence inference algorithm, and
##therefore it is often useful to estimate the error rates from a
##(sufficiently large) subset of the data

##Use the function learnerrors. If multithread is set to TRUE,
##multithreading is enabled and the number of available threads is
##automatically determined

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)


##In order to verify that the error rates have been reasonably
##well-estimated, we inspect the fit between the observed error rates
##(black points) and the fitted error rates (black lines)

##These figures will show the frequencies of each type of transition
##as a function of the quality.

plotErrors(errF)

plotErrors(errR)

## We now need to infer ASVs, amplicon sequence variants. (This is an
## alternative that DADA2 offers to the traditional method of
## clustering OTUs) For more information describing the
## high-resolution methods DADA2 uses to infer ASVs, refer to Benjamin
## J Callahan et al. 2016 You can run this step and simultaneously
## remove substitution and indel errors The dada function takes as
## input dereplicated amplicon sequencing reads and returns the
## inferred composition of samples.Put another way, dada removes all
## sequencing errors to reveal the members of the sequenced community
## This method can run in two different modes: Independent inference
## by sample (pool=FALSE), and inference from the pooled sequencing
## reads from all samples (pool=TRUE).

## Independent inference has the advantage that computation time is
## linear in the number of samples, and memory requirements are flat
## with the number of samples. This allows scaling out to datasets of
## almost unlimited size.  Pooled inference is more computationally
## taxing, and can become intractable for datasets of tens of millions
## of reads. However, pooling improves the detection of rare variants
## that were seen just once or twice in an individual sample but many
## times across all samples.  As this dataset is not particularly
## large, we perform pooled inference. As of version 1.2,
## multithreading can now be activated with the arguments multithread
## = TRUE, which substantially speeds this step.  Do I need to add
## "pool=TRUE" to these arguements?

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

##Inspect in the dada-class returned
dadaFs[[1]]
dadaRs[[1]]

##Read the output

#### E.g.
#### 128 sample sequences were inferred from 1979 input unique sequences.
#### Key parameters: OMEGA_A = 1e-40, BAND_SIZE = 16, USE_QUALS = TRUE

##This means The DADA2 algorithm inferred 128 real sequence variants from the 1979 unique sequences in the first sample.
##The dada-class object contains multiple diagnostics about the quality of each inferred sequence variant(see help("dada-class") for some info).

##We now merge together the inferred forward and reverse sequences, removing paired sequences that do not perfectly overlap as a final control against residual errors.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

##Make a sequence table
##Unsure how to modify this for our data
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))

##We need to do some more filtering because the sequence inference algorithm used by DADA2 does not remove chimeras () errors from incomplete PCR product)
##The sequence table we have includes chimeras, which we can remove by comparing each inferred sequence to others in the table, and removing those that can be reproduced by stitching together two more abundant sequences

eqtabNoC <- removeBimeraDenovo(seqtabAll)

##How about removing chloroplast reads, rare reads, and DNA from controls?


##Compare data to a standard database of bacteria and label them. You can use these labeled sequences to build a de novo phylogenetic tree to relate the sample sequences.
##Let's use the XXX Database from the XXX Institute

##





##Downstream analyses. There are many helpful packages in R for analyzing and visualizing microbiome data.
##We will use the "Microbiome", "MicrobiomeSeq", and "Phyloseq", which is especially useful for comparative analyses, dimensionality reduction methods and visualization, and plotting with ggplot2 graphics
