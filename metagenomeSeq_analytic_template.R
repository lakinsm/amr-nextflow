## Template script for starting a new statistical analysis
## using output from AmrPlusPlus: kraken & AMR

## The files you want to use for input to this (for the MEG group analyses)
## are the AMR_analytic_matrix.csv and kraken_analytic_matrix.csv.  The AMR
## matrix is identical to the Gene.csv matrix, however the kraken analytic
## matrix is not due to the way that reads get classified using each of
## these methods.

## So you should have pulled these files from the output of the nextflow pipeline
## and you are now performing this analysis on your local machine.  We will assume
## that you've set your working directory to where these files are located on your
## local machine and that you have installed the metagenomeSeq package.

## For the AMR analysis, you will also need to download the megares_annotations.csv
## file from the MEGARes website; the annotation file must be from the same version
## of the database as the file you used in the AmrPlusPlus pipeline, i.e. the headers
## must match between the annotation file and the database file.

require(metagenomeSeq)

# Load the data and MEGARes annotations
kraken <- newMRexperiment(read.table('kraken_analytic_matrix.csv', header=T, row.names=1, sep=','))
amr <- newMRexperiment(read.table('AMR_analytic_matrix.csv', header=T, row.names=1, sep=','))
annotations <- read.csv('megares_annotations.csv', header=T)


# Calculate normalization factors on the analytic data
cumNorm(kraken)
cumNorm(amr)


# Extract the normalized counts into a new data frame for aggregation
kraken_norm <- data.frame(MRcounts(kraken, norm=T))
amr_norm <- data.frame(MRcoefs(amr, norm=T))


# Aggregate the normalized counts for AMR using the annotations data frame




