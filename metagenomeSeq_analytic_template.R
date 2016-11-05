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
require(data.table)

set.seed(154)  # Seed the RNG, necessary for reproducibility

# Load the data and MEGARes annotations
kraken <- newMRexperiment(read.table('kraken_analytic_matrix.csv', header=T, row.names=1, sep=','))
amr <- newMRexperiment(read.table('AMR_analytic_matrix.csv', header=T, row.names=1, sep=','))
annotations <- data.table(read.csv('megares_annotations.csv', header=T))
setkey(annotations, header)  # Data tables are SQL objects with optional primary keys


# Calculate normalization factors on the analytic data.
# We use Cumulative Sum Scaling as implemented in metagenomeSeq.
# You will likely get a warning about this step, but it's safe to ignore
cumNorm(kraken)
cumNorm(amr)


# Extract the normalized counts into a new data frame for aggregation
kraken_norm <- data.table(MRcounts(kraken, norm=T))
amr_norm <- data.table(MRcounts(amr, norm=T))


# Aggregate the normalized counts for AMR using the annotations data table, SQL
# outer join, and aggregation with vectorized lapply
amr_norm[, header :=( rownames(amr) ), ]
setkey(amr_norm, header)
amr_norm <- annotations[amr_norm]  # left outer join


# Group the AMR data by level for analysis
amr_class <- amr_norm[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
amr_class_analytic <- newMRexperiment(counts=amr_class[, .SD, .SDcols=!'class'])
rownames(amr_class_analytic) <- amr_class$class

amr_mech <- amr_norm[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
amr_mech_analytic <- newMRexperiment(counts=amr_mech[, .SD, .SDcols=!'mechanism'])
rownames(amr_mech_analytic) <- amr_mech$mechanism

amr_group <- amr_norm[, lapply(.SD, sum), by='group', .SDcols=!c('header', 'mechanism', 'class')]
amr_group_analytic <- newMRexperiment(counts=amr_group[, .SD, .SDcols=!'group'])
rownames(amr_group_analytic) <- amr_group$group

amr_gene_analytic <- amr
rm(amr)


# Aggregate the kraken data using the rownames:
# this set of commands splits the rownames into their taxonomic levels and
# fills empty values with NA.  We then join that taxonomy data table with
# the actual data and aggregate using lapply as before.
kraken_taxonomy <- data.table(id=rownames(kraken))
setDT(kraken_taxonomy)[, c('Domain',
                           'Phylum',
                           'Class',
                           'Order',
                           'Family',
                           'Genus',
                           'Species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]
setkey(kraken_taxonomy, id)
kraken_norm[, id :=(rownames(kraken)), ]
setkey(kraken_norm, id)
kraken_norm <- kraken_taxonomy[kraken_norm]  # left outer join


# Group the kraken data by level for analysis, removing NA entries
kraken_domain <- kraken_norm[!is.na(Domain), lapply(.SD, sum), by='Domain', .SDcols=!1:8]
kraken_domain_analytic <- newMRexperiment(counts=kraken_domain[, .SD, .SDcols=!'Domain'])
rownames(kraken_domain_analytic) <- kraken_domain$Domain

kraken_phylum <- kraken_norm[!is.na(Phylum), lapply(.SD, sum), by='Phylum', .SDcols=!1:8]
kraken_phylum_analytic <- newMRexperiment(counts=kraken_phylum[, .SD, .SDcols=!'Phylum'])
rownames(kraken_phylum_analytic) <- kraken_phylum$Phylum

kraken_class <- kraken_norm[!is.na(Class), lapply(.SD, sum), by='Class', .SDcols=!1:8]
kraken_class_analytic <- newMRexperiment(counts=kraken_class[, .SD, .SDcols=!'Class'])
rownames(kraken_class_analytic) <- kraken_class$Class

kraken_order <- kraken_norm[!is.na(Order), lapply(.SD, sum), by='Order', .SDcols=!1:8]
kraken_order_analytic <- newMRexperiment(counts=kraken_order[, .SD, .SDcols=!'Order'])
rownames(kraken_order_analytic) <- kraken_order$Order

kraken_family <- kraken_norm[!is.na(Family), lapply(.SD, sum), by='Family', .SDcols=!1:8]
kraken_family_analytic <- newMRexperiment(counts=kraken_family[, .SD, .SDcols=!'Family'])
rownames(kraken_family_analytic) <- kraken_family$Family

kraken_genus <- kraken_norm[!is.na(Genus), lapply(.SD, sum), by='Genus', .SDcols=!1:8]
kraken_genus_analytic <- newMRexperiment(counts=kraken_genus[, .SD, .SDcols=!'Genus'])
rownames(kraken_genus_analytic) <- kraken_genus$Genus

kraken_species <- kraken_norm[!is.na(Species), lapply(.SD, sum), by='Species', .SDcols=!1:8]
kraken_species_analytic <- newMRexperiment(counts=kraken_species[, .SD, .SDcols=!'Species'])
rownames(kraken_species_analytic) <- kraken_species$Species

rm(kraken)


# At this point, you will be working with the "analytic" MRexperiment objects.
# There should be 4 for the AMR data and 8 for kraken.  You likely won't be using
# them all, but you can if you wish.

# The following code produces PCA, PCoA, and NMDS plots by level for both AMR
# and kraken data.  By default, these graphs get output to your working directory.











