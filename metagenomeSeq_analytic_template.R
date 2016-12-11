## R script template for the production of basic qualitative
## statistics on the AMR and kraken output from AMR++

## Author: Steven Lakin


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


###################
## User Controls ##
###################
## Hopefully, this section should be the only code you need to modify.
## However, you can look into the code in further sections if you need
## to change other, more subtle variables in the exploratory or
## statistical functions.

# Set your working directory to the main folder for analysis:
setwd('..')

# Set the output directory for graphs:
graph_output_dir = 'graphs'


# Set the output directory for statistics:
stats_output_dir = 'stats'


# Where is the metadata file stored on your machine?
metadata_filepath = 'metadata.csv'


# In which column of the metadata file are the sample IDs stored?
sample_column_id = 'ID'


# The following is a list of column names in your metadata.csv file that you want
# to use for EXPLORATORY analysis (NMDS, PCA, alpha rarefaction, barplots)
exploratory_vars = c('Location', 'Type')


# The following is a list of strings specifying the formula for your
# model matrices (will be passed to model.matrix) that you want
# to use for STATISTICAL analyses.  An example is given below.
stats_model_matrices = list(
    '~ 0 + Type',
    '~ 0 + Type'
)


# The following is a list of character lists, where each string includes the variables that
# you desire contrasts for a given model matrix.  For instance, the second vector
# in this list of vectors corresponds to the second model matrix in the previous
# list of model matrices.  An example is given below.
stats_contrast_vars = list(
    list('TypeVar1 - TypeVar2',
      'TypeVar1 - TypeVar3',
      'TypeVar2 - TypeVar3'),

    list('TypeVar1 - TypeVar2',
         'TypeVar1 - TypeVar3',
         'TypeVar2 - TypeVar3')
)


# If you wish to use RANDOM EFFECTS in the statistical analysis,
# specify a single variable for each model matrix in the above model matrix
# list.  If no random effect is desired, then put NA.
stats_random_effect = c(NA, 'Location')


# Optional vector of names for your statistical model matrices.  These will be
# used to name the files for the statistical output.
stats_analysis_names = c('TypeNoRandom', 'TypeLocationRandom')


# The following are a list of column subsets to be evaluated for
# each statistical analysis.  If you desire a subset of the data
# to be analyzed, then put the sample metadata type and its
# value here for each analysis, else NA.  For example:
stats_analysis_subsets = list(
    list(NA, NA),
    list('Type', 'Var1')
)



####################
## Automated Code ##
####################
## Modify this as necessary, though you shouldn't need to for basic use.


## Check to see that the user input the stats lists correctly
if(length(unique(sapply(list(stats_analysis_subsets,
                             stats_model_matrices,
                             stats_contrast_vars,
                             stats_random_effect,
                             stats_analysis_names), FUN=length))) != 1) {
    stop('The length of all \"stats\" variables in the user input must be the same,
         i.e. the length of each input should equal the number of analyses performed.')
}

# Source the utility functions file, which should be in the scripts folder with this file
source('scripts/meg_utility_functions.R')

require(metagenomeSeq)
require(data.table)
require(ggplot2)
require(vegan)

set.seed(154)  # Seed the RNG, necessary for reproducibility


# We usually filter out genes with wild-type potential.  If you want to include these
# in your analysis, comment this vector out
snp_regex = c('ACRR',
              'CATB',
              'CLS',
              'DFRC',
              'DHFR',
              'DHFRIII',
              'DHFRIX',
              'EMBA',
              'embB',
              'EMBB',
              'EMBC',
              'EMBR',
              'ETHA',
              'FOLP',
              'GIDB',
              'GYRA',
              'gyrB',
              'GYRB',
              'INHA',
              'INIA',
              'INIC',
              'KASA',
              'LIAFSR',
              'LMRA',
              'MARR',
              'MEXR',
              'MEXZ',
              'mprF',
              'MPRF',
              'NDH',
              'omp36',
              'OMP36',
              'OMPF',
              'OPRD',
              'PARC',
              'parE',
              'PARE',
              'PGSA',
              'phoP',
              'PHOP',
              'PNCA',
              'POR',
              'PORB',
              'RAMR',
              'rpoB',
              'RPOB',
              'RPOC',
              'RPSL',
              'SOXS',
              'tetR',
              'TETR',
              'TLYA',
              'TUFAB')


##########################
## Import & Format Data ##
##########################
## These files should be standard for all analyses, as they are
## the output matrices from AMR++ nextflow.  Additionally,
## you will need to obtain the most recent megares annotations file
## from megares.meglab.org


# If subdirs for stats and exploratory variables don't exist, create them
ifelse(!dir.exists(file.path(graph_output_dir)), dir.create(file.path(graph_output_dir), mode='777'), FALSE)
ifelse(!dir.exists(file.path(stats_output_dir)), dir.create(file.path(stats_output_dir), mode='777'), FALSE)

for( dtype in c('AMR', 'Microbiome') ) {
    ifelse(!dir.exists(file.path(graph_output_dir, dtype)),
           dir.create(file.path(graph_output_dir, dtype), mode='777'), FALSE)
    
    for( v in 1:length(exploratory_vars) ) {
        ifelse(!dir.exists(file.path(graph_output_dir, dtype, exploratory_vars[v])),
               dir.create(file.path(graph_output_dir, dtype, exploratory_vars[v]), mode='777'), FALSE)
    }
    
    ifelse(!dir.exists(file.path(stats_output_dir, dtype)),
           dir.create(file.path(stats_output_dir, dtype), mode='777'), FALSE)
    
    for( a in 1:length(stats_analysis_names) ) {
        ifelse(!dir.exists(file.path(stats_output_dir, dtype, stats_analysis_names[a])),
               dir.create(file.path(stats_output_dir, dtype, stats_analysis_names[a]), mode='777'), FALSE)
    }
}



# Load the data, MEGARes annotations, and metadata
kraken <- newMRexperiment(read.table('kraken_analytic_matrix.csv', header=T, row.names=1, sep=','))
amr <- newMRexperiment(read.table('AMR_analytic_matrix.csv', header=T, row.names=1, sep=','))
annotations <- data.table(read.csv('megares_annotations.csv', header=T))
setkey(annotations, header)  # Data tables are SQL objects with optional primary keys

metadata <- read.csv(metadata_filepath, header=T)
metadata[, sample_column_id] <- make.names(metadata[, sample_column_id])


# Calculate normalization factors on the analytic data.
# We use Cumulative Sum Scaling as implemented in metagenomeSeq.
# You will likely get a warning about this step, but it's safe to ignore
cumNorm(kraken)
cumNorm(amr)


# Extract the normalized counts into data tables for aggregation
kraken_norm <- data.table(MRcounts(kraken, norm=T))
kraken_raw <- data.table(MRcounts(kraken, norm=F))
amr_norm <- data.table(MRcounts(amr, norm=T))
amr_raw <- data.table(MRcounts(amr, norm=F))


# Aggregate the normalized counts for AMR using the annotations data table, SQL
# outer join, and aggregation with vectorized lapply
amr_norm[, header :=( rownames(amr) ), ]
setkey(amr_norm, header)
amr_norm <- annotations[amr_norm]  # left outer join

amr_raw[, header :=( rownames(amr) ), ]
setkey(amr_raw, header)
amr_raw <- annotations[amr_raw]  # left outer join

# Remove groups that correspond to potentially wild-type genes
amr_raw <- amr_raw[!(group %in% snp_regex), ]
amr_norm<- amr_norm[!(group %in% snp_regex), ]


# Group the AMR data by level for analysis
amr_class <- amr_norm[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
amr_class_analytic <- newMRexperiment(counts=amr_class[, .SD, .SDcols=!'class'])
rownames(amr_class_analytic) <- amr_class$class

amr_class_raw <- amr_raw[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')]
amr_class_raw_analytic <- newMRexperiment(counts=amr_class_raw[, .SD, .SDcols=!'class'])
rownames(amr_class_raw_analytic) <- amr_class_raw$class

amr_mech <- amr_norm[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
amr_mech_analytic <- newMRexperiment(counts=amr_mech[, .SD, .SDcols=!'mechanism'])
rownames(amr_mech_analytic) <- amr_mech$mechanism

amr_mech_raw <- amr_raw[, lapply(.SD, sum), by='mechanism', .SDcols=!c('header', 'class', 'group')]
amr_mech_raw_analytic <- newMRexperiment(counts=amr_mech_raw[, .SD, .SDcols=!'mechanism'])
rownames(amr_mech_raw_analytic) <- amr_mech_raw$mechanism

amr_group <- amr_norm[, lapply(.SD, sum), by='group', .SDcols=!c('header', 'mechanism', 'class')]
amr_group_analytic <- newMRexperiment(counts=amr_group[, .SD, .SDcols=!'group'])
rownames(amr_group_analytic) <- amr_group$group

amr_group_raw <- amr_raw[, lapply(.SD, sum), by='group', .SDcols=!c('header', 'mechanism', 'class')]
amr_group_raw_analytic <- newMRexperiment(counts=amr_group_raw[, .SD, .SDcols=!'group'])
rownames(amr_group_raw_analytic) <- amr_group_raw$group

amr_gene_analytic <- newMRexperiment(
    counts=amr_norm[!(group %in% snp_regex),
                    .SD, .SDcols=!c('header', 'class', 'mechanism', 'group')])
amr_gene_raw_analytic <- newMRexperiment(
    counts=amr_raw[!(group %in% snp_regex),
                   .SD, .SDcols=!c('header', 'class', 'mechanism', 'group')])

rownames(amr_gene_analytic) <- amr_norm$header
rownames(amr_gene_raw_analytic) <- amr_raw$header


# Make long data frame for plotting with ggplot2
amr_melted_analytic <- rbind(melt_dt(MRcounts(amr_class_analytic), 'Class'),
                             melt_dt(MRcounts(amr_mech_analytic), 'Mechanism'),
                             melt_dt(MRcounts(amr_group_analytic), 'Group'),
                             melt_dt(MRcounts(amr_gene_analytic), 'Gene'))
amr_melted_raw_analytic <- rbind(melt_dt(MRcounts(amr_class_raw_analytic), 'Class'),
                                 melt_dt(MRcounts(amr_mech_raw_analytic), 'Mechanism'),
                                 melt_dt(MRcounts(amr_group_raw_analytic), 'Group'),
                                 melt_dt(MRcounts(amr_gene_raw_analytic), 'Gene'))



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

kraken_raw[, id :=(rownames(kraken)), ]
setkey(kraken_raw, id)
kraken_raw <- kraken_taxonomy[kraken_raw]  # left outer join


# Group the kraken data by level for analysis, removing NA entries
kraken_domain <- kraken_norm[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:8]
kraken_domain_analytic <- newMRexperiment(counts=kraken_domain[, .SD, .SDcols=!'Domain'])
rownames(kraken_domain_analytic) <- kraken_domain$Domain

kraken_domain_raw <- kraken_raw[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:8]
kraken_domain_raw_analytic <- newMRexperiment(counts=kraken_domain_raw[, .SD, .SDcols=!'Domain'])
rownames(kraken_domain_raw_analytic) <- kraken_domain_raw$Domain

kraken_phylum <- kraken_norm[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:8]
kraken_phylum_analytic <- newMRexperiment(counts=kraken_phylum[, .SD, .SDcols=!'Phylum'])
rownames(kraken_phylum_analytic) <- kraken_phylum$Phylum

kraken_phylum_raw <- kraken_raw[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:8]
kraken_phylum_raw_analytic <- newMRexperiment(counts=kraken_phylum_raw[, .SD, .SDcols=!'Phylum'])
rownames(kraken_phylum_raw_analytic) <- kraken_phylum_raw$Phylum

kraken_class <- kraken_norm[!is.na(Class) & Class != 'NA', lapply(.SD, sum), by='Class', .SDcols=!1:8]
kraken_class_analytic <- newMRexperiment(counts=kraken_class[, .SD, .SDcols=!'Class'])
rownames(kraken_class_analytic) <- kraken_class$Class

kraken_class_raw <- kraken_raw[!is.na(Class) & Class != 'NA', lapply(.SD, sum), by='Class', .SDcols=!1:8]
kraken_class_raw_analytic <- newMRexperiment(counts=kraken_class_raw[, .SD, .SDcols=!'Class'])
rownames(kraken_class_raw_analytic) <- kraken_class_raw$Class

kraken_order <- kraken_norm[!is.na(Order) & Order != 'NA', lapply(.SD, sum), by='Order', .SDcols=!1:8]
kraken_order_analytic <- newMRexperiment(counts=kraken_order[, .SD, .SDcols=!'Order'])
rownames(kraken_order_analytic) <- kraken_order$Order

kraken_order_raw <- kraken_raw[!is.na(Order) & Order != 'NA', lapply(.SD, sum), by='Order', .SDcols=!1:8]
kraken_order_raw_analytic <- newMRexperiment(counts=kraken_order_raw[, .SD, .SDcols=!'Order'])
rownames(kraken_order_raw_analytic) <- kraken_order_raw$Order

kraken_family <- kraken_norm[!is.na(Family) & Family != 'NA', lapply(.SD, sum), by='Family', .SDcols=!1:8]
kraken_family_analytic <- newMRexperiment(counts=kraken_family[, .SD, .SDcols=!'Family'])
rownames(kraken_family_analytic) <- kraken_family$Family

kraken_family_raw <- kraken_raw[!is.na(Family) & Family != 'NA', lapply(.SD, sum), by='Family', .SDcols=!1:8]
kraken_family_raw_analytic <- newMRexperiment(counts=kraken_family_raw[, .SD, .SDcols=!'Family'])
rownames(kraken_family_raw_analytic) <- kraken_family_raw$Family

kraken_genus <- kraken_norm[!is.na(Genus) & Genus != 'NA', lapply(.SD, sum), by='Genus', .SDcols=!1:8]
kraken_genus_analytic <- newMRexperiment(counts=kraken_genus[, .SD, .SDcols=!'Genus'])
rownames(kraken_genus_analytic) <- kraken_genus$Genus

kraken_genus_raw <- kraken_raw[!is.na(Genus) & Genus != 'NA', lapply(.SD, sum), by='Genus', .SDcols=!1:8]
kraken_genus_raw_analytic <- newMRexperiment(counts=kraken_genus_raw[, .SD, .SDcols=!'Genus'])
rownames(kraken_genus_raw_analytic) <- kraken_genus_raw$Genus

kraken_species <- kraken_norm[!is.na(Species) & Species != 'NA', lapply(.SD, sum), by='Species', .SDcols=!1:8]
kraken_species_analytic <- newMRexperiment(counts=kraken_species[, .SD, .SDcols=!'Species'])
rownames(kraken_species_analytic) <- kraken_species$Species

kraken_species_raw <- kraken_raw[!is.na(Species) & Species != 'NA', lapply(.SD, sum), by='Species', .SDcols=!1:8]
kraken_species_raw_analytic <- newMRexperiment(counts=kraken_species_raw[, .SD, .SDcols=!'Species'])
rownames(kraken_species_raw_analytic) <- kraken_species_raw$Species


# Make long data frame for plotting with ggplot2
kraken_melted_analytic <- rbind(melt_dt(MRcounts(kraken_domain_analytic), 'Domain'),
                                melt_dt(MRcounts(kraken_phylum_analytic), 'Phylum'),
                                melt_dt(MRcounts(kraken_class_analytic), 'Class'),
                                melt_dt(MRcounts(kraken_order_analytic), 'Order'),
                                melt_dt(MRcounts(kraken_family_analytic), 'Family'),
                                melt_dt(MRcounts(kraken_genus_analytic), 'Genus'),
                                melt_dt(MRcounts(kraken_species_analytic), 'Species'))
kraken_melted_raw_analytic <- rbind(melt_dt(MRcounts(kraken_domain_raw_analytic), 'Domain'),
                                    melt_dt(MRcounts(kraken_phylum_raw_analytic), 'Phylum'),
                                    melt_dt(MRcounts(kraken_class_raw_analytic), 'Class'),
                                    melt_dt(MRcounts(kraken_order_raw_analytic), 'Order'),
                                    melt_dt(MRcounts(kraken_family_raw_analytic), 'Family'),
                                    melt_dt(MRcounts(kraken_genus_raw_analytic), 'Genus'),
                                    melt_dt(MRcounts(kraken_species_raw_analytic), 'Species'))

# Ensure that the metadata entries match the factor order of the MRexperiments
metadata <- data.table(metadata[match(colnames(MRcounts(amr_class_analytic)), metadata[, sample_column_id]), ])
setkeyv(metadata, sample_column_id)


# Vector of objects for iteration and their names
AMR_analytic_data <- c(amr_class_analytic,
                       amr_mech_analytic,
                       amr_group_analytic,
                       amr_gene_analytic)
AMR_analytic_names <- c('Class', 'Mechanism', 'Group', 'Gene')
AMR_raw_analytic_data <- c(amr_class_raw_analytic,
                           amr_mech_raw_analytic,
                           amr_group_raw_analytic,
                           amr_gene_raw_analytic)
AMR_raw_analytic_names <- c('Class', 'Mechanism', 'Group', 'Gene')
kraken_analytic_data <- c(kraken_domain_analytic,
                          kraken_phylum_analytic,
                          kraken_class_analytic,
                          kraken_order_analytic,
                          kraken_family_analytic,
                          kraken_genus_analytic,
                          kraken_species_analytic)
kraken_analytic_names <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
kraken_raw_analytic_data <- c(kraken_domain_raw_analytic,
                              kraken_phylum_raw_analytic,
                              kraken_class_raw_analytic,
                              kraken_order_raw_analytic,
                              kraken_family_raw_analytic,
                              kraken_genus_raw_analytic,
                              kraken_species_raw_analytic)
kraken_raw_analytic_names <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')


for( l in 1:length(AMR_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(AMR_analytic_data[[l]])), metadata[[sample_column_id]])
    pData(AMR_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(AMR_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(AMR_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(AMR_analytic_data[[l]])))
    rownames(fData(AMR_analytic_data[[l]])) <- rownames(MRcounts(AMR_analytic_data[[l]]))
}

for( l in 1:length(AMR_raw_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(AMR_raw_analytic_data[[l]])), metadata[[sample_column_id]])
    pData(AMR_raw_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(AMR_raw_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(AMR_raw_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(AMR_raw_analytic_data[[l]])))
    rownames(fData(AMR_raw_analytic_data[[l]])) <- rownames(MRcounts(AMR_raw_analytic_data[[l]]))
}


for( l in 1:length(kraken_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(kraken_analytic_data[[l]])), metadata[[sample_column_id]])
    pData(kraken_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(kraken_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(kraken_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_analytic_data[[l]])))
    rownames(fData(kraken_analytic_data[[l]])) <- rownames(MRcounts(kraken_analytic_data[[l]]))
}

for( l in 1:length(kraken_raw_analytic_data) ) {
    sample_idx <- match(colnames(MRcounts(kraken_raw_analytic_data[[l]])), metadata[[sample_column_id]])
    pData(kraken_raw_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
    rownames(pData(kraken_raw_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
    fData(kraken_raw_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_raw_analytic_data[[l]])))
    rownames(fData(kraken_raw_analytic_data[[l]])) <- rownames(MRcounts(kraken_raw_analytic_data[[l]]))
}


#############################################
## Exploratory Analyses: Alpha Rarefaction ##
#############################################
for( v in 1:length(exploratory_vars) ) {
    # AMR
    meg_alpha_rarefaction(data_list=AMR_raw_analytic_data,
                          data_names=AMR_raw_analytic_names,
                          metadata=metadata,
                          sample_var=sample_column_id,
                          group_var=exploratory_vars[v],
                          outdir=paste(graph_output_dir, 'AMR', exploratory_vars[v],
                                       sep='/', collapse=''),
                          data_type='AMR')
    
    # Microbiome
    meg_alpha_rarefaction(data_list=kraken_raw_analytic_data,
                          data_names=kraken_raw_analytic_names,
                          metadata=metadata,
                          sample_var=sample_column_id,
                          group_var=exploratory_vars[v],
                          outdir=paste(graph_output_dir, 'Microbiome', exploratory_vars[v],
                                       sep='/', collapse=''),
                          data_type='Microbiome')
}


######################################
## Exploratory Analyses: Ordination ##
######################################
for( v in 1:length(exploratory_vars) ) {
    # AMR NMDS
    meg_ordination(data_list = AMR_analytic_data,
                   data_names = AMR_analytic_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = exploratory_vars[v],
                   outdir = paste(graph_output_dir, 'AMR', exploratory_vars[v],
                                  sep='/', collapse=''),
                   data_type = 'AMR',
                   method = 'NMDS')
    
    # AMR PCA
    meg_ordination(data_list = AMR_analytic_data,
                   data_names = AMR_analytic_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = exploratory_vars[v],
                   outdir = paste(graph_output_dir, 'AMR', exploratory_vars[v],
                                  sep='/', collapse=''),
                   data_type = 'AMR',
                   method = 'PCA')
    
    # Microbiome NMDS
    meg_ordination(data_list = kraken_analytic_data,
                   data_names = kraken_analytic_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = exploratory_vars[v],
                   outdir = paste(graph_output_dir, 'Microbiome', exploratory_vars[v],
                                  sep='/', collapse=''),
                   data_type = 'Microbiome',
                   method = 'NMDS')
    
    # Microbiome PCA
    meg_ordination(data_list = kraken_analytic_data,
                   data_names = kraken_analytic_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = exploratory_vars[v],
                   outdir = paste(graph_output_dir, 'Microbiome', exploratory_vars[v],
                                  sep='/', collapse=''),
                   data_type = 'Microbiome',
                   method = 'PCA')
}


####################################
## Exploratory Analyses: Heatmaps ##
####################################

# AMR Heatmaps for each level
for( v in 1:length(exploratory_vars) ) {
    for( l in 1:length(AMR_analytic_names) ) {
        meg_heatmap(melted_data=amr_melted_analytic,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_vars[v],
                    level_var=AMR_analytic_names[l],
                    outdir=paste(graph_output_dir, 'AMR', exploratory_vars[v],
                                 sep='/', collapse=''),
                    data_type='AMR')
    }
}

# Microbiome
for( v in 1:length(exploratory_vars) ) {
    for( l in 1:length(kraken_analytic_names) ) {
        meg_heatmap(melted_data=kraken_melted_analytic,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_vars[v],
                    level_var=kraken_analytic_names[l],
                    outdir=paste(graph_output_dir, 'Microbiome', exploratory_vars[v],
                                 sep='/', collapse=''),
                    data_type='Microbiome')
    }
}


####################################
## Exploratory Analyses: Barplots ##
####################################

# AMR
for( v in 1:length(exploratory_vars) ) {
    for( l in 1:length(AMR_analytic_names) ) {
        suppressWarnings(
            meg_barplot(melted_data=amr_melted_analytic,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_vars[v],
                    level_var=AMR_analytic_names[l],
                    outdir=paste(graph_output_dir, 'AMR', exploratory_vars[v],
                                 sep='/', collapse=''),
                    data_type='AMR')
        )
    }
}

# Microbiome
for( v in 1:length(exploratory_vars) ) {
    for( l in 1:length(kraken_analytic_names) ) {
        suppressWarnings(
            meg_barplot(melted_data=kraken_melted_analytic,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_vars[v],
                    level_var=kraken_analytic_names[l],
                    outdir=paste(graph_output_dir, 'Microbiome', exploratory_vars[v],
                                 sep='/', collapse=''),
                    data_type='Microbiome')
        )
    }
}


##########################
## Statistical Analyses ##
##########################
for( a in 1:length(stats_model_matrices) ) {
    meg_fitZig(data_list=AMR_analytic_data,
               data_names=AMR_analytic_names,
               metadata=metadata,
               zero_mod=model.matrix(~1 + log(libSize(amr))),
               data_mod=stats_model_matrices[a],
               filter_min_threshold=0.15,
               contrast_list=stats_contrast_vars[[a]],
               random_effect_var=stats_random_effect[a],
               outdir=paste(stats_output_dir, 'AMR', stats_analysis_names[a],
                            sep='/', collapse=''),
               analysis_name=stats_analysis_names[a],
               analysis_subset=stats_analysis_subsets[[a]],
               data_type='AMR',
               pval=0.05,
               top_hits=100)
    
    meg_fitZig(data_list=kraken_analytic_data,
               data_names=kraken_analytic_names,
               metadata=metadata,
               zero_mod=model.matrix(~1 + log(libSize(kraken))),
               data_mod=stats_model_matrices[a],
               filter_min_threshold=0.15,
               contrast_list=stats_contrast_vars[[a]],
               random_effect_var=stats_random_effect[a],
               outdir=paste(stats_output_dir, 'Microbiome', stats_analysis_names[a],
                            sep='/', collapse=''),
               analysis_name=stats_analysis_names[a],
               analysis_subset=stats_analysis_subsets[[a]],
               data_type='Microbiome',
               pval=0.05,
               top_hits=100)
}







