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
metadata_filepath = 'project_metadata.csv'


# Name of the megares annotation file used for this project
megares_annotation_filename = 'megares_annotations_v1.01.csv'


# In which column of the metadata file are the sample IDs stored?
sample_column_id = 'ID'


# The following is a list of analyses based on variables in 
# your metadata.csv file that you want
# to use for EXPLORATORY analysis (NMDS, PCA, alpha rarefaction, barplots)
exploratory_analyses = list(
    # Analysis 1
    # Description: 
    list(
        name = 'Variable1',
        subsets = list(),
        exploratory_var = 'Variable1'
    ),
    
    # Analysis 2
    # Description: 
    list(
        name = 'Variable2_Variable1_Subset',
        subsets = list('Variable1 == Value1'),
        exploratory_var = 'Variable2'
    )
)


# Each analyses you wish to perform should have its own list in the following
# statistical_analyses list.  A template is provided to get you started.
# Multiple analyses, subsets, and contrasts are valid, but only one random
# effect can be used per analysis.  The contrasts of interest must have their
# parent variable in the model matrix equation.  Contrasts are named by
# parent variable then child variable without a space inbetween, for example:
# PVar1Cvar1 where the model matrix equation is ~ 0 + Pvar1.
statistical_analyses = list(
  # Analysis 1
  # Description: 
  list(
    name = 'Variable1',
    subsets = list(),
    model_matrix = '~ 0 + Variable1 + Variable2',
    contrasts = list('Variable1Value1 - Variable1Value2'),
    random_effect = NA
  ),
  
  # Analysis 2
  # Description: 
  list(
    name = 'Variable2_Variable1_Subset',
    subsets = list('Variable1 == Value1'),
    model_matrix = '~ 0 + Variable2',
    contrasts = list('Variable2Value1 - Variable2Value2',
                     'Variable2Value1 - Variable2Value3',
                     'Variable2Value2 - Variable2Value3'),
    random_effect = NA
  )
)






####################
## Automated Code ##
####################
## Modify this as necessary, though you shouldn't need to for basic use.

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
    
    for( v in 1:length(exploratory_analyses) ) {
        ifelse(!dir.exists(file.path(graph_output_dir, dtype, exploratory_analyses[[v]]$name)),
               dir.create(file.path(graph_output_dir, dtype, exploratory_analyses[[v]]$name), mode='777'), FALSE)
    }
    
    ifelse(!dir.exists(file.path(stats_output_dir, dtype)),
           dir.create(file.path(stats_output_dir, dtype), mode='777'), FALSE)
    
    for( a in 1:length(statistical_analyses) ) {
        ifelse(!dir.exists(file.path(stats_output_dir, dtype, statistical_analyses[[a]]$name)),
               dir.create(file.path(stats_output_dir, dtype, statistical_analyses[[a]]$name), mode='777'), FALSE)
    }
}

ifelse(!dir.exists(file.path('amr_matrices')), dir.create(file.path('amr_matrices'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('kraken_matrices')), dir.create(file.path('kraken_matrices'), mode='777'), FALSE)

ifelse(!dir.exists(file.path('amr_matrices/sparse_normalized')), dir.create(file.path('amr_matrices/sparse_normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('amr_matrices/normalized')), dir.create(file.path('amr_matrices/normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('amr_matrices/raw')), dir.create(file.path('amr_matrices/raw'), mode='777'), FALSE)

ifelse(!dir.exists(file.path('kraken_matrices/sparse_normalized')), dir.create(file.path('kraken_matrices/sparse_normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('kraken_matrices/normalized')), dir.create(file.path('kraken_matrices/normalized'), mode='777'), FALSE)
ifelse(!dir.exists(file.path('kraken_matrices/raw')), dir.create(file.path('kraken_matrices/raw'), mode='777'), FALSE)



# Load the data, MEGARes annotations, and metadata
temp_kraken <- read.table('kraken_analytic_matrix.csv', header=T, row.names=1, sep=',')
kraken <- newMRexperiment(temp_kraken[rowSums(temp_kraken) > 0, ])
amr <- newMRexperiment(read.table('AMR_analytic_matrix.csv', header=T, row.names=1, sep=','))
annotations <- data.table(read.csv(megares_annotation_filename, header=T))
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
for( v in 1:length(exploratory_analyses) ) {
    # AMR
    meg_alpha_rarefaction(data_list=AMR_raw_analytic_data,
                          data_names=AMR_raw_analytic_names,
                          metadata=metadata,
                          sample_var=sample_column_id,
                          group_var=exploratory_analyses[[v]]$exploratory_var,
                          analysis_subset=exploratory_analyses[[v]]$subsets,
                          outdir=paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                       sep='/', collapse=''),
                          data_type='AMR')
    
    # Microbiome
    meg_alpha_rarefaction(data_list=kraken_raw_analytic_data,
                          data_names=kraken_raw_analytic_names,
                          metadata=metadata,
                          sample_var=sample_column_id,
                          group_var=exploratory_analyses[[v]]$exploratory_var,
                          analysis_subset=exploratory_analyses[[v]]$subsets,
                          outdir=paste(graph_output_dir, 'Microbiome', exploratory_analyses[[v]]$name,
                                       sep='/', collapse=''),
                          data_type='Microbiome')
}


######################################
## Exploratory Analyses: Ordination ##
######################################
for( v in 1:length(exploratory_analyses) ) {
    # AMR NMDS
    meg_ordination(data_list = AMR_analytic_data,
                   data_names = AMR_analytic_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = exploratory_analyses[[v]]$exploratory_var,
                   analysis_subset=exploratory_analyses[[v]]$subsets,
                   outdir = paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                  sep='/', collapse=''),
                   data_type = 'AMR',
                   method = 'NMDS')
    
    # AMR PCA
    meg_ordination(data_list = AMR_analytic_data,
                   data_names = AMR_analytic_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = exploratory_analyses[[v]]$exploratory_var,
                   analysis_subset=exploratory_analyses[[v]]$subsets,
                   outdir = paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                  sep='/', collapse=''),
                   data_type = 'AMR',
                   method = 'PCA')
    
    # Microbiome NMDS
    meg_ordination(data_list = kraken_analytic_data,
                   data_names = kraken_analytic_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = exploratory_analyses[[v]]$exploratory_var,
                   analysis_subset=exploratory_analyses[[v]]$subsets,
                   outdir = paste(graph_output_dir, 'Microbiome', exploratory_analyses[[v]]$name,
                                  sep='/', collapse=''),
                   data_type = 'Microbiome',
                   method = 'NMDS')
    
    # Microbiome PCA
    meg_ordination(data_list = kraken_analytic_data,
                   data_names = kraken_analytic_names,
                   metadata = metadata,
                   sample_var = sample_column_id,
                   hull_var = exploratory_analyses[[v]]$exploratory_var,
                   analysis_subset=exploratory_analyses[[v]]$subsets,
                   outdir = paste(graph_output_dir, 'Microbiome', exploratory_analyses[[v]]$name,
                                  sep='/', collapse=''),
                   data_type = 'Microbiome',
                   method = 'PCA')
}


####################################
## Exploratory Analyses: Heatmaps ##
####################################

# AMR Heatmaps for each level
for( v in 1:length(exploratory_analyses) ) {
    for( l in 1:length(AMR_analytic_names) ) {
        meg_heatmap(melted_data=amr_melted_analytic,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_analyses[[v]]$exploratory_var,
                    level_var=AMR_analytic_names[l],
                    analysis_subset=exploratory_analyses[[v]]$subsets,
                    outdir=paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                 sep='/', collapse=''),
                    data_type='AMR')
    }
}

# Microbiome
for( v in 1:length(exploratory_analyses) ) {
    for( l in 1:length(kraken_analytic_names) ) {
        meg_heatmap(melted_data=kraken_melted_analytic,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_analyses[[v]]$exploratory_var,
                    level_var=kraken_analytic_names[l],
                    analysis_subset=exploratory_analyses[[v]]$subsets,
                    outdir=paste(graph_output_dir, 'Microbiome', exploratory_analyses[[v]]$name,
                                 sep='/', collapse=''),
                    data_type='Microbiome')
    }
}


####################################
## Exploratory Analyses: Barplots ##
####################################

# AMR
for( v in 1:length(exploratory_analyses) ) {
    for( l in 1:length(AMR_analytic_names) ) {
        suppressWarnings(
            meg_barplot(melted_data=amr_melted_analytic,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_analyses[[v]]$exploratory_var,
                    level_var=AMR_analytic_names[l],
                    analysis_subset=exploratory_analyses[[v]]$subsets,
                    outdir=paste(graph_output_dir, 'AMR', exploratory_analyses[[v]]$name,
                                 sep='/', collapse=''),
                    data_type='AMR')
        )
    }
}

# Microbiome
for( v in 1:length(exploratory_analyses) ) {
    for( l in 1:length(kraken_analytic_names) ) {
        suppressWarnings(
            meg_barplot(melted_data=kraken_melted_analytic,
                    metadata=metadata,
                    sample_var=sample_column_id,
                    group_var=exploratory_analyses[[v]]$exploratory_var,
                    level_var=kraken_analytic_names[l],
                    analysis_subset=exploratory_analyses[[v]]$subsets,
                    outdir=paste(graph_output_dir, 'Microbiome', exploratory_analyses[[v]]$name,
                                 sep='/', collapse=''),
                    data_type='Microbiome')
        )
    }
}


##########################
## Statistical Analyses ##
##########################
for( a in 1:length(statistical_analyses) ) {
    meg_fitZig(data_list=AMR_analytic_data,
               data_names=AMR_analytic_names,
               metadata=metadata,
               zero_mod=model.matrix(~1 + log(libSize(amr))),
               data_mod=statistical_analyses[[a]]$model_matrix,
               filter_min_threshold=0.15,
               contrast_list=statistical_analyses[[a]]$contrasts,
               random_effect_var=statistical_analyses[[a]]$random_effect,
               outdir=paste(stats_output_dir, 'AMR', statistical_analyses[[a]]$name,
                            sep='/', collapse=''),
               analysis_name=statistical_analyses[[a]]$name,
               analysis_subset=statistical_analyses[[a]]$subsets,
               data_type='AMR',
               pval=0.1,
               top_hits=1000)
    
    meg_fitZig(data_list=kraken_analytic_data,
               data_names=kraken_analytic_names,
               metadata=metadata,
               zero_mod=model.matrix(~1 + log(libSize(kraken))),
               data_mod=statistical_analyses[[a]]$model_matrix,
               filter_min_threshold=0.15,
               contrast_list=statistical_analyses[[a]]$contrasts,
               random_effect_var=statistical_analyses[[a]]$random_effect,
               outdir=paste(stats_output_dir, 'Microbiome', statistical_analyses[[a]]$name,
                            sep='/', collapse=''),
               analysis_name=statistical_analyses[[a]]$name,
               analysis_subset=statistical_analyses[[a]]$subsets,
               data_type='Microbiome',
               pval=0.1,
               top_hits=1000)
}


########################
## Output of matrices ##
########################
write.csv(make_sparse(amr_class, 'class', c('class')), 'amr_matrices/sparse_normalized/AMR_Class_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_class, 'amr_matrices/normalized/AMR_Class_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_class_raw, 'amr_matrices/raw/AMR_Class_Raw.csv', sep=',', row.names = F, col.names = T)


write.csv(make_sparse(amr_mech, 'mechanism', c('mechanism')), 'amr_matrices/sparse_normalized/AMR_Mechanism_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_mech, 'amr_matrices/normalized/AMR_Mechanism_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_mech_raw, 'amr_matrices/raw/AMR_Mechanism_Raw.csv', sep=',', row.names = F, col.names = T)

write.csv(make_sparse(amr_group, 'group', c('group')), 'amr_matrices/sparse_normalized/AMR_Group_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_group, 'amr_matrices/normalized/AMR_Group_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_mech_raw, 'amr_matrices/raw/AMR_Group_Raw.csv', sep=',', row.names = F, col.names = T)

write.csv(make_sparse(amr_norm, 'header', c('header', 'class', 'mechanism', 'group')),
          'amr_matrices/sparse_normalized/AMR_Gene_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_norm, 'amr_matrices/normalized/AMR_Gene_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_raw, 'amr_matrices/raw/AMR_Gene_Raw.csv', sep=',', row.names = F, col.names = T)


write.csv(make_sparse(kraken_domain, 'Domain', c('Domain')),
          'kraken_matrices/sparse_normalized/kraken_Domain_Sparse_Normalized.csv',
          row.names=T)
write.table(kraken_domain, 'kraken_matrices/normalized/kraken_Domain_Normalized.csv', sep=',', row.names=F, col.names=T)
write.table(kraken_domain_raw, 'kraken_matrices/raw/kraken_Domain_Raw.csv', sep=',', row.names=F, col.names=T)

write.csv(make_sparse(kraken_phylum, 'Phylum', c('Phylum')),
          'kraken_matrices/sparse_normalized/kraken_Phylum_Sparse_Normalized.csv',
          row.names=T)
write.table(kraken_phylum, 'kraken_matrices/normalized/kraken_Phylum_Normalized.csv', sep=',', row.names=F, col.names=T)
write.table(kraken_phylum_raw, 'kraken_matrices/raw/kraken_Phylum_Raw.csv', sep=',', row.names=F, col.names=T)

write.csv(make_sparse(kraken_class, 'Class', c('Class')),
          'kraken_matrices/sparse_normalized/kraken_Class_Sparse_Normalized.csv',
          row.names=T)
write.table(kraken_class, 'kraken_matrices/normalized/kraken_Class_Normalized.csv', sep=',', row.names=F, col.names=T)
write.table(kraken_class_raw, 'kraken_matrices/raw/kraken_Class_Raw.csv', sep=',', row.names=F, col.names=T)

write.csv(make_sparse(kraken_order, 'Order', c('Order')),
          'kraken_matrices/sparse_normalized/kraken_Order_Sparse_Normalized.csv',
          row.names=T)
write.table(kraken_order, 'kraken_matrices/normalized/kraken_Order_Normalized.csv', sep=',', row.names=F, col.names=T)
write.table(kraken_order_raw, 'kraken_matrices/raw/kraken_Order_Raw.csv', sep=',', row.names=F, col.names=T)

write.csv(make_sparse(kraken_family, 'Family', c('Family')),
          'kraken_matrices/sparse_normalized/kraken_Family_Sparse_Normalized.csv',
          row.names=T)
write.table(kraken_family, 'kraken_matrices/normalized/kraken_Family_Normalized.csv', sep=',', row.names=F, col.names=T)
write.table(kraken_family_raw, 'kraken_matrices/raw/kraken_Family_Raw.csv', sep=',', row.names=F, col.names=T)

write.csv(make_sparse(kraken_genus, 'Genus', c('Genus')),
          'kraken_matrices/sparse_normalized/kraken_Genus_Sparse_Normalized.csv',
          row.names=T)
write.table(kraken_genus, 'kraken_matrices/normalized/kraken_Genus_Normalized.csv', sep=',', row.names=F, col.names=T)
write.table(kraken_genus_raw, 'kraken_matrices/raw/kraken_Genus_Raw.csv', sep=',', row.names=F, col.names=T)

write.csv(make_sparse(kraken_species, 'Species', c('Species')),
          'kraken_matrices/sparse_normalized/kraken_Species_Sparse_Normalized.csv',
          row.names=T)
write.table(kraken_species, 'kraken_matrices/normalized/kraken_Species_Normalized.csv', sep=',', row.names=F, col.names=T)
write.table(kraken_species_raw, 'kraken_matrices/raw/kraken_Species_Raw.csv', sep=',', row.names=F, col.names=T)




