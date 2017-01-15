###############
## Functions ##
###############
## Utility functions that can be optionally used later for producing
## graphs and performing reshaping operations

# Call variables from parent scope
`..` <- function (..., .env = sys.parent(2)) {
    get(deparse(substitute(...)), env = .env)
}

# Misc reshape function for data table
melt_dt <- function(D, level_id) {
    temp <- melt(D, variable.name='Sample', value.name='Normalized_Count')
    names(temp) <- c('Name', 'ID', 'Normalized_Count')
    temp <- data.table(cbind(rep(level_id, nrow(temp)), temp))
    names(temp)[1] <- 'Level_ID'
    return(temp)
}

# Calculate the points for convex hulls around data for ordination plots
meg_find_hulls <- function(x) x[chull(x$Ord1, x$Ord2),]

# Function that returns species count, rarefied species count, and alpha diversity measures
# for each sample in the m x n matrix, m = features, n = samples
alpha_rarefaction <- function(X, minlevel, method='invsimpson') {
    S <- specnumber(X, MARGIN=2)
    raremax <- min(colSums(X))
    if( raremax < minlevel ) raremax <- minlevel
    Srare <- rarefy(X, raremax, MARGIN=2)
    Xrare <- t(rrarefy(t(X), raremax))
    alphadiv <- diversity(Xrare, index=method, MARGIN=2)
    return(list(raw_species_abundance=S,
                rarefied_species_abundance=Srare,
                rarefied_data=Xrare,
                alphadiv=alphadiv))
}


heatmap_select_top_counts <- function(X, group_var, sample_var, n) {
    return(X[, tail(.SD, n), by=c(group_var, sample_var)])
}

bar_select_top_counts <- function(X, group_var, n) {
    return(X[, tail(.SD, n), by=group_var])
}


make_sparse <- function(df, rownames, excludes, filter_min_threshold=0.15) {
    local_df <- df[, .SD, .SDcols=!excludes]
    filter_threshold <- quantile(rowSums(local_df), 0.15)
    if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
    chosen <- which(rowSums(local_df) >= filter_threshold )
    ret <- as.data.frame(local_df[chosen, ])
    rownames(ret) <- df[[rownames]][chosen]
    return(ret)
}

data_subset <- function(data_obj, subsets) {
  local_meta <- data.table(pData(data_obj))
  local_subset <- c()          
    for( c in 1:length(subsets) ) {
      
      conditional_terms <- unlist(strsplit(subsets[[c]], ' '))
      conditional_string <- paste('local_meta[[\'', conditional_terms[1],
                                  '\']] ', conditional_terms[2],
                                  ' \'', conditional_terms[3], '\'',
                                  sep='', collapse='')
      if(length(local_subset) > 0) {
        local_subset <- intersect(local_subset, which(eval(parse(text=conditional_string))))
      }
      else {
        local_subset <- which(eval(parse(text=conditional_string)))
      }
  }
  return(data_obj[, local_subset])
}


data_subset_long <- function(data_obj, subsets) {
    local_subset <- c()          
    for( c in 1:length(subsets) ) {
        
        conditional_terms <- unlist(strsplit(subsets[[c]], ' '))
        conditional_string <- paste('data_obj[[\'', conditional_terms[1],
                                    '\']] ', conditional_terms[2],
                                    ' \'', conditional_terms[3], '\'',
                                    sep='', collapse='')
        if(length(local_subset) > 0) {
            local_subset <- intersect(local_subset, which(eval(parse(text=conditional_string))))
        }
        else {
            local_subset <- which(eval(parse(text=conditional_string)))
        }
    }
    return(data_obj[local_subset, ])
}


# Function for computing ordination plots with convex hulls
# using ggplot2.  You will have to specify the facet variable,
# i.e. which experimental design variable determines the grouping of points.
#
# Function arguments:
#   data_list: a list containing the MRexperiment objects for all levels
#   data_names: a character vector of level names for each MRexperiment in data_list
#   metadata: a data.table of metadata for each sample
#   sample_var: the column name in metadata that specifies the sample IDs (must match MRexperiment columns)
#   hull_var: the metadata column name by which to group data points
#   outdir: the file path of the directory for output of files
#   data_type: the data type being computed on, e.g. AMR or Microbiome
#   method: choice of 'NMDS' or 'PCA' for method of ordination
#
#   outputs:  
meg_ordination <- function(data_list,
                           data_names,
                           metadata,
                           sample_var,
                           hull_var,
                           analysis_subset,
                           outdir,
                           data_type,
                           method='NMDS') {
    all_ord <- data.table(ID=character(),
                          Level_ID=character(),
                          Ord1=numeric(),
                          Ord2=numeric(),
                          Group_Var=character())
    all_hulls <- data.table(Group_Var=character(),
                            Ord1=numeric(),
                            Ord2=numeric(),
                            Level_ID=character())
    setkey(all_ord, ID)
    
    local_obj <- data_list
    for( l in 1:length(local_obj) ) {
      
        if(length(analysis_subset) > 0) {
          local_obj[[l]] <- data_subset(local_obj[[l]], analysis_subset)
        }
        local_meta <- metadata
        local_meta <- local_meta[local_meta[[sample_var]] %in% colnames(MRcounts(local_obj[[l]])), ]
        
        # Transpose the matrix for NMDS (groups are now in rows and features in columns)
        t_data <- t(MRcounts(local_obj[[l]]))
        local_meta <- local_meta[rowSums(t_data) > 0, ]
        t_data <- t_data[rowSums(t_data) > 0, ]
        t_data <- t_data[which(!is.na(local_meta[[sample_var]]) & local_meta[[sample_var]] != 'NA'), colSums(t_data) > 0]
        
        if( method == 'NMDS' ) {
            # Set parallel to whatever your computer can support in terms of CPU count
            ord.res <- metaMDS(t_data, autotransform=F, parallel=7, trymax=49)
            ord_points <- data.table(ord.res$points)
            names(ord_points) <- c('Ord1', 'Ord2')
            ord_points[, ID :=( rownames(ord.res$points) )]
        }
        else if( method == 'PCA' ) {
            ord.res <- prcomp(t_data, center=T, scale=T)
            
            # Format to include metadata for ggplot2
            ord_points <- data.table(ord.res$x[, 1:2])
            names(ord_points) <- c('Ord1', 'Ord2')
            ord_points[, ID :=( rownames(ord.res$x) )]
        }
        else {
            stop('method must be either NMDS or PCA')
        }
        
        ord_points[, Level_ID :=( rep(data_names[l], nrow(ord_points)) )]
        setkey(ord_points, ID)
        ord_points <- metadata[ord_points]
        ord_points <- ord_points[, .SD, .SDcols=c(sample_var, hull_var, 'Ord1', 'Ord2', 'Level_ID')]
        names(ord_points)[2] <- 'Group_Var'
        all_ord <- rbind(all_ord, ord_points)
        
        hulls <- ord_points[, meg_find_hulls(.SD), .SDcols=c('Ord1', 'Ord2'), by=Group_Var]
        hulls[, Level_ID :=( rep(data_names[l], nrow(hulls)) )]
        hulls[, .SD, .SDcols=!'Level_ID']
        all_hulls <- rbind(all_hulls, hulls)
        
        # Plot graphs with convex hulls
        g_ord <- ggplot(data=ord_points, aes(Ord1, Ord2, color=Group_Var, fill=Group_Var)) +
            geom_point(size=2.5) + geom_polygon(data=hulls, aes(x=Ord1, y=Ord2, color=Group_Var, fill=Group_Var),
                                                alpha=0.2, show.legend=F)
        g_ord <- g_ord +
            ggtitle(paste(method, ' for ', data_type, ' ', 'by ', hull_var, '\nAnnotation Level: ',
                          data_names[l],
                          sep='',
                          collapse='')) +
            labs(color=hull_var) +
            guides(fill=F) +
            theme(strip.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.text.x=element_blank(),
                  axis.title.x=element_text(size=24),
                  axis.title.y=element_text(size=24, hjust=0.5),
                  #legend.position="right",
                  legend.title=element_text(size=20),
                  legend.text=element_text(size=18),
                  plot.title=element_text(size=24, hjust=0.5))
        if( method == 'NMDS' ) {
            g_ord <- g_ord + xlab('MDS1') + ylab('MDS2')
        }
        else if( method == 'PCA' ) {
            g_ord <- g_ord + xlab('PC1') + ylab('PC2')
        }
        
        # Open the graphics device at the specified location and figure size
        png(filename=paste(outdir, '/', method, '_', hull_var, '_',
                           data_names[l], '.png', sep='', collapse=''), 
            width=1024, height=768)
        
        print(g_ord)
        
        # Turn off graphics device to save the graphic
        dev.off()
    }
    
    all_ord <- within(all_ord, Level_ID
                      <- factor(Level_ID, levels=data_names,
                                ordered=T))
    all_hulls <- within(all_hulls, Level_ID
                        <- factor(Level_ID, levels=data_names,
                                  ordered=T))

    g_all_ord <- ggplot(data=all_ord, aes(Ord1, Ord2, color=Group_Var, fill=Group_Var)) +
        geom_point(size=3) + 
        geom_polygon(data=all_hulls, aes(x=Ord1, y=Ord2, color=Group_Var, fill=Group_Var),
                     alpha=0.2, show.legend=F) +
        facet_wrap(~Level_ID, nrow=2)
    g_all_ord <- g_all_ord +
        ggtitle(paste(method, ' for ', data_type, ' by ', hull_var, sep='', collapse='')) +
        labs(color=hull_var) +
        guides(fill=F) +
        theme(strip.text.x=element_text(size=26),
              axis.text.y=element_blank(),
              axis.text.x=element_blank(),
              axis.title.x=element_text(size=26),
              axis.title.y=element_text(size=26, hjust=0.5),
              #legend.position="right",
              legend.title=element_text(size=24, hjust=0.5),
              legend.text=element_text(size=20),
              plot.title=element_text(size=30, hjust=0.5))
    
    png(filename=paste(outdir, '/', method, '_', hull_var, '_',
                       'AllLevels.png', sep='', collapse=''),
        width=1024, height=768)
    print(g_all_ord)
    dev.off()
}


meg_heatmap <- function(melted_data,
                        metadata,
                        sample_var,
                        group_var,
                        level_var,
                        analysis_subset,
                        outdir,
                        data_type) {
    tile_subset <- melted_data[Level_ID == level_var, ]
    colnames(tile_subset)[colnames(tile_subset) == 'ID'] <- sample_var
    setkeyv(tile_subset, sample_var)
    
    tile_subset <- metadata[tile_subset]
    tile_subset <- tile_subset[!is.na(tile_subset[[group_var]]), ]
    
    if(length(analysis_subset) > 0) {
        tile_subset <- data_subset_long(tile_subset, analysis_subset)
    }
    
    sample_order <- unique(tile_subset[order(group_var), sample_var])
    tile_subset <- within(tile_subset, sample_var
                                <- factor(sample_var,
                                          levels=sample_order,
                                          ordered=T))
    
    setkey(tile_subset, Normalized_Count)
    tile_subset <- tile_subset[, sum(Normalized_Count),
                                by=c(group_var, sample_var, 'Name')]
    names(tile_subset)[length(names(tile_subset))] <- 'Normalized_Count'
    
    numselect <- 20
    tile_names <- heatmap_select_top_counts(tile_subset, group_var,
                                     sample_var, numselect)
    name_count <- length(unique(tile_names$Name))
    while(name_count > 20) {
        numselect <- numselect - 1
        tile_names <- heatmap_select_top_counts(tile_subset, group_var,
                                         sample_var, numselect)
        name_count <- length(unique(tile_names$Name))
    }
    
    tile_subset <- tile_subset[Name %in% tile_names$Name, ]
    
    tile <- ggplot(tile_subset, aes_string(x=sample_var, y='Name')) +
        geom_tile(aes(fill=log2(Normalized_Count+1))) +
        facet_wrap(as.formula(paste('~', group_var)), switch='x', scales = 'free_x', nrow = 1) +
        theme(panel.background=element_rect(fill="black", colour="black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text.x=element_text(size=20),
              axis.text.y=element_text(size=20),
              axis.text.x=element_blank(),
              axis.title.x=element_text(size=20),
              axis.title.y=element_blank(),
              legend.position="bottom",
              legend.title=element_text(size=20),
              panel.margin=unit(0.1, "lines"),
              plot.title=element_text(size=24, hjust=0.5),
              plot.margin=unit(c(0,0,0,0), "cm")) +
        xlab(paste('Samples by ', group_var, sep='', collapse='')) +
        scale_fill_gradient(low="black", high="cyan") +
        labs(fill= 'Log2 Normalized Count') +
        ggtitle(paste(data_type, ' ', level_var, ' Normalized Counts by ', group_var, '\n',
                      sep='', collapse=''))
    png(filename=paste(outdir, '/', level_var, '_', group_var, '_',
                       'Heatmap.png', sep='', collapse=''), width=1400, height=700)
    print(tile)
    dev.off()
}


meg_alpha_rarefaction <- function(data_list,
                                  data_names,
                                  metadata,
                                  sample_var,
                                  group_var,
                                  analysis_subset,
                                  outdir,
                                  data_type) {
    all_alphadiv <- data.table(ID=character(),
                               Level=character(),
                               Value=numeric())
    names(all_alphadiv)[1] <- sample_var
    
    all_species_raw <- data.table(ID=character(),
                                  Level=character(),
                                  Value=numeric())
    names(all_species_raw)[1] <- sample_var
    
    all_species_rare <- data.table(ID=character(),
                                   Level=character(),
                                   Value=numeric())
    names(all_species_rare)[1] <- sample_var
    
    local_data <- data_list
    
    for( l in 1:length(local_data) ) {
        
        if(length(analysis_subset) > 0) {
            local_data[[l]] <- data_subset(local_data[[l]], analysis_subset)
        }
        
        sample_counts <- colSums(MRcounts(local_data[[l]]))
        if(min(sample_counts) == 0) {
            non_zero_sample <- min(sample_counts[sample_counts > 0])
        }
        else {
            non_zero_sample <- 0
        }
        
        local_obj <- alpha_rarefaction(MRcounts(local_data[[l]]), minlevel = non_zero_sample)
        temp <- data.table(ID=names(local_obj$alphadiv),
                           Level=rep(data_names[l],
                                     length(local_obj$alphadiv)),
                           Value=as.numeric(local_obj$alphadiv))
        names(temp)[1] <- sample_var
        all_alphadiv <- rbind(all_alphadiv, temp)
        
        temp <- data.table(ID=names(local_obj$raw_species_abundance),
                   Level=rep(data_names[l],
                             length(local_obj$raw_species_abundance)),
                   Value=as.numeric(local_obj$raw_species_abundance))
        names(temp)[1] <- sample_var
        all_species_raw <- rbind(all_species_raw, temp)
        
        temp <- data.table(ID=names(local_obj$rarefied_species_abundance),
                   Level=rep(data_names[l],
                             length(local_obj$rarefied_species_abundance)),
                   Value=as.numeric(local_obj$rarefied_species_abundance))
        names(temp)[1] <- sample_var
        all_species_rare <- rbind(all_species_rare, temp)
    }
    
    all_alphadiv <- within(all_alphadiv, Level
                           <- factor(Level, levels=data_names,
                                     ordered=T))
    all_species_raw <- within(all_species_raw, Level
                              <- factor(Level, levels=data_names,
                                        ordered=T))
    all_species_rare <- within(all_species_rare, Level
                               <- factor(Level, levels=data_names,
                                         ordered=T))
    setkeyv(all_alphadiv, sample_var)
    setkeyv(all_species_raw, sample_var)
    setkeyv(all_species_rare, sample_var)
    
    all_alphadiv <- metadata[all_alphadiv]
    all_species_raw <- metadata[all_species_raw]
    all_species_rare <- metadata[all_species_rare]
    
    all_alphadiv <- all_alphadiv[!is.na(all_alphadiv[[group_var]]), ]
    all_species_raw <- all_species_raw[!is.na(all_species_raw[[group_var]]), ]
    all_species_rare <- all_species_rare[!is.na(all_species_rare[[group_var]]), ]
        
    alphadiv_type_sums <- all_alphadiv[Level==data_names[2], median(Value), by=group_var]
    alphadiv_value_labels <- as.character(alphadiv_type_sums[[group_var]][order(alphadiv_type_sums$V1, decreasing=T)])
    
    species_raw_type_sums <- all_species_raw[Level==data_names[2], median(Value), by=group_var]
    species_raw_value_labels <- as.character(species_raw_type_sums[[group_var]][order(species_raw_type_sums$V1,
                                                                                decreasing=T)])
    species_rare_type_sums <- all_species_rare[Level==data_names[2], median(Value), by=group_var]
    species_rare_value_labels <- as.character(species_rare_type_sums[[group_var]][order(species_rare_type_sums$V1,
                                                                                       decreasing=T)])
    
    all_alphadiv[[group_var]] <- factor(all_alphadiv[[group_var]],
                                     levels=alphadiv_value_labels, ordered=T)
    
    all_species_raw[[group_var]] <- factor(all_species_raw[[group_var]],
                              levels=species_raw_value_labels, ordered=T)
    
    all_species_rare[[group_var]] <- factor(all_species_rare[[group_var]],
                               levels=species_rare_value_labels, ordered=T)
    
    png(filename=paste(outdir, '/', data_type, '_alphadiversity_by_', group_var, '.png',
                       sep='', collapse=''),
        width=1024, height=768)
    g_alphadiv <- ggplot(data=all_alphadiv, aes_string(x=group_var,
                             y='Value',
                             color=group_var)) +
        geom_boxplot(size=1) + 
        facet_wrap(~Level, nrow=2, scales='free_y')
    g_alphadiv <- g_alphadiv +
        ggtitle(paste('Alpha Diversity by ', group_var, ' for Rarefied data\nInverse Simpson Index',
                      sep='', collapse='')) +
        ylab('Inverse Simpson\'s Index\n') +
        xlab(paste('\n', group_var, sep='', collapse='')) +
        theme(strip.text.x=element_text(size=26),
              axis.text.y=element_text(size=20),
              axis.text.x=element_blank(),
              axis.title.x=element_text(size=26),
              axis.title.y=element_text(size=26, vjust=1),
              #legend.position="right",
              legend.title=element_text(size=24, vjust=1),
              legend.text=element_text(size=20),
              plot.title=element_text(size=30, hjust=0.5))
    print(g_alphadiv)
    dev.off()
    
    png(filename=paste(outdir, '/', data_type, '_raw_richness_by_', group_var, '.png',
                       sep='', collapse=''),
        width=1024, height=768)
    g_sraw <- ggplot(data=all_species_raw, aes_string(group_var, 'Value', color=group_var)) +
        geom_boxplot(size=1) + 
        facet_wrap(~Level, nrow=2, scales='free_y')
    g_sraw <- g_sraw +
        ggtitle(paste('Species Richness by ', group_var, ' for Raw data\nInverse Simpson Index',
                      sep='', collapse='')) +
        ylab('Unique Species\n') +
        xlab(paste('\n', group_var, sep='', collapse='')) +
        theme(strip.text.x=element_text(size=26),
              axis.text.y=element_text(size=20),
              axis.text.x=element_blank(),
              axis.title.x=element_text(size=26),
              axis.title.y=element_text(size=26, vjust=1),
              #legend.position="right",
              legend.title=element_text(size=24, vjust=1),
              legend.text=element_text(size=20),
              plot.title=element_text(size=30, hjust=0.5))
    print(g_sraw)
    dev.off()
    
    png(filename=paste(outdir, '/', data_type, '_rarefied_richness_by_', group_var, '.png',
                       sep='', collapse=''),
        width=1024, height=768)
    g_srare <- ggplot(data=all_species_rare, aes_string(group_var, 'Value', color=group_var)) +
        geom_boxplot(size=1) + 
        facet_wrap(~Level, nrow=2, scales='free_y')
    g_srare <- g_srare +
        ggtitle(paste('Species Richness by ', group_var, ' for Rarefied data\nInverse Simpson Index',
                      sep='', collapse='')) +
        ylab('Unique Species\n') +
        xlab(paste('\n', group_var, sep='', collapse='')) +
        theme(strip.text.x=element_text(size=26),
              axis.text.y=element_text(size=20),
              axis.text.x=element_blank(),
              axis.title.x=element_text(size=26),
              axis.title.y=element_text(size=26, vjust=1),
              #legend.position="right",
              legend.title=element_text(size=24, vjust=1),
              legend.text=element_text(size=20),
              plot.title=element_text(size=30, hjust=0.5))
    print(g_srare)
    dev.off()
}


meg_barplot <- function(melted_data,
                        metadata,
                        sample_var,
                        group_var,
                        level_var,
                        analysis_subset,
                        outdir,
                        data_type) {
    setkeyv(melted_data, sample_var)
    melted_data <- metadata[melted_data]
    
    if(length(analysis_subset) > 0) {
        bar_subset <- data_subset_long(melted_data, analysis_subset)
    }
    else {
        bar_subset <- melted_data
    }
    
    bar_subset <- data.table(bar_subset[Level_ID == level_var &
                                                 !is.na(bar_subset[[group_var]]),
                                                    .SD,
                                                    .SDcols=c(sample_var,
                                                              group_var,
                                                              'Name',
                                                              'Normalized_Count')])
    setkeyv(bar_subset, sample_var)
    bar_subset <- metadata[bar_subset]
    bar_subset[[sample_var]] <- factor(bar_subset[[sample_var]],
                                          levels=unique(bar_subset[[sample_var]][order(bar_subset[[group_var]])]),
                                          ordered=T)
    
    setkey(bar_subset, Normalized_Count)
    bar_subset[, sample_number:=(length(unique(bar_subset[[sample_var]]))), by=c(group_var, 'Name')]
    bar_subset <- unique(bar_subset[, sum(Normalized_Count) / sample_number,
                                                  by=c(group_var, 'Name')])
    
    numselect <- 11
    bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
    name_count <- length(unique(bar_names$Name))
    while(name_count > 11) {
        numselect <- numselect - 1
        bar_names <- bar_select_top_counts(bar_subset, group_var, numselect)
        name_count <- length(unique(bar_names$Name))
    }
    
    bar_subset <- bar_subset[Name %in% bar_names$Name, ]
    
    
    names(bar_subset)[length(names(bar_subset))] <- 'Normalized_Count'
    source_sums <- tapply(bar_subset[['Normalized_Count']],
                          bar_subset[[group_var]], sum)
    source_labels <- names(source_sums)[order(source_sums, decreasing=T)]
    bar_subset[[group_var]] <- factor(bar_subset[[group_var]],
                                          levels=source_labels, ordered=T)
    
    
    meg_bar <- ggplot(bar_subset, aes_string(x=group_var, y='Normalized_Count', fill='Name')) +
        geom_bar(stat='identity') + 
        scale_fill_brewer(palette="Spectral") +
        theme(strip.text.x=element_text(size=18),
              axis.text.y=element_text(size=20),
              axis.text.x=element_text(size=22, vjust=1, hjust=1, angle=33),
              axis.title.x=element_text(size=24),
              axis.title.y=element_text(size=24),
              legend.title=element_text(size=20),
              legend.text=element_text(size=18),
              plot.title=element_text(size=26, hjust=0.25)) +
        xlab(group_var) +
        ylab('Mean of Normalized Count\n') +
        ggtitle(paste('Mean ', data_type, ' ', level_var, ' Normalized Count by ', group_var, '\n',
                      sep='', collapse=''))
    png(filename=paste(outdir, '/', data_type, '_', level_var, '_BarPlot_by_', group_var, '.png',
                       sep='', collapse=''), width=1024, height=768)
    print(meg_bar)
    dev.off()
}


meg_fitZig <- function(data_list,
                       data_names,
                       metadata,
                       zero_mod,
                       data_mod,
                       filter_min_threshold,
                       contrast_list,
                       random_effect_var,
                       outdir,
                       analysis_name,
                       analysis_subset,
                       data_type,
                       pval=0.1,
                       top_hits=200) {
    settings <- zigControl(maxit=50, verbose=F)
    
    local_obj <- data_list
    res <- list()
    for( l in 1:length(local_obj) ) {
        
        filter_threshold <- quantile(rowSums(MRcounts(local_obj[[l]])), 0.15)
        if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
        local_obj[[l]] <- local_obj[[l]][which(rowSums(MRcounts(local_obj[[l]])) >= filter_threshold ), ]
        
        if(length(analysis_subset) > 0) {
            local_obj[[l]] <- data_subset(local_obj[[l]], analysis_subset)
        }
        
        
        col_selection <- as.integer(which(colSums(MRcounts(local_obj[[l]]) > 0) > 1))
        local_obj[[l]] <- local_obj[[l]][, col_selection]

        mod_select <- model.matrix(eval(parse(text=data_mod)), data=pData(local_obj[[l]]))
        zero_mod_select <- zero_mod[col_selection, ]
        
        cumNorm(local_obj[[l]])  # This is a placeholder for metagenomeSeq; we don't actually use these values
        
        tryCatch(
            {
                if( is.na(random_effect_var) ) {
                    res[[l]] <- fitZig(obj=local_obj[[l]],
                                       mod=mod_select,
                                       zeroMod=zero_mod_select,
                                       control=settings,
                                       useCSSoffset=F)
                }
                else {
                    res[[l]] <- fitZig(obj=local_obj[[l]],
                                       mod=mod_select,
                                       zeroMod=zero_mod_select,
                                       control=settings,
                                       useCSSoffset=F,
                                       useMixedModel=T,
                                       block=pData(local_obj[[l]])[, random_effect_var])
                }
            },
            error=function(e) {
                print(paste('Model failed to converge for ', data_type, ' ', data_names[l], ' ', analysis_name,
                      sep='', collapse=''))
            },
            finally={
                if( length(res) != l ) {
                    next
                }
            }
        )

        local_contrasts <- contrast_list
        local_contrasts[[length(local_contrasts)+1]] <- res[[l]]$fit$design
        names(local_contrasts)[length(local_contrasts)] <- 'levels'
        
        contrast_matrix <- do.call(makeContrasts, local_contrasts)
        colnames(contrast_matrix) <- make.names(contrast_list)
        
        contrast_fit <- contrasts.fit(res[[l]]$fit, contrast_matrix)        
        contrast_fit <- eBayes(contrast_fit)
        
        stats_results <- data.table(
            Node.Name=character(),
            Contrast=character(),
            logFC=numeric(),
            CI.L=numeric(),
            CI.R=numeric(),
            AveExpr=numeric(),
            t=numeric(),
            P.Value=numeric(),
            adj.P.Val=numeric(),
            B=numeric()
        )
        
        for( c in 1:ncol(contrast_fit$contrasts) ) {
            tophits <- topTable(contrast_fit, p.value=pval, confint=T,
                                number=top_hits, sort.by='AveExpr', coef=c)
            
            if( nrow(tophits) > 0) {
                temp_res <- data.table(
                    Node.Name=rownames(tophits),
                    Contrast=rep(colnames(contrast_fit$contrasts)[c], nrow(tophits))
                )
                temp_res <- cbind(temp_res, tophits)
                stats_results <- rbind(stats_results, temp_res)
            }
            else {
                print(paste('No significant results for', data_type,
                            data_names[l], analysis_name,
                            colnames(contrast_fit$contrasts)[c],
                            sep=' ', collapse=''))
            }
        }
                
        if( nrow(stats_results) > 0 ) {
            write.csv(stats_results,
                      file=paste(outdir, '/', analysis_name, '_', data_type, '_',
                                                data_names[l], '_',
                                                contrast_list[1], '_Model_Contrasts.csv',
                                                sep='', collapse=''),
                      quote=F, row.names=F)
        }
    }
}






