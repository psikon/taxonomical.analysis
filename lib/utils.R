#' rm.taxa
#'
#' remove a taxa with a specific tax id from phyloseq object 
#' 
#'@param phyloseq   phyloseq object
#'@param taxa       taxa to remove
#'
#'@return phyloseq object
#'@export
#'
rm.taxa <- function(phyloseq, 
                    taxa) {
    phyloseq <- prune_species(!grepl(paste0("\\<", taxa, "\\>"),
                                     labels(otu_table(phyloseq))[[1]]), 
                              phyloseq)
    return(phyloseq)
}

#' rm.underscore
#'
#' remove underscore syntax from tax_table taxonomies
#' 
#'@param phyloseq   phyloseq object
#'
#'@return phyloseq object
#'@export
#'
rm.underscore <- function(phyloseq) {
    # extract the tax_table
    data <- tax_table(phyloseq)
    if(any(grepl("__", data))) {
        # use only the contents of string after the 'x__' syntax
        data <- substr(data, 4, length(data))
        # subset empty strings with a userdefined pattern
        data[which(data == "")] <- "not classified"
    }
    # add the modified tax_table to phyloseq object
    tax_table(phyloseq) <- data
    return(phyloseq)
}

#' get.phyloLevels
#'
#' count levels of taxomomies
#' 
#'@param phyloseq   phyloseq object
#'@param sample     sample for count
#'@param ranks      taxonomical ranks
#'@param absolute   absolute or percantage values
#'
#'@return data.frame
#'@keywords internal
#'
# get the number of levels of a specific OTU
get.phyloLevels <- function(phyloseq, sample, ranks, absolute) {
    # get all otus with hits in otu_table from phyloseq object
    id <- which(get_taxa(phyloseq, sample) > 0)
    if (length(id) == 0) {
        # when no hits 
        otus <- NULL
    } else {
        # seperate taxomomies of otus 
        otus <- tax_table(phyloseq)[id]
    }
    # get number of each taxonomic level 
    res <-  sapply(ranks, function(x){
        get.levelCount(otus, x)
    })
    if(!absolute) {
        # convert values to percent 
        res <- as.percent(res)
    }
    return(res)
}

#' get.levelCount
#'
#' count the taxonomical levels of of a OTU
#' 
#'@param tax_table   phyloseq tax_table
#'@param level       specific taxonomical level
#'
#'@return level
#'@keywords internal
#'
get.levelCount <- function(tax_table, level) {
    if(is.null(tax_table)) {
        # no otus found
        x <- 0 
    } else {
        # seperate x__ string from taxonomies
        x <- substring(as.character(tax_table[, level]), 4)
        # convert empty levels to NA
        x[x == ""] <- NA
        # get length of levels without NA
        x <- length(compactNA(x))
    }
    return(x)
}

#' as.percent
#'
#' convert a value to percantage
#' 
#'@param x value
#'
#'@return percentage
#'@keywords internal
#'
as.percent <- function(x) {
    # w/p = g/100
    # first value in df is 100%
    g <- x[1]
    # get percentage for other values
    y <- sapply(x, function(w){
        round(w * 100 / g, 2)
    })
    names(y) <- names(x)
    return(y)
}

#' get.free
#'
#' get only sample with HoldingCondition free_living
#' 
#'@param phyloseq   phyloseq object
#'
#'@return phyloseq object
#'@keywords internal
#'
get.free <- function(phyloseq) {
    return(subset_samples(phyloseq, 
                          HoldingCondition == "free living"))
}

#' get.aqua
#'
#' get only sample with HoldingCondition mariculture
#' 
#'@param phyloseq   phyloseq object
#'
#'@return phyloseq object
#'@keywords internal
#'
get.aqua <- function(phyloseq) {
    return(subset_samples(phyloseq, 
                          HoldingCondition == "mariculture"))
}

#' last.rank
#'
#' get the last rank of a linage from tax_table
#' 
#'@param x   tax_table
#'
#'@return last rank
#'@keywords internal
#'
last.rank <- function(x) {
    x <- tail(colnames(x[, sapply(strsplit(x, "*__"), function(y){
        length(y) > 1
    })]), n = 1)
    return(x)
}
#' last.taxa
#'
#' get only the last taxa from tax_table
#' 
#'@param x  tax_table
#'
#'@return last taxa
#'@keywords internal
#'
last.taxa <- function(x) {
    rank <- tail(colnames(x[, sapply(strsplit(x, "*__"), function(y){
        length(y)>1
    })]), n = 1)
    x <- unlist(str_split(x[,rank],"__"))[2]
    return(x)
}

#' get.richness theme
#'
#' return a predefined theme for richness plot
#'
#'@return ggplot2 theme
#'@keywords internal
#'
get.richnessTheme <- function() {
    ggtheme.alpha <- theme_bw() + 
        theme(legend.position = "top",
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 12),
              axis.text = element_text(size = 8),
              axis.title = element_text(size = 12),
              strip.text = element_text(size = 12),
              plot.margin = unit(c(0.025,0.025,.025,0.025), "npc"),
              axis.text.x = element_text(face = "italic", 
                                         size=rel(1), 
                                         angle = 30, 
                                         hjust = 1, 
                                         vjust = 1),
              axis.title = element_text(size=rel(1), 
                                        lineheight=1.5),
              legend.key = element_rect(colour = "white")
              )
    return(ggtheme.alpha)
}

#' countArtifacts
#'
#' count the number of Artifacts in phyloseq object. 
#' An artifact is defined by mapping only to the 
#' superkingdom level in the linage
#' 
#'@param phyloseq   phyloseq object
#'@param sample     specific sample
#'
#'@return number of Artifacts
#'@keywords internal
#'
countArtifacts <- function(phyloseq, sample) {
    # extract all bacterial taxa
    with <- subset_taxa(phyloseq, superkingdom == 'k__Bacteria')
    # remove taxa with maximal mapping to superkingdom
    without <- rm.taxa(with, "2")
    # calculate number of atrifacts and return them
    return(countReads(with,
                      as.character(sample)) - countReads(without, 
                                                         as.character(sample)))
}

#' countReads
#'
#' count the total number of taxonomical hits in a sample
#' 
#'@param phyloseq   phyloseq object
#'@param sample     specific sample
#'
#'@return number of taxonomical hits
#'@keywords internal
#'
countReads <- function(phyloseq, sample) {
    value <- sum(count(get_taxa(phyloseq, as.character(sample)))$x * count(get_taxa(phyloseq, as.character(sample)))$freq)
    return(value)
}
#' getOTUs
#'
#' count the number of OTUs in a sample
#' 
#'@param phyloseq   phyloseq object
#'@param sample     specific sample
#'
#'@return number OTUs
#'@keywords internal
#'
getOTUs <- function(phyloseq, sample) {
    return(length(which(otu_table(phyloseq)[, as.character(sample)] > 0) == TRUE))
}

#' traverse
#'
#' recursion function for newick string creation
#' 
#'@description from: Martin Turjak
#'             at: http://stackoverflow.com/questions/15343338/        
#'
#'@param a   
#'@param i     
#'@param innerl
#'
#'@return newick string
#'@keywords internal
#'
traverse <- function(a, i, innerl){
    if(i < (ncol(df))){
        alevelinner <- as.character(unique(df[which(as.character(df[ , i]) == a), i+1]))
        desc <- NULL
        if(length(alevelinner) == 1) (newickout <- traverse(alevelinner, i+1, innerl))
        else {
            for(b in alevelinner) desc <- c(desc, traverse(b, i+1, innerl))
            il <- NULL; if(innerl == TRUE) il <- a
            (newickout <- paste("(", 
                                paste(desc, collapse = ","), 
                                ")", il, sep = ""))
        }
    }
    else { (newickout <- a) }
}

#' df2newick
#'
#' data.frame to newick function
#' 
#'@description from: Martin Turjak
#'             at: http://stackoverflow.com/questions/15343338/        
#'
#'@param df   
#'@param innerlabel     
#'
#'@return newick string
#'@keywords internal
#'
df2newick <- function(df, innerlabel = FALSE) {
    # determine root
    alevel <- as.character(unique(df[ , 1]))
    newick <- NULL
    # traverse throug data.frame and create newick string
    for(x in alevel) newick <- c(newick, traverse(x, 1, innerlabel))
    # add outer strings
    (newick <- paste("(", 
                     paste(newick, collapse = ","), ");", 
                     sep = ""))
    return(newick)
}
#' psmelt
#'
#' overwrite to fix bug in plot_bar       
#'
#'@param physeq phyloseq object       
#'
#'@return mdf
#'@keywords internal
#'
psmelt <- function(physeq) 
{
    if (!inherits(physeq, "phyloseq")) {
        rankNames = NULL
        sampleVars = NULL
    }
    else {
        rankNames = rank_names(physeq, FALSE)
        sampleVars = sample_variables(physeq, FALSE)
    }
    reservedVarnames = c("Sample", "Abundance", "OTU")
    type1aconflict = intersect(reservedVarnames, sampleVars)
    if (length(type1aconflict) > 0) {
        wh1a = which(sampleVars %in% type1aconflict)
        new1a = paste0("sample_", sampleVars[wh1a])
        warning("The sample variables: \n", paste(sampleVars[wh1a], 
                                                  collapse = ", "), "\n have been renamed to: \n", 
                paste0(new1a, collapse = ", "), "\n", "to avoid conflicts with special phyloseq plot attribute names.")
        colnames(sample_data(physeq))[wh1a] <- new1a
    }
    type1bconflict = intersect(reservedVarnames, rankNames)
    if (length(type1bconflict) > 0) {
        wh1b = which(rankNames %in% type1bconflict)
        new1b = paste0("taxa_", rankNames[wh1b])
        warning("The rank names: \n", paste(rankNames[wh1b], 
                                            collapse = ", "), "\n have been renamed to: \n", 
                paste0(new1b, collapse = ", "), "\n", "to avoid conflicts with special phyloseq plot attribute names.")
        colnames(tax_table(physeq))[wh1b] <- new1b
    }
    type2conflict = intersect(sampleVars, rankNames)
    if (length(type2conflict) > 0) {
        wh2 = which(sampleVars %in% type2conflict)
        new2 = paste0("sample_", sampleVars[wh2])
        warning("The sample variables: \n", paste0(sampleVars[wh2], 
                                                   collapse = ", "), "\n have been renamed to: \n", 
                paste0(new2, collapse = ", "), "\n", "to avoid conflicts with taxonomic rank names.")
        colnames(sample_data(physeq))[wh2] <- new2
    }
    otutab = otu_table(physeq)
    if (!taxa_are_rows(otutab)) {
        otutab <- t(otutab)
    }
    mdf = melt(as(otutab, "matrix"), value.name = "Abundance")
    colnames(mdf)[1] <- "OTU"
    colnames(mdf)[2] <- "Sample"
    colnames(mdf)[3] <- "Abundance"
    mdf$OTU <- as.character(mdf$OTU)
    mdf$Sample <- as.character(mdf$Sample)

    if (!is.null(sampleVars)) {
        sdf = data.frame(sample_data(physeq), stringsAsFactors = FALSE)
        sdf$Sample <- sample_names(physeq)
        mdf <- merge(mdf, sdf, by.x = "Sample")
    }
    if (!is.null(rankNames)) {
        TT = access(physeq, "tax_table")
        TT <- TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
        tdf = data.frame(TT, OTU = taxa_names(physeq), stringsAsFactors = FALSE)
        mdf <- merge(mdf, tdf, by.x = "OTU")
    }
    mdf = mdf[order(mdf$Abundance, decreasing = TRUE), ]
    return(mdf)
}
