# remove a taxa with a specific tax id from phyloseq object
rm.taxa <- function(phyloseq, taxa) {
    phyloseq <- prune_species(!grepl(paste0("\\<", taxa, "\\>"),
                                     labels(otu_table(phyloseq))[[1]]), phyloseq)
    return(phyloseq)
}

# remove the undercore syntax from tax_table originating from the metaR or RDP process 
# for beatiful plotting
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

# count the taxonomical levels of of a OTU
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

# convert a absolute number to percent value
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

# get only the free living samples
get.free <- function(phyloseq) {
    return(subset_samples(phyloseq, HoldingCondition == "free living"))
}

# get only the mariculture samples
get.aqua <- function(phyloseq) {
    return(subset_samples(phyloseq, HoldingCondition == "mariculture"))
}

# get the last rank of a linage from tax_table
last.rank <- function(x) {
    x <- tail(colnames(x[, sapply(strsplit(x, "*__"), function(y){
        length(y) > 1
    })]), n = 1)
    return(x)
}

# get only the last taxa from a tax_table
last.taxa <- function(x) {
    rank <- tail(colnames(x[, sapply(strsplit(x, "*__"), function(y){
        length(y)>1
    })]), n = 1)
    x <- unlist(str_split(x[,rank],"__"))[2]
    return(x)
}

# user defined theme for richness plot
get.richnessTheme <- function() {
    ggtheme.alpha <- theme_bw() + theme(
                                        legend.position = "top",
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

# count the number of Artifacts in phyloseq object. An artifact is defined by mapping 
# only to the superkingdom level in the linage
countArtifacts <- function(phyloseq, sample) {
    # extract all bacterial taxa
    with <- subset_taxa(phyloseq, superkingdom == 'k__Bacteria')
    # remove taxa with maximal mapping to superkingdom
    without <- rm.taxa(with, "2")
    # calculate number of atrifacts and return them
    return(countReads(with, as.character(sample)) - countReads(without, as.character(sample)))
}

# count the total number of taxonomical hits of a sample
countReads <- function(phyloseq, sample) {
    value <- sum(count(get_taxa(phyloseq, as.character(sample)))$x * count(get_taxa(phyloseq, as.character(sample)))$freq)
    return(value)
}
# get the number of OTUs of a sample
getOTUs <- function(phyloseq, sample) {
    return(length(which(otu_table(phyloseq)[, as.character(sample)] > 0) == TRUE))
}

# from: Martin Turjak
# at: http://stackoverflow.com/questions/15343338/
# recursion function for newick string creation
traverse <- function(a,i,innerl){
    if(i < (ncol(df))){
        alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
        desc <- NULL
        if(length(alevelinner) == 1) (newickout <- traverse(alevelinner,i+1,innerl))
        else {
            for(b in alevelinner) desc <- c(desc,traverse(b,i+1,innerl))
            il <- NULL; if(innerl==TRUE) il <- a
            (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
        }
    }
    else { (newickout <- a) }
}

# from: Martin Turjak
# at: http://stackoverflow.com/questions/15343338/
# data.frame to newick function
df2newick <- function(df, innerlabel = FALSE) {
    # determine root
    alevel <- as.character(unique(df[,1]))
    newick <- NULL
    # traverse throug data.frame and create newick string
    for(x in alevel) newick <- c(newick, traverse(x, 1, innerlabel))
    # add outer strings
    (newick <- paste("(",paste(newick,collapse=","),");",sep=""))
    return(newick)
}

# overwrite to fix bug in plot_bar
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

# overwrite to fix bug in plot_bar
plot_bar <- function (physeq, x = "Sample", y = "Abundance", fill = NULL, 
                      title = NULL, facet_grid = NULL) {
    mdf = psmelt(physeq)
    p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
    p = p + geom_bar(stat = "identity", position = "stack", color = "black")
    p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
    if (!is.null(facet_grid)) {
        p <- p + facet_grid(facet_grid)
    }
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    return(p)
}

