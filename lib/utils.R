rm.taxa <- function(phyloseq, taxa) {
    phyloseq <- prune_species(!grepl(paste0("\\<", taxa, "\\>"),
                                     labels(otu_table(phylo))[[1]]), phyloseq)
    return(phyloseq)
}

rm.underscore <- function(phyloseq) {
    data <- tax_table(phyloseq)
    if(any(grepl("__", data))) {
        data <- substr(data, 4, length(data))
        data[which(data == "")] <- "not classified"
    }
    tax_table(phyloseq) <- data
    return(phyloseq)
}


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

get.free <- function(phyloseq) {
    return(subset_samples(phyloseq, HoldingCondition == "free living"))
}

get.aqua <- function(phyloseq) {
    return(subset_samples(phyloseq, HoldingCondition == "mariculture"))
}

last.rank <- function(x) {
    x <- tail(colnames(x[, sapply(strsplit(x, "*__"), function(y){
        length(y) > 1
    })]), n = 1)
    return(x)
}

last.taxa <- function(x) {
    rank <- tail(colnames(x[, sapply(strsplit(x, "*__"), function(y){
        length(y)>1
    })]), n = 1)
    x <- unlist(str_split(x[,rank],"__"))[2]
    return(x)
}

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

countArtifacts <- function(phyloseq, sample) {
    with <- subset_taxa(phyloseq, superkingdom == 'k__Bacteria')
    without <- rm.taxa(with, "2")
    countReads(with, as.character(sample)) - countReads(without, as.character(sample))
}

countReads <- function(phyloseq, sample) {
    value <- sum(count(get_taxa(phyloseq, as.character(sample)))$x * count(get_taxa(phyloseq, as.character(sample)))$freq)
    return(value)
}

getOTUs <- function(phyloseq, sample) {
    return(length(which(otu_table(phyloseq)[, as.character(sample)] > 0) == TRUE))
}

