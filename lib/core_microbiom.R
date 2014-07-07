# get the core microbiome from a phyloseq object by only select the OTUs, 
# that have hits in all samples in the otu_table, the other components of 
# the phyloseq object will then be reduced and a new phyloseq object created
get.coreMicrobiome <- function(phyloseq) {
    # split up the phyloseq object in the 3 seperate parts
    otu_tbl <- otu_table(phyloseq)
    sample_tbl <- sample_data(phyloseq)
    tax_tbl <- tax_table(phyloseq)
    # find all OTUs contained in every sample
    no_zero = apply(otu_tbl, 1, function(row) all(row !=0 ))
    # get desired subset the OTUs 
    otu_tbl <- otu_tbl[no_zero,]
    # extract tax_ids of the remaining taxa
    remain_taxa <- rownames(otu_tbl)
    # subset the desired taxa from tax_table to make both equal
    tax_tbl <- tax_tbl[rownames(tax_tbl) %in%  remain_taxa]
    # build a new phyloseq object
    phyloseq <- phyloseq(otu_tbl, sample_tbl, tax_tbl)
    return(phyloseq)
}

# create a venn diagramm with the two core microbiome of the 
# free living and the mariculture HoldingCondition. 
# The plot can only be saved in a .tiff file!
plot.coreMicrobiome <- function(phyloseq, file) {
    free <- get.coreMicrobiome(get.free(phyloseq))
    aqua <- get.coreMicrobiome(get.aqua(phyloseq))
    p <- venn.diagram(x = list(
        free = rownames(otu_table(free)),
        aqua = rownames(otu_table(aqua))),
        filename = file, main = "Core Microbiome",
        main.cex = 2, main.fontface = "bold",
        col = c("red","blue"), fill = c("red","blue"),
        lwd = 2, alpha = 0.3, cex = 1.0, 
        cat.cex = 1.0, cat.fontface = "bold")
    return(p)
}
# gather informations about important features of the core microbiome and save 
# them in a data.frame. This will be:
# tax id          - taxonomical identification number of the NCBI taxonomy database
# scientific.name - name of the OTU in the NCBI taxonomy database
# tax.rank        - taxonomic rank in the NCBI taxonomy database
# count           - number of occurrance in every sample of the core microbiome
get.coreTable <- function(phyloseq, file = NULL, col.names = F, row.names = F) {
    # build a data.frame with desired inormations
    data <- as.data.frame(t(sapply(seq_along(1:nrow(tax_table(phyloseq))), function(x){
        c("tax id" = rownames(otu_table(phyloseq)[x]),
          "scientific.name" = sub("_"," ",last.taxa(tax_table(phyloseq)[x])),
          "tax.rank" = last.rank(tax_table(phyloseq)[x]),
          "count"= as.vector(otu_table(phyloseq)[x]))
    })))
    # adjust column names
    colnames(data) <- c("tax.id", "scientific.name", "tax.rank", 
                        paste("counts", colnames(otu_table(phyloseq))))
    # save the data frame in a tab separeted file
    if(!is.null(file)) {
        write.table(data, file, sep = "\t", 
                    quote = F, row.names = row.names,
                    col.names = col.names)
    }
    return(data) 
}

# create a specialized structure for the generation of rarefaction curves for the 
# core otus - this part seperate the core otus from the phyloseq object
coreOtus <- function(x, n = NULL, min_fraction = 1) {
    n_samples <- nsamples(x)
    ## if min_fraction = 1, OTUs must be observed across all samples to count as core
    min_samples_obs <- ceiling(n_samples * min_fraction)
    ## purge all taxa that don't occur in min_samples_obs (core OTUSs)
    x.core <- prune_taxa(apply(otu_table(x), 1, function(x) sum(x > 0) >= min_samples_obs), x)
    
    ## The rest is just looking for a common name for these core OTUs
    ## sort OTUs by decreasing abundance
    otus_by_abundance <- sort(taxa_sums(x.core), decreasing = TRUE)
    tx <- tax_table(x.core)[names(otus_by_abundance), ]
    ## optionally select only the n most abundant core OTUs
    if (!is.null(n)) {
        tx <- tx[1:n, ]
    }
    tx <- as.matrix(tx@.Data)
    tx[tx[, "species"] == "unclassified", "species"] <- "sp."
    tx[tx == "unclassified"] <- ""
    otu <- gsmisc::trim(paste0(tx[, "genus"], " ", tx[, "species"]), trim = "^ sp\\.$")
    rank <- 5 ## Count down from family level
    while (any(!nzchar(otu)) && rank > 0) {
        otu[!nzchar(otu)] <- tx[!nzchar(otu), rank]
        rank <- rank - 1
    }
    names(otu) <- rownames(tx)
    return(structure(otu, abundance = otus_by_abundance[rownames(tx)]))
}

# create a specialized structure for the generation of rarefaction curves for the 
# core otus - this part generate a number of rarefaction objects from the same core otus 
# to make the rarefaction curves smoother
rarefy.coreOtus <- function(phyloseq, n = 1, raremax = NULL, nsteps = 1) {
    stopifnot(is(phyloseq, "phyloseq"))
    ## extract and transpose the OTU table from the phyloseq object 
    phyloseq.mat <- t(phyloseq@otu_table@.Data)
    ## default to rarefying to smallest sample sizemax(colSums(otu_table(pgg_amate)))
    if (is.null(raremax)) {
        raremax <- min(rowSums(phyloseq.mat))
    }
    stepsize <- trunc(raremax/nsteps)
    sample <- seq(from = stepsize, to = raremax, by = stepsize)
    xr <- vector("list", length(sample))
    for (i in seq_along(sample)) {
        xr[[i]] <- replicate(n, {
            cat(".") # progress
            otu_table(phyloseq) <- otu_table(t(Rrarefy(phyloseq.mat, sample[i])), taxa_are_rows = TRUE)
            core_otus(phyloseq)
        }, simplify = FALSE)
        cat("\n") # progress
    }
    xr <- structure(xr, nreads = sample)
    return(xr)
}

# create a plot, showing the rarefaction curves for the core otus
plot.rarifyCoreOtus <- function(phyloseq, n = 25, steps = 20, 
                                file = NULL, title = "Core OTUs - Rarefaction curves") {
    # create 3 seperate core otus objects with a n number of replicates 
    rare.core.all <- rarefy.coreOtus(phyloseq, n = n, nsteps = steps)
    rare.core.free <- rarefy.coreOtus(get.free(phyloseq), n = n, nsteps = steps)
    rare.core.aqua <- rarefy.coreOtus(get.aqua(phyloseq), n = n, nsteps = steps)
    
    #### rarefaction curves ####
    make.coreOTU.df <- function(x = rare.core.all, group) {
        mean.core <- sapply(x, function(x) mean(sapply(x, length)))
        sd.core <- sapply(x, function(x) sd(sapply(x, length)))
        data.frame(
            stringsAsFactors = FALSE,
            group = group,
            core = mean.core,
            upper = mean.core + 1.96 * sd.core,
            lower = mean.core - 1.96 * sd.core,
            reads = attr(x, "nreads")
        )
    }
    
    # create a data.frame from the core otus
    df.all <- make.coreOTU.df(rare.core.all, "All species")
    df.free <- make.coreOTU.df(rare.core.free, "free living")
    df.aqua <- make.coreOTU.df(rare.core.aqua, "mariculture")
    # combine the three data.frames to one data.frame usable by ggplot2
    core.df <- do.call("rbind", list(df.all, df.free, df.aqua))
    # insert a variable for seperation in legend and color
    core.df$group <- factor(core.df$group, labels = c("E. fuscoguttatus - all", 
                                                      "E. fuscoguttatus - free living",
                                                      "E. fuscoguttatus - mariculture"))
    # create the ggplot
    p <- ggplot(core.df, aes(x = reads, y = core, colour = group)) +
                geom_line(size = 1) +
                xlab("Number of reads") +
                ylab("Number of shared OTUs") +
                theme(legend.position = c(0.65, 0.15)) + 
                scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
                scale_y_continuous(breaks = pretty_breaks(n = 5)) +
                ggtitle(title)
    # save plot in file
    if(!is.null(file)) ggsave(file)
    return(p)

}
