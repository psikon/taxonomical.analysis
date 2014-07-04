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


core_otus <- function(x, n = NULL, min_fraction = 1) {
    n_samples <- nsamples(x)
    ## if min_fraction = 1, OTUs must be observed across all samples to count as core
    min_samples_obs <- ceiling(n_samples*min_fraction)
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
    structure(otu, abundance = otus_by_abundance[rownames(tx)])
}

rarefy_core_otus <- function(x, n = 1, raremax = NULL, nsteps = 1) {
    stopifnot(is(x, "phyloseq"))
    ## extract and transpose the OTU table from the phyloseq object x
    x.mat <- t(x@otu_table@.Data)
    ## default to rarefying to smallest sample sizemax(colSums(otu_table(pgg_amate)))
    if (is.null(raremax)) {
        raremax <- min(rowSums(x.mat))
    }
    stepsize <- trunc(raremax/nsteps)
    sample <- seq(from = stepsize, to = raremax, by = stepsize)
    xr <- vector("list", length(sample))
    for (i in seq_along(sample)) {
        xr[[i]] <- replicate(n, {
            cat(".") # progress
            otu_table(x) <- otu_table(t(Rrarefy(x.mat, sample[i])), taxa_are_rows = TRUE)
            core_otus(x)
        }, simplify = FALSE)
        cat("\n") # progress
    }
    xr <- structure(xr, nreads = sample)
    xr
}

rarify_core_otus <- function() {
    n <- 25
    nsteps <- 20
    core_all <- core_otus(bakteria)
    rare_core_all <- rarefy_core_otus(bakteria, n = n, nsteps = nsteps)
    save(rare_core_all, file = "cache/rare_core_all.rda")
    
    free <- get_free(bakteria)
    max(colSums(otu_table(free)))
    core_free <- core_otus(free)
    rare_core_free <- rarefy_core_otus(free, n = n, nsteps = nsteps)
    
    aqua <- get_aqua(bakteria)
    max(colSums(otu_table(aqua)))
    core_aqua <- core_otus(aqua)
    rare_core_aqua <- rarefy_core_otus(aqua, n = n, nsteps = nsteps)
    
    #### rarefaction curves ####
    make_core_otu_df <- function(x = rare_core_all, group) {
        mean_core <- sapply(x, function(x) mean(sapply(x, length)))
        sd_core <- sapply(x, function(x) sd(sapply(x, length)))
        data.frame(
            stringsAsFactors = FALSE,
            group = group,
            core = mean_core,
            upper = mean_core + 1.96*sd_core,
            lower = mean_core - 1.96*sd_core,
            reads = attr(x, "nreads")
        )
    }
    
    
    df.all <- make_core_otu_df(rare_core_all, "All species")
    df.free <- make_core_otu_df(rare_core_free,"free living")
    df.aqua <- make_core_otu_df(rare_core_aqua,"mariculture")
    
    
    core.df <- do.call("rbind", list(df.all, df.free, df.aqua))
    
    core.df$group <- factor(core.df$group, labels = c("E. fuscoguttatus - all", 
                                                      "E. fuscoguttatus - free living",
                                                      "E. fuscoguttatus - mariculture"))
    ggplot(core.df, aes(x = reads, y = core, colour = group)) +
        geom_line(size = 1) +
        xlab("Number of reads") +
        ylab("Number of shared OTUs") +
        theme(legend.position = c(0.65, 0.15)) + scale_x_continuous(breaks = pretty_breaks(n=10))
    
}
