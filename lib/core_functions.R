#' select_core
#'
#' get the core microbiome from a phyloseq object 
#' 
#'@description Only the OTU's with hits in all samples 
#'              in the otu_table, will be selected. The 
#'              other components of the phyloseq object 
#'              will then be reduced and a new phyloseq 
#'              object created.
#' 
#'@param phyloseq   phyloseq object
#'
#'@return phyloseq object
#'@export
#'
select_core <- function(phyloseq) {
    # split up the phyloseq object in the 3 seperate parts
    otu_tbl <- otu_table(phyloseq)
    sample_tbl <- sample_data(phyloseq)
    tax_tbl <- tax_table(phyloseq)
    
    # find all OTUs contained in every sample
    no_zero = apply(otu_tbl, 1, function(row) all(row !=0 ))
    
    # get desired subset of the OTU's 
    otu_tbl <- otu_tbl[no_zero,]
    
    # extract tax_ids of the remaining taxa
    remain_taxa <- rownames(otu_tbl)
    
    # subset the desired taxa from tax_table to make both equal
    tax_tbl <- tax_tbl[rownames(tax_tbl) %in%  remain_taxa]
    
    # build a new phyloseq object
    phyloseq <- phyloseq(otu_tbl, sample_tbl, tax_tbl)
    
    return(phyloseq)
}

#' coreOtus
#'
#' create a specialized structure for rarefaction 
#' 
#'@description create a specialized structure for the 
#'             generation of rarefaction curves for the 
#'             core otus - this part seperate the core otus 
#'             from the phyloseq object         
#'
#'@param x              phyloseq object
#'@param n              most abundant OTU's (NULL == all OTU's)
#'@param min_fraction   specify how often an OTU must be 
#'                      observed across the samples 
#'                      ( 1 == all samples)
#'                      
#'@keyword internal
#'@return internal structure
#'
coreOtus <- function(x, 
                     n = NULL, 
                     min_fraction = 1,
                     rank) {
    
    n_samples <- nsamples(x)
    # if min_fraction = 1, OTUs must be observed across 
    # all samples to count as core
    min_samples_obs <- ceiling(n_samples * min_fraction)
    # purge all taxa that don't occur in 
    # min_samples_obs (core OTUSs)
    x.core <- prune_taxa(apply(otu_table(x),
                               1, 
                               function(x) sum(x > 0) >= min_samples_obs), x)
    
    # The rest is just looking for a common name for these 
    # core OTUs
    # sort OTUs by decreasing abundance
    otus_by_abundance <- sort(taxa_sums(x.core), decreasing = TRUE)
    tx <- tax_table(x.core)[names(otus_by_abundance), ]
    
    # optionally select only the n most abundant core OTUs
    if (!is.null(n)) {
        tx <- tx[1:n, ]
    }
    tx <- as.matrix(tx@.Data)
    tx[tx[, "species"] == "unclassified", "species"] <- "sp."
    tx[tx == "unclassified"] <- ""
    otu <- gsmisc::trim(paste0(tx[, "genus"], " ", 
                               tx[, "species"]), 
                        trim = "^ sp\\.$")
    # Count down from family level
    while (any(!nzchar(otu)) && rank > 0) {
        otu[!nzchar(otu)] <- tx[!nzchar(otu), rank]
        rank <- rank - 1
    }
    # adjust row.names and return structure
    names(otu) <- rownames(tx)
    return(structure(otu, abundance = otus_by_abundance[rownames(tx)]))
}

#' rarefy.coreOtus 
#'
#' calculate rarefaction for OTU's of the core microbiome
#'
#'@description create a specialized structure for the generation 
#'             of rarefaction curves for the core otus - this part 
#'             generate a number of rarefaction objects from the 
#'             same core otus to make the rarefaction curves smoother
#'
#'@param x              phyloseq object
#'@param n              most abundant OTU's (NULL == all OTU's)
#'@param raremax        max number of rarefaction steps
#'@param nsteps         steps for rarefaction to make curves smoother
#'                      
#'@keyword internal
#'@return internal structure
#'
rarefy.coreOtus <- function(phyloseq, 
                            n = 1, 
                            raremax = NULL, 
                            nsteps = 1) {
    # check if given object is phyloseq object
    stopifnot(is(phyloseq, "phyloseq"))
    
    # extract and transpose the OTU table from the phyloseq object 
    phyloseq.mat <- t(phyloseq@otu_table@.Data)
    
    # default to rarefying to smallest sample 
    # sizemax(colSums(otu_table(pgg_amate)))
    if (is.null(raremax)) {
        raremax <- min(rowSums(phyloseq.mat))
    }
    
    # get stepsize
    stepsize <- trunc(raremax/nsteps)
    # create loop index for every sample
    sample <- seq(from = stepsize, to = raremax, by = stepsize)
    xr <- vector("list", length(sample))
    # calculate rarefaction for every sample
    for (i in seq_along(sample)) {
        xr[[i]] <- replicate(n, {
            cat(".") # progress
            otu_table(phyloseq) <- otu_table(t(Rrarefy(phyloseq.mat, 
                                                       sample[i])), 
                                             taxa_are_rows = TRUE)
            coreOtus(phyloseq)
        }, 
        simplify = FALSE)
        cat("\n") # progress
    }
    # create and return specialised structure
    xr <- structure(xr, nreads = sample)
    return(xr)
}

core_otus <- function(x, n = NULL, min_fraction = 1, rank) {
    n_samples <- nsamples(x)
    min_samples_obs <- ceiling(n_samples*min_fraction)
    x <- prune_taxa(apply(otu_table(x), 1, function(x) sum(x > 0) >= min_samples_obs), x)
    ## sort by abundance
    otus_by_abundance <- sort(taxa_sums(x), TRUE)
    tx <- tax_table(x)[names(otus_by_abundance), ]
    if (!is.null(n)) {
        tx <- tx[1:n, ]
    }
    tx <- as.data.frame(tx, stringsAsFactors = FALSE)
    tx[tx$species == "unclassified", "species"] <- "sp."
    tx[tx == "unclassified"] <- ""
    otu <- rmisc::trim(paste0(tx[, "genus"], " ", tx[, "species"]), trim = "^ sp\\.$")
    while (any(are_empty(otu)) && rank > 0) {
        otu[are_empty(otu)] <- tx[are_empty(otu), rank]
        rank <- rank - 1
    }
    names(otu) <- rownames(tx)
    structure(otu, abundance = otus_by_abundance[rownames(tx)])
}

compactZero <- function(x) {
    is.zero = function(x) x == 0
    x[!vapply(x, is.zero, FALSE)]
}

most_abundant_taxa <- function(x, n) {
    names(sort(taxa_sums(x), decreasing=TRUE)[1:n])
}
