#' plot.rareCurve
#'
#' calculate and plor rarefaction curves for phyloseq object
#'
#'@description create rarefaction curves for a phyloseq object, 
#'             where the number of OTUs against the number of 
#'             used reads is used to show a hopefully a 
#'             saturation over all samples 
#' 
#'@param phyloseq       phyloseq object
#'@param stepsize       steps for rarefaction
#'@param file           output file
#'
#'@return ggplot2 object
#'@export
#'
plot.rareCurve <- function(phyloseq, 
                           stepsize, 
                           file = NULL) {
    
    # parse important variables from phyloseq object
    
    # transposed abundance table (taxa are cols) 
    abundance <- t(as.data.frame(otu_table(phyloseq)))
    # total number of hits for every samples
    hits <- sample_sums(phyloseq)
    # total number of taxa for every sample
    species <- specnumber(abundance)
    # number of samples
    num_samples <- nrow(sample_data(phyloseq))
    
    # calculate data for rarefaction curve
    out <- lapply(seq_len(num_samples), function(i) {
        n <- seq(1, hits[i], by = stepsize)
        if(n[length(n)] != hits[i])
            n <- c(n, hits[i])
        drop(rarefy(abundance[i, ], n))
    })
    
    # convert the object for ggplot2 purposes
    
    # 1. convert all length of lists in the object to 
    # same length by adding na values
    
    # determine length of every list in the object
    lengths <- sapply(out, length)
    # create lists with na values for every list
    na_list <- Map(function(x) rep(NA, x), x = max(lengths) - lengths)
    # map these na_list to out object
    out2 <- Map(function(a, b) c(a, b), a = out, 
                b = na_list)
    # and create a matrix of the list object
    res <- do.call("cbind", out2)
    
    # adjust colnames for legend 
    #colnames(res) <- paste0("S", 1:length(out))
    colnames(res) <- sample_data(phyloseq)$SampleName
    
    # convert res matrix to ggplot2 data.frame structure
    res2 <- na.omit(melt(res))
    
    # convert first column from factor to numeric
    res2$X1 <- as.numeric(substring(as.character(res2$X1), 
                                    first = 2))
    
    # adjust colnames for easier access
    colnames(res2) <- c("seqs", "SampleName", "value")
    
    # determine minimal values for rarefy to even depth starting
    # by get list of all samples
    samples <- as.character(unique(res2$SampleName))
    
    # find last value for every sample in seq column
    v <- sapply(samples, function(x) {
        max(compactNA(res2[which(res2$SampleName == x), ]$seqs))
    })
    
    # find last value for every sample in value column 
    # (representing taxa)
    h <- sapply(samples, function(x) {
        max(compactNA(res2[which(res2$SampleName == x), ]$value))
    })
    
    # get number of taxa for sample with minimal sequences
    h <- as.numeric(h[which.min(v)])
    
    # get sample with minimal number of sequences
    v <- as.numeric(v[which.min(v)])
    res2 <- cbind(res2, h, v)
    
    # create the plot 
    p <- ggplot(res2, aes(x = seqs, 
                          y = value, 
                          colour = SampleName)) +
        geom_line(size = 1) +
        scale_x_continuous("\nNumber of Sequences", 
                           breaks = pretty_breaks(n = 10)) +
        scale_y_continuous("\nObserved species", 
                           breaks = pretty_breaks(n = 10)) +
        ggtitle("Rarefaction Curves") + theme_bw() +
        geom_hline(aes(yintercept = h), color = "grey") +
        geom_vline(aes(xintercept = v), color = "grey")
    # save plot to file
    if(!is.null(file)) ggsave(file)
    return(p)
}
#' rarify.phyloseq
#'
#' wrapper for the phyloseq::rarefy_even_depth function 
#' with needed parameters set
#' 
#'@param phyloseq   phyloseq object
#'@param rngseed    random seed set
#'@param replace    sample with or without replacement
#'@param trimOTUs   trim OTUs that are no longer observed 
#'
#'@return rarified phyloseq object
#'@keywords internal
#'
rarify.phyloseq <- function(phyloseq, 
                            rngseed = 1234, 
                            replace = F, 
                            trimOTUs = T) {
    
    return(rarefy_even_depth(phyloseq, 
                             min(sample_sums(phyloseq)),
                             rngseed = rngseed, 
                             replace = replace, 
                             trimOTUs = trimOTUs, 
                             verbose = F))
}
#' Rrarefy
#'
#' overwritten Rrarefy function from vegan, 
#' used for rarefaction of the core microbiome
#'
#'@param x       data
#'@param sample  name of sample
#'
#'@return data
#'@keywords internal
#' 
Rrarefy <- function(x, sample) {
    replace <- if (sample > min(rowSums(x))) TRUE else FALSE
    sample <- rep(sample, length = nrow(x))
    colnames(x) <- colnames(x, do.NULL = FALSE)
    nm <- colnames(x)
    for (i in seq_len(nrow(x))) {
        row <- sample(rep(nm, times = x[i, ]), size = sample[i], replace = replace)
        row <- table(row)
        ind <- names(row)
        x[i, ] <- 0
        x[i, ind] <- row
    }
    return(x)
}
