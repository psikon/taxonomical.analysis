# from phyloseq extension tutorial (http://joey711.github.io/phyloseq-extensions/DESeq.html)
# import a phyloseq object to DESeq
phyloseq_to_DESeq = function(physeq, designFactor, fitType = "local", locfunc = median, 
                             ...) {
    # Enforce Orientation
    if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
    }
    # Convert to matrix, round up to nearest integer
    x = ceiling(as(otu_table(physeq), "matrix")) + 1L
    # Add taxonomy data, if present.
    if (!is.null(tax_table(physeq, FALSE))) {
        taxADF = as(data.frame(as(tax_table(physeq), "matrix")), "AnnotatedDataFrame")
    } else {
        taxADF = NULL
    }
    # Add sample data if present
    if (!is.null(sample_data(physeq, FALSE))) {
        samplesADF = as(data.frame(sample_data(physeq)), "AnnotatedDataFrame")
    } else {
        samplesADF = NULL
    }
    # Initalize the count data sets.
    if (identical(length(designFactor), 1L)) {
        if (designFactor %in% sample_variables(physeq)) {
            designFactor <- get_variable(physeq, designFactor)
        } else {
            stop("You did not provide an appropriate `designFactor` argument. Please see documentation.")
        }
    }
    cds = newCountDataSet(x, conditions = designFactor, phenoData = samplesADF, 
                          featureData = taxADF)
    # First, estimate size factors, then estimate dispersion.
    cds = estimateSizeFactors(cds, locfunc = locfunc)
    # Now dispersions/variance estimation, passing along additional options
    cds = estimateDispersions(cds, fitType = fitType, ...)
    return(cds)
}
# from phyloseq extension tutorial (http://joey711.github.io/phyloseq-extensions/edgeR.html)
# import a phyloseq object to edgeR
phyloseq_to_edgeR = function(physeq, group, method = "RLE", ...) {
    # Enforce orientation.
    if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
    }
    x = as(otu_table(physeq), "matrix")
    # Add one to protect against overflow, log(0) issues.
    x = x + 1
    # Check `group` argument
    if (identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1) {
        # Assume that group was a sample variable name (must be categorical)
        group = get_variable(physeq, group)
    }
    # Define gene annotations (`genes`) as tax_table
    taxonomy = tax_table(physeq, errorIfNULL = FALSE)
    if (!is.null(taxonomy)) {
        taxonomy = data.frame(as(taxonomy, "matrix"))
    }
    # Now turn into a DGEList
    y = DGEList(counts = x, group = group, genes = taxonomy, remove.zeros = TRUE, 
                ...)
    # Calculate the normalization factors
    z = calcNormFactors(y, method = method)
    # Check for division by zero inside `calcNormFactors`
    if (!all(is.finite(z$samples$norm.factors))) {
        stop("Something wrong with edgeR::calcNormFactors on this data,\n         non-finite $norm.factors, consider changing `method` argument")
    }
    # Estimate dispersions
    return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

# plot the variance structure of the OTUs to discover the right variance value for deseq
plot.variance <- function(phyloseq, x.text = "log10(variance)", title = "Variance of OTUs" ) {
  hist(log10(apply(otu_table(phyloseq), 1, var)), xlab = x.text, main = title)
}

# import a phyloseq object to deseq and calculate the differential appearing OTUs
deseq.differential.OTUs <- function(phyloseq, variance_threshold = 4, alpha = 0.2) {
    # remove underscore syntax from phyloseq tax_table
    phyloseq <- rm.underscore(phyloseq)
    # remove OTUs not mapping criteria
    keepOTUs <- apply(otu_table(phyloseq), 1, var) > variance_threshold
    data <- prune_taxa(keepOTUs, phyloseq)
    #create deseq cds object
    cds <- suppressWarnings(phyloseq_to_DESeq(data, 
                                              "HoldingCondition", 
                                              sharingMode = "gene-est-only"))
    # calculate test statistics
    res <- nbinomTest(cds, "free living", "mariculture")
    # create significance table
    sigtab <- res[which(res$padj < alpha),]
    if (nrow(sigtab) == 0) {
        warning("Significance table is empty! May be alpha is to low?")
        return(NULL)
    } else {
        sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(data)[sigtab$id, ], "matrix"))
        sigtab <- sigtab[order(sigtab$padj),]
        return(sigtab)
    }
}

# import a phyloseq object to deseq2 and calculate the differential appearing OTUs
deseq2.differential.OTUs <- function(phyloseq, alpha) {
    # convert phyloseq object to deseq2 object
    deseq2 <- suppressWarnings(phyloseq_to_deseq2(rm.underscore(phyloseq), ~HoldingCondition))
    # calculate statistical test
    deseq2 <- DESeq(deseq2, test = "Wald", fitType = "local")
    # extract results from object
    res <- results(deseq2)
    # create significance table
    sigtab <- res[which(res$padj <= alpha), ]
    if (nrow(sigtab) == 0) {
        warning("Significance table is empty! May be alpha is to low?")
        return(NULL)
    } else {
        sigtab <-  cbind(as(sigtab,"data.frame"), 
                        as(tax_table(rm.underscore(phyloseq))[rownames(sigtab), ], "matrix")) 
        return(sigtab)
    }
}

# import a phyloseq object to edgeR and calculate the differential appearing OTUs
edgeR.differential.OTUs <- function(phyloseq, alpha) {
    # convert ophyloseq object to edgeR object
    dge <- phyloseq_to_edgeR(rm.underscore(phyloseq), "HoldingCondition")
    # calculate statistics
    et <- exactTest(dge)
    # convert the object to usable structure
    tt <- topTags(et, n = nrow(dge$table), adjust.method= "BH", sort.by = "PValue")
    res <- tt@.Data[[1]]
    sigtab = res[(res$FDR < alpha),]
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(rm.underscore(phyloseq))[rownames(sigtab), 
                                                                    ], "matrix"))
    return(sigtab)   
}

# plot the log2FoldChange of the differential appearing OTUs, calculated by deseq, deseq2 or edgeR object 
plot.differential.OTUs <- function(sigtab, file = NULL, origin = "deseq") {
    theme_set(theme_bw())
    scale_fill_discrete <- function(palname = "Set1", ...) {
        scale_fill_brewer(palette = palname, ...)
    }
    if(is.null(sigtab)) {
        stop('significance table is empty')
    }
    if (origin == "deseq") {
        # Phylum order
        x <- tapply(sigtab$log2FoldChange, sigtab$phylum, function(x)
            max(x))
        x <- sort(x,TRUE)
        sigtab$phylum <- factor(as.character(sigtab$phylum), levels = names(x))
        
        # Genus order
        x <- tapply(sigtab$log2FoldChange, sigtab$genus, function(x)
            max(x))
        x <- sort(x,TRUE)
        sigtab$genus <- factor(as.character(sigtab$genus), levels = names(x))
    } else if (origin == "edgeR") {
        sigtab = subset(sigtab, !is.na(genus))
        # Phylum order
        x = tapply(sigtab$logFC, sigtab$phylum, function(x) max(x))
        x = sort(x, TRUE)
        sigtab$phylum = factor(as.character(sigtab$phylum), levels = names(x))
        # Genus order
        x = tapply(sigtab$logFC, sigtab$genus, function(x) max(x))
        x = sort(x, TRUE)
        sigtab$genus = factor(as.character(sigtab$genus), levels = names(x))
        colnames(sigtab)[8] <- "log2FoldChange" 
    }
    # create the plot
    p <- ggplot(sigtab, aes(x = genus, y = log2FoldChange, color = phylum)) + 
         geom_point(size = 3) + 
         theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
    # save plot to file
    if(!is.null(file)) ggsave(file)
    return(p)
}
