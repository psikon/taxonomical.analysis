# from phyloseq extension tutorial (http://joey711.github.io/phyloseq-extensions/DESeq.html)
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

phyloseq_to_edgeR = function(physeq, group, method = "RLE", ...) {
    require("edgeR")
    require("phyloseq")
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

plot_variance <- function(phyloseq, x.text = "log10(variance)", title = "Variance of OTUs" ) {
  hist(log10(apply(otu_table(phyloseq),1, var)), xlab = x.text, main = title)
}

deseq_differential_OTUs <- function(phyloseq, variance_threshold, alpha) {
    keepOTUs <- apply(otu_table(phyloseq), 1, var) > variance_threshold
    data <- prune_taxa(keepOTUs, phyloseq)
    #create cds
    cds <- suppressWarnings(phyloseq_to_DESeq(data, 
                                              "HoldingCondition", 
                                              sharingMode = "gene-est-only"))
    res <- nbinomTest(cds, "free living", "mariculture")
    sigtab <- res[which(res$padj < alpha),]
    sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(data)[sigtab$id, ], "matrix"))
    sigtab <- sigtab[order(sigtab$padj),]
    sigtab
}

plot_differences <- function(sigtab) {
    theme_set(theme_bw())
    scale_fill_discrete <- function(palname = "Set1", ...) {
        scale_fill_brewer(palette = palname, ...)
    }
    
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
    
    ggplot(sigtab, 
           aes(x = genus, y = log2FoldChange, color = phylum)) + 
        geom_point(size = 3) + 
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
}

deseq2_differntial_OTUs <- function(phyloseq, alpha) {
    deseq2 <- suppressWarnings(phyloseq_to_deseq2(phyloseq, ~HoldingCondition))
    deseq2 <- DESeq(deseq2, test = "Wald", fitType="local")
    res <- results(deseq2)
    sigtab <- res[which(res$padj < alpha),]
    sigtab = cbind(as(sigtab,"data.frame"), 
                   as(tax_table(phyloseq)[rownames(sigtab), ], "matrix")) 
    sigtab
}
# deseq2 <- deseq2_differntial_OTUs(bakteria, 0.3)
# plot_differences(deseq2)
# dge <- phyloseq_to_edgeR(bakteria, "HoldingCondition")
# dge
# et = exactTest(dge)
# et
# tt = topTags(et, n = nrow(dge$table),adjust.method= "BH", sort.by = "PValue")
# res = tt@.Data[[1]]
# res
# alpha = 0.003
# sigtab = res[(res$FDR < alpha),]
# sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(bakteria)[rownames(sigtab), 
#                                                                ], "matrix"))
# dim(sigtab)
# head(sigtab)
# 
# theme_set(theme_bw())
# scale_fill_discrete <- function(palname = "Set1", ...) {
#     scale_fill_brewer(palette = palname, ...)
# }
# sigtabgen = subset(sigtab, !is.na(genus))
# # Phylum order
# x = tapply(sigtabgen$logFC, sigtabgen$phylum, function(x) max(x))
# x = sort(x, TRUE)
# sigtabgen$phylum = factor(as.character(sigtabgen$phylum), levels = names(x))
# # Genus order
# x = tapply(sigtabgen$logFC, sigtabgen$genus, function(x) max(x))
# x = sort(x, TRUE)
# sigtabgen$genus = factor(as.character(sigtabgen$genus), levels = names(x))
# ggplot(sigtabgen, aes(x = genus, y = logFC, color = phylum)) + geom_point(size = 3) + 
#     theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))

