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