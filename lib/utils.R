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

#' create.QueryIDFile
#'
#' extract from taxonomyReportDB object the query_ids with 
#' successfull annotation and create file of it
#'
#'@param path       path for output file
#'@param position   position of samples in db structure
#'                      
#'@keyword internal
#'@return list of files
#'
create.QueryIDFile <- function(path = "data/functional", 
                               position = NULL) {
    
    # get all taxonomyReportDB's
    db <- get.DBConnection.new(get.metadata.list())
    
    # reduce to desired dbs
    if (!is.null(position)) {
        db <- db[position]
    }
    
    files <- lapply(db, function(x) {
        # create filename for index file
        file <- paste0(path, 
                       "/", 
                       as.character(x$.metadata["SampleName"]), ".txt")
        # extract the query_ids
        query_ids <- db_query(conn(x),
                              "Select query_def from query", 1)
        # write query_ids to file line by line
        write.table(query_ids, 
                    file = file, 
                    quote = FALSE, 
                    row.names = F, 
                    col.names = F)
        # return filepath
        file
    })
    return(files)
}

#' create.FastaFromTaxonomyReportDB
#'
#' extract from fasta input file the reads with annotation 
#' result in database 
#'
#'@param path           path for output file
#'@param fasta.files    list of input fasta files
#'@param position       position of samples in db structure
#'                      
#'@keyword internal
#'
create.FastaFromTaxonomyReportDB <- function(path = "data/functional", 
                                             fasta.files, 
                                             position = NULL) {
    # create the index files
    index.files <- create.QueryIDFile(path = path, 
                                      position = position)
    for (i in seq_along(index.files)) {
        # resolve filename
        name <- strsplit(strsplit(index.files[[i]], 
                                  split = "/")[[1]][3], '.', 1)[[1]][1]
        message(paste0("create fasta file for: ", name))
        # construct output path
        output <- paste0(path, "/", name, ".fasta")
        # generate commandline string
        cmd <- paste0("python src/fastaExtractor.py --header ", 
                      index.files[[i]], " --fasta ", 
                      fasta.files[i], " --output ",
                      output)
        # run src/fastaExtractor.py
        system(cmd, wait = T)
    }
}

#' phyloseq_to_DESeq
#'
#' import a phyloseq object to DESeq 
#' 
#'@description from phyloseq extension tutorial 
#'              (http://joey711.github.io/phyloseq-extensions/DESeq.html)
#' 
#'@param physeq         phyloseq object
#'@param designFactor   test variable in metadata
#'@param fitType        see documentation
#'@param locfunc        see documentation
#'
#'@return DESeq object
#'@export
#'
phyloseq_to_DESeq = function(physeq, 
                             designFactor, 
                             fitType = "local", 
                             locfunc = median, 
                             ...) {
    
    # Enforce Orientation
    if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
    }
    
    # Convert to matrix, round up to nearest integer
    x = ceiling(as(otu_table(physeq), "matrix")) + 1L
    
    # Add taxonomy data, if present.
    if (!is.null(tax_table(physeq, FALSE))) {
        taxADF = as(data.frame(as(tax_table(physeq), 
                                  "matrix")), 
                    "AnnotatedDataFrame")
    } else {
        taxADF = NULL
    }
    
    # Add sample data if present
    if (!is.null(sample_data(physeq, FALSE))) {
        samplesADF = as(data.frame(sample_data(physeq)), 
                        "AnnotatedDataFrame")
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
    
    cds = newCountDataSet(x, 
                          conditions = designFactor, 
                          phenoData = samplesADF, 
                          featureData = taxADF)
    # First, estimate size factors, then estimate dispersion.
    cds = estimateSizeFactors(cds, 
                              locfunc = locfunc)
    # Now dispersions/variance estimation, 
    # passing along additional options
    cds = estimateDispersions(cds, 
                              fitType = fitType, ...)
    
    return(cds)
}
#' phyloseq_to_edgeR
#'
#' import a phyloseq object to edgeR 
#' 
#'@description from phyloseq extension tutorial 
#'              (http://joey711.github.io/phyloseq-extensions/edgeR.html)
#' 
#'@param physeq         phyloseq object
#'@param group          test variable in metadata
#'@param method         see documentation
#'@param ...            see documentation
#'
#'@return DESeq object
#'@export
#'
phyloseq_to_edgeR = function(physeq, 
                             group, 
                             method = "RLE", 
                             ...) {
    # Enforce orientation.
    if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
    }
    x = as(otu_table(physeq), "matrix")
    
    # Add one to protect against overflow, log(0) issues.
    x = x + 1
    
    # Check `group` argument
    if (identical(all.equal(length(group), 1),
                  TRUE) & nsamples(physeq) > 1) {
        # Assume that group was a sample variable name 
        # (must be categorical)
        group = get_variable(physeq, group)
    }
    
    # Define gene annotations (`genes`) as tax_table
    taxonomy = tax_table(physeq, errorIfNULL = FALSE)
    if (!is.null(taxonomy)) {
        taxonomy = data.frame(as(taxonomy, "matrix"))
    }
    
    # Now turn into a DGEList
    y = DGEList(counts = x, 
                group = group, 
                genes = taxonomy, 
                remove.zeros = TRUE, 
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


