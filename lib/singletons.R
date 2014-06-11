get_singletons <- function(phyloseq, num_samples = 1, remove_samples = TRUE) {
    # get all rows not belonging to the core microbiom
    non_zero <- rownames(otu_table(get_core_microbiom(phyloseq)))
    # reduce otu_table
    otu_tbl <- otu_table(phyloseq)[-which(rownames(otu_table(phyloseq)) %in% non_zero)]
    # reduce tax_table
    tax_tbl <- tax_table(phyloseq)[-which(rownames(tax_table(phyloseq)) %in% non_zero)]
    # get sample_data
    sample_tbl <- sample_data(phyloseq)
    # define max number of zeros
    max <- nrow(sample_data(phyloseq)) - num_samples
    # find the singletons
    data <- apply(otu_tbl, 1, function(row) {
        if(table(row)["0"] >= max) {
            row
        }
    })
    # format back to matrix and remove NULL entries
    data <- do.call(rbind, data[!sapply(data, is.null)])
    # get the names of the candidates
    names <- rownames(data)
    # reduce tax_table to candidates
    tax_tbl <- tax_tbl[rownames(tax_tbl) %in%  names]
    # reduce the phyloseq object containing only samples with present singletons
    if(remove_samples) {
        # get a list of samples without a hit
        samples <- which(unlist(lapply(colnames(data), 
                                       function(column) all(data[, column] == 0 ))
                                ) == TRUE)
        # remove the non present samples from otu_table
        if (length(samples) != 0) data <- data[ , -samples]
        # remove the non present samples from sample data
        if (length(samples) != 0) sample_tbl <- sample_tbl[-samples, ]
    }
    # create a new phyloseq object from the results 
    result <- phyloseq(otu_table(data, taxa_are_rows = T),
                       sample_tbl,
                       tax_tbl)
    result
}

get_singleton_list <- function(phyloseq, file = NULL, col.names = F, row.names = F) {
    # build a data frame with desired inormations
    data <- as.data.frame(t(sapply(seq_along(1:nrow(tax_table(phyloseq))), function(x) {
        c("tax id" = rownames(otu_table(phyloseq)[x]),
          "scientific.name" = sub("_"," ",last_taxa(tax_table(phyloseq)[x])),
          "tax.rank" = last_rank(tax_table(phyloseq)[x]),
          "samples" = paste(colnames(otu_table(phyloseq)[x][, otu_table(phyloseq)[x] > 0]), 
                            sep = ",", collapse = ","),
          "counts"= paste(as.vector(otu_table(phyloseq)[x][, otu_table(phyloseq)[x] > 0]), 
                          sep = ",", collapse = ","))
    })))
    # adjust column names
    colnames(data) <- c("tax.id", "scientific.name", "tax.rank", 
                        "samples", "counts")
    # save the data frame in a tab separeted file
    if(!is.null(file)) {
        write.table(data,file, sep="\t",quote = F, row.names=row.names,col.names=col.names)
    }
    data 
}
