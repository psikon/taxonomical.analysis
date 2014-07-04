# wrapper function to create a biom file from a list of taxonomyReportDB objects
generate.biomFile <- function(name, data = list(sample60.old, sample62.old, sample66.old,
                                                sample68.old, sample70.old, sample76.old, 
                                                sample78.old, sample80.old, sample82.old)){
    # generate biom object from databases and write it to disk
    biom <- to_biom(data)
    write_biom(biom, biom_file = paste0("data/", name, ".biom"))
    return(paste0("data/", name, ".biom"))
}

# wrapper to import a biom file to a phyloseq object with some modifications
generate.phyloseq <- function(name) {
    # create phyloseq object
    phyloseq <- import_biom(paste0("data/", name, ".biom"), parallel = TRUE)
    # adjust colnames of the tax_table to map metaR specifications
    colnames(tax_table(phyloseq)) = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    phyloseq
}