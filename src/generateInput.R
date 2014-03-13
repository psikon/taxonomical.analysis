generateInput <- function(data = list(), name) {
    if (is.empty(data)) {
        message("Error: list of samples is empty!")
        phylo <- data
    } else {
        # get all metadata objects
        source("src/metadata.R")
        # connect to all databses
        source("src/connections.R")
        # generate biom object from databases and write it to disk
        biom <- to_biom(data)
        write_biom(biom, biom_file = paste0("data/",name,".biom"))
        # create phyloseq object
        phylo <- import_biom(paste0("data/",name,".biom"), parallel = TRUE)
        # adjuxt colnames of the tax_table to map metaR specifications
        colnames(tax_table(phylo)) = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    }
    return(phylo)
}