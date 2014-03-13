generateBiomFile <- function(name){
    # connect to all databses
    source("src/connections.R")
    data <- list(sample60, sample62, 
                 sample66, sample68, sample70, 
                 sample76, sample78, 
                 sample80, sample82)
    # generate biom object from databases and write it to disk
    biom <- to_biom(data)
    write_biom(biom, biom_file = paste0("data/",name,".biom"))
    # clean up metadata
    rm(metadata60, metadata62, metadata64, metadata66, metadata68, metadata70, 
       metadata74, metadata76, metadata78, metadata80, metadata82, envir= globalenv())
    # clean up connections
    rm(sample60, sample62, sample66, sample68, sample70, sample76,
       sample78, sample80, sample82, envir= globalenv())
    return(paste0("data/",name,".biom"))
}

generatePhyloseq <- function(name) {
    # create phyloseq object
    phylo <- import_biom(paste0("data/",name,".biom"), parallel = TRUE)
    # adjust colnames of the tax_table to map metaR specifications
    colnames(tax_table(phylo)) = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    phylo
}