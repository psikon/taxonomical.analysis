getMetadata <- function(){
    source("src/metadata.R")
}

getConnections <- function() {
    source("src/connections.R")
}

generateBiomFile <- function(name, data = list(sample60.old, sample62.old, sample66.old,
                                               sample68.old, sample70.old, sample76.old, 
                                               sample78.old, sample80.old, sample82.old)){
    # generate biom object from databases and write it to disk
    biom <- to_biom(data)
    write_biom(biom, biom_file = paste0("data/",name,".biom"))
    return(paste0("data/",name,".biom"))
}

generatePhyloseq <- function(name) {
    # create phyloseq object
    phylo <- import_biom(paste0("data/",name,".biom"), parallel = TRUE)
    # adjust colnames of the tax_table to map metaR specifications
    colnames(tax_table(phylo)) = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    phylo
}