remove_taxa <- function(phylo, taxa) {
    phylo <- prune_species(!grepl(paste0("\\<",taxa,"\\>"),labels(otu_table(phylo))[[1]]), phylo)
    phylo
}

remove_Underscore <- function(phyloseq) {
    data <- phyloseq@tax_table@.Data
    if(any(grepl("__", data))) {
        data <- substr(data, 4, length(data))
        data[which(data == "")] <- "undefined"
    }
    phyloseq@tax_table@.Data <- data
    phyloseq
}


