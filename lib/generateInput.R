#' generate.biomFile
#'
#' wrapper function converting taxonomyReportDB obeject to biom format
#' 
#'@param name           output file
#'@param data           list of txonomyReportDB's
#'
#'@return path to output file
#'@export
generate.biomFile <- function(name, 
                              data = list(sample60.old, sample62.old, 
                                          sample66.old, sample68.old, 
                                          sample70.old, sample76.old, 
                                          sample78.old, sample80.old, 
                                          sample82.old)){
    
    # generate biom object from databases and write it to disk
    biom <- to_biom(data)
    
    # write the biom file to disk
    write_biom(biom, biom_file = paste0("data/", name, ".biom"))
    # return path
    return(paste0("data/", name, ".biom"))
}
#' generate.phyloseq
#'
#' wrapper function importing biom file to phyloseq object
#' 
#'@param name   output file
#'
#'@return phyloseq object
#'@export
#'
generate.phyloseq <- function(name) {
    
    # create phyloseq object
    phyloseq <- import_biom(paste0("data/", name, ".biom"), 
                            parallel = TRUE)
    
    # adjust colnames of the tax_table to map metaR specifications
    colnames(tax_table(phyloseq)) = c("superkingdom", "phylum", 
                                      "class", "order", 
                                      "family", "genus", 
                                      "species")
    return(phyloseq)
}