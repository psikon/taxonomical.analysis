#' generate.biomFile
#'
#' wrapper function converting taxonomyReportDB obeject to biom format
#' 
#'@param name           output file
#'@param data           list of txonomyReportDB's
#'
#'@return path to output file
#'@export
generate.biomFile <- function(file, data){
    
    # generate biom object from databases and write it to disk
    biom <- to_biom(data)
    
    # write the biom file to disk
    write_biom(biom, biom_file = file)
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
generate.phyloseq <- function(file) {
    
    # create phyloseq object
    phyloseq <- import_biom(file, parallel = TRUE)
    
    # adjust colnames of the tax_table to map metaR specifications
    colnames(tax_table(phyloseq)) = c("superkingdom", "phylum", 
                                      "class", "order", 
                                      "family", "genus", 
                                      "species")
    return(phyloseq)
}