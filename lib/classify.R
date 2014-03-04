
intersectTaxa <- function(phylo) {
    intersect(otu_table(phylo)[,1],otu_table(phylo)[,2])
}







