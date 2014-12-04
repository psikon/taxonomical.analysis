library(ProjectTemplate)
load.project()
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
registerDoParallel(cores = 20)

phylo.new <- generate.phyloseq("bacterial.new")
# subset the phyloseq object based on taxonomy expressions
fungi <- subset_taxa(phylo.new, superkingdom == "k__Eukaryota")
# only artefacts
virus <- subset_taxa(phylo.new, superkingdom == "k__Viruses")

bakteria <- subset_taxa(phylo.new, superkingdom == 'k__Bacteria')
# remove contamination
bakteria <- rm.taxa(bakteria, "2")
# filter less abundant taxa
#bakteria <- prune_taxa(taxa_sums(bakteria) > 3, bakteria)
bakteria
