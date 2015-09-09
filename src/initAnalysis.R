library(grid)
library(ggplot2)
library(ProjectTemplate)
load.project()
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
registerDoParallel(cores = 20)

phylo.bac <- generate.phyloseq("data/bacterial.biom")
# assign new SampleNames
sample_data(phylo.bac)$SampleName <- c("free1", "free2", "free3",
                                      "free4", "free5", "mari1",
                                      "mari2", "mari3", "mari4",
                                      "mari5", "mari6", "mari7") 
bacteria <- subset_taxa(phylo.bac, superkingdom == 'k__Bacteria')
# remove contamination
bacteria <- subset_samples(bacteria, SampleName!="free4")
bacteria <- rm.taxa(bacteria, "2")
bacteria


phylo.euk <- generate.phyloseq("data/eukaryotic.biom")
sample_data(phylo.euk)$SampleName <- c("free1", "free2", "free3",
                                       "free4", "free5", "mari1",
                                       "mari2", "mari3", "mari4",
                                       "mari5", "mari6", "mari7") 
eukaryotic <- subset_taxa(phylo.euk, superkingdom == 'k__Eukaryota')
eukaryotic

