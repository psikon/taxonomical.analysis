library(ProjectTemplate)
load.project()
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
registerDoParallel(cores = 20)
# in zeile 46 hinzufügen
# colnames(mdf)[3] <- "Abundance"
fixInNamespace(psmelt, pos = "package:phyloseq")

############################################################
################# initialize the data ######################
############################################################

# generate old dataset
path <- generate.biomFile("bacterial.old", data = get.DBConnection.old(get.metadata.list()))
phylo.old <- generate.phyloseq("bacterial.old")

# generate a biom file
path <- generate.biomFile("bacterial.new", data = get.DBConnection.new(get.metadata.list()))
# init the phyloseq object
phylo.new <- generate.phyloseq("bacterial.new")

#############################################################
######## Evaluation of different pipeline steps #############
#############################################################
load.project()

plot.taxaResolution(phylo.new,
                    file = "graphs/evaluation/taxa_res.abs.pdf",
                    absolute = TRUE, sep = TRUE, 
                    length_group1 = 4, length_group2 = 7,
                    title = "Taxonomical Resolution per samples \n(abs)")
plot.taxaResolution(phylo.new,
                    file = "graphs/evaluation/taxa_res.perc.pdf",
                    absolute = FALSE, sep = TRUE, 
                    length_group1 = 4, length_group2 = 7,
                    title = "Taxonomical Resolution per samples \n(percent)")

plot.groupedAbundance(phylo.new, 
                      file = "graphs/evaluation/abundance.abs.pdf",
                      absolute = TRUE, sep = TRUE,
                      length_group1 = 4, length_group2 = 7,
                      title = "Abundance in defined groups per sample\n(absolute)")
plot.groupedAbundance(phylo.new, 
                      file = "graphs/evaluation/abundance.perc.pdf",
                      absolute = FALSE, sep = TRUE,
                      length_group1 = 4, length_group2 = 7,
                      title = "Abundance in defined groups per sample\n(percent)")

plot.DBcount(data = get.DBConnection.new(get.metadata.list()) , 
             names = c("sample 60", "sample 64", "sample 66","sample 68", "sample 70", 
                       "sample 72", "sample 74", "sample 76", "sample 78", 
                       "sample 80", "sample 82"), 
             file = "graphs/evaluation/database_counts.pdf",
             sep = TRUE, length_group1 = 4, length_group2 = 7,
             title = "Hits in taxonomyReportDB per sample")

############################################################
################## General separations #####################
############################################################
load.project()

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

#############################################################
################### Rarefaction Curves ######################
#############################################################
load.project()

# plot rarefaction curves
plot.rareCurve(bakteria, stepsize = 20, file = "graphs/rarefaction/rareCurve.pdf")

# rarefaction curve for mariculture samples 70 - 82
plot.rareCurve(get.aqua(bakteria), stepsize = 20, 
               file = "graphs/rarefaction/mariculture.rareCurve.pdf")

# rarefaction curve for free living samples 60 - 68
plot.rareCurve(get.free(bakteria), stepsize = 20,
               file = "graphs/rarefaction/free_living.rareCurve.pdf")


###############################################################
########### Most abundant Taxa per Habitat/Sample ############
###############################################################
load.project()

# create the rarefied phyloseq object from bakteria 
# corresponding to the minimum number of reads
rare.bak <- rarify.phyloseq(bakteria, rngseed = 1234, 
                            replace = TRUE, trimOTUs = TRUE)

# Plot the most abundant Taxa per Habitat with a rarefied dataset
plot.mostAbundant.habitat(rare.bak, level = "phylum", threshold = 0.01,
                          file = "graphs/most_abundant/habitat.phylum.pdf")
plot.mostAbundant.habitat(rare.bak, level = "class", threshold = 0.01,
                          file = "graphs/most_abundant/habitat.class.pdf")
plot.mostAbundant.habitat(rare.bak, level = "order", threshold = 0.01,
                          file = "graphs/most_abundant/habitat.order.pdf")
plot.mostAbundant.habitat(rare.bak, level = "family", threshold = 0.01,
                          file = "graphs/most_abundant/habitat.family.pdf")
plot.mostAbundant.habitat(rare.bak, level = "genus", threshold = 0.01,
                          file = "graphs/most_abundant/habitat.genus.pdf")

# plot most abundant taxa per sample with a rarefied dataset
plot.mostAbundant.sample(rare.bak, level = "phylum", threshold = 0.01,
                         file = "graphs/most_abundant/sample.phylum.pdf")
plot.mostAbundant.sample(rare.bak, level = "class", threshold = 0.01,
                         file = "graphs/most_abundant/sample.class.pdf")
plot.mostAbundant.sample(rare.bak, level = "order", threshold = 0.01,
                         file = "graphs/most_abundant/sample.order.pdf")
plot.mostAbundant.sample(rare.bak, level = "family", threshold = 0.01,
                         file = "graphs/most_abundant/sample.family.pdf")
plot.mostAbundant.sample(rare.bak, level = "genus", threshold = 0.01,
                         file = "graphs/most_abundant/sample.genus.pdf")

#############################################################
################# Richness per Sample #######################
#############################################################
load.project()

# plot richness
plot.overview.richness(rare.bak, file = "graphs/richness/richness.samples.pdf")

# get a list of richniess indices per sample
richness <- get.richness(rare.bak) 
richness
#############################################################
################### Core Microbiome #########################
#############################################################
load.project()

# create Core Microbiome
all.core <- get.coreMicrobiome(tax_glom(rare.bak, "genus"))
free.core <- get.coreMicrobiome(tax_glom(get.free(rare.bak), "genus"))
aqua.core <- get.coreMicrobiome(tax_glom(get.aqua(rare.bak), "genus"))

# create lists of all OTUs in Core Microbiome
all.core.table <- get.coreTable(all.core, "lists/core/core.all.txt")
free.core.table <- get.coreTable(free.core, "lists/core/core.free.list.txt")
aqua.core.table <- get.coreTable(aqua.core, "lists/core/core.aqua.list.txt")

# make a Venn Diagram free living core microbiome vs. mariculture core microbiome
plot.coreMicrobiome(tax_glom(rare.bak, "genus"), file = "graphs/core/core.venn.tiff")
# plota rarefaction curve for core otus
plot.rarifyCoreOtus(rare.bak, n = 50, steps = 20, 
                    file = "graphs/core/core.rarefactionCurve.pdf",
                    title = "Development of Core OTUs \n(Rarefaction Curve)")

 # plot the content of the different core microbiome at phylum level
plot_bar(all.core, file = "graphs/core/all.core.phylum.pdf", 
         level = "phylum", title = "Core Microbiome at\n phylum level")
plot.taxonomy.graph(all.core,file = "graphs/core/all.core.content.pdf", level= "genus")

plot_bar(free.core, file = "graphs/core/free.core.phylum.pdf", 
         level = "phylum", title = "Core Microbiome of free living samples \nat phylum level")
plot.taxonomy.graph(free.core,file = "graphs/core/free.core.content.pdf",level= "genus")

plot_bar(aqua.core, file = "graphs/core/aqua.core.phylum.pdf", 
         level = "phylum", title = "Core Microbiome of mariculture samples \nat phylum level")
plot.taxonomy.graph(aqua.core,file = "graphs/core/aqua.core.content.pdf",level= "genus")

# create phyloseq objects containing only singletons
singleton.one <- get.singleton.phyloseq(tax_glom(rare.bak,"genus"),
                                        num_samples = 1,
                                        rm_samples = TRUE)
singleton.two <- get.singleton.phyloseq(tax_glom(rare.bak,"genus"),
                                        num_samples = 2,
                                        rm_samples = TRUE)

# plot the content of the different singletons at phylum level 
plot_bar(singleton.one, file = "graphs/core/single.one.phylum.pdf",
         level = "phylum", title = "Singletons on phylum level\n(occurence: max 1 sample)")
plot_bar(singleton.two, file = "graphs/core/single.two.phylum.pdf",
         level = "phylum", title = "Singletons on phylum level\n(occurence: max 2 samples")

# save the singletons as a defined data.frame and a tab seperated file
singleton.one.list <- get.singleton.list(phyloseq = singleton.one,
                                         file = "lists/singletons/singletons.one_sample.txt")
singleton.two.list <- get.singleton.list(phyloseq = singleton.two,
                                         file = "lists/singletons/singletons.two_sample.txt")

##################################################################
############### Ordination between samples/OTUs ##################
##################################################################
load.project()

# ordination distance between samples
plot.ordination.samples(rare.bak,  
                        file = "graphs/ordination/ordination.sample.phylum.MDS.pdf",
                        method = "NMDS", distance = "bray",
                        hits = 5, num_samples = 0.5, num_best = 5,
                        level = "phylum", title = "Ordination of Samples\n(phylum level)")

# ordination between OTUs
plot.ordination.OTUs(rare.bak, file = "graphs/ordination/ord.otus.phylum.pdf",
                     hits = 5, num_samples = 0.5, num_best = 5, 
                     level = "phylum", title = "Ordination of OTUs\n(phylum level)", 
                     facet = T, sep = 4)

##################################################################
########### Differential OTUs deseq, deseq2, edgeR ###############
##################################################################
load.project()

# deseq1
plot.differential.OTUs(deseq.differential.OTUs(rare.bak, 6, alpha = 0.73), 
                                   file = "graphs/differential/deseq.differential.pdf", 
                                   origin = "deseq")
# deseq2
plot.differential.OTUs(deseq2.differential.OTUs(rare.bak, alpha = 0.998), 
                       file = "graphs/differential/deseq2.differential.pdf", 
                       origin = "deseq")
# edgeR
plot.differential.OTUs(edgeR.differential.OTUs(rare.bak, alpha = 0.2),
                       file = "graphs/differential/edgeR.differential.pdf",
                       origin = "edgeR")

plot_heatmap(rm.underscore(tax_glom(rare.bak,"class")),
             sample.label = "SampleName",
             sample.order = "HoldingCondition",
             taxa.label = "class",)

##################################################################
######### Function Analysis with blast2GO ########################
##################################################################
load.project()

create.FastaFromTaxonomyReportDB(path = "data/functional", position = c(1,4,5),
                                 fasta.files = c("data/functional/sample60.classify.fasta",
                                                 "data/functional/sample68.classify.fasta",
                                                 "data/functional/sample70.classify.fasta"))




