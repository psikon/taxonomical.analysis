library(ProjectTemplate)
load.project()
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
registerDoParallel(cores = 20)
# in zeile 46 hinzufügen
# colnames(mdf)[3] <- "Abundance"
fixInNamespace(psmelt, pos = "package:phyloseq")


#######################
# init the data       #
#######################

# load all databases and metadata
getMetadata()
getConnections()
# generate old dataset
path <- generateBiomFile("bacterial.old")
phylo.old <- generatePhyloseq("bacterial.old")

# generate a biom file
path <- generateBiomFile("bacterial.new", data = list(sample60.new, 
                                                      sample64.new,
                                                      sample68.new, sample70.new,
                                                      sample72.new, sample74.new,
                                                      sample76.new,
                                                      sample80.new, sample82.new))
# init the phyloseq object
phylo.new <- generatePhyloseq("bacterial.new")

#############################################################

# Evaluation of different pipeline steps

source("evaluation_report.R")

#######################
# General separations #
#######################

# subset the phyloseq object based on taxonomy expressions
fungi <- subset_taxa(phylo.new, superkingdom == "k__Eukaryota")
# only artefacts
virus <- subset_taxa(phylo.new, superkingdom == "k__Viruses")

bakteria <- subset_taxa(phylo.new, superkingdom == 'k__Bacteria')
# remove contamination
bakteria <- remove_taxa(bakteria, "2")
# filter less abundant taxa
#bakteria <- prune_taxa(taxa_sums(bakteria) > 3, bakteria)
bakteria

#############################################################
load.project()

# plot richness

# eliminate species without any abundance
bak <- prune_species(speciesSums(bakteria) > 0, bakteria)
bak

plot_richness_overview(merge_samples(bak,group="Environment"), file = "graphs/richness/overview.environment.pdf")
plot_richness_overview(bak,file="graphs/richness/overview.all_sample.pdf")
plot_richness_shannon(bak)
plot_richness_simpson(bak)
plot_richness_Chao(bak)
plot_richness_ACE(bak)
richness <- get_richness(bak) 

#############################################################
load.project()

# Core Microbioms

all_core <- get_core_microbiom(bakteria)
all_core.table <- get_core_table(all_core,"graphs/core microbiom/all_list.txt")
free_core <- get_core_microbiom(get_free(bakteria))
free_core.table <- get_core_table(free_core,"graphs/core microbiom/free_list.txt")
aqua_core <- get_core_microbiom(get_aqua(bakteria))
aqua_core.table <- get_core_table(aqua_core,"graphs/core microbiom/aqua_list.txt")
plot_core_venn(bakteria,"graphs/core microbiom/free_vs_aqua_venn")

pdf("graphs/all_core.pdf")
    plot_bar(tax_glom(all_core,"class"),fill="class")
dev.off()
pdf("graphs/aqua_core.pdf")
    plot_bar(tax_glom(aqua_core,"class"),fill="class")
dev.off()
# Singletons
single_one_sample <- get_singletons(phyloseq = bakteria,
                                    num_samples = 1,
                                    remove_samples = TRUE)
single_two_sample <- get_singletons(phyloseq = bakteria,
                                    num_samples = 2,
                                    remove_samples = TRUE)


get_singleton_list(phyloseq = single_one_sample,
                   file = "singletons.one_sample.txt")
get_singleton_list(phyloseq = single_two_sample,
                   file = "singletons.two_sample.txt")

##################################################################

# ordination distance between samples
ord.samples <- ordination_preprocess(bak, hits = 5, num_samples = 0.5, 
                                     level="phylum", num_best = 5)
plot_ordination_Samples(ord.samples, file = "graphs/ordination/ord.sample.phylum5.pdf",
                        level = "phylum", title = "Ordination of Samples by 5 best phyla")

# ordination between OTUs
ord.otus <- ordination_preprocess(bak, hits = 5, num_samples= 0.5, level = "class", num_best = 10)
plot_ordination_OTUs(ord.otus, file = "graphs/ordination/ord.otus.class10.pdf",
                     level = "class", title = "Ordination of OTUs by 10 best class", 
                     facet = T, sep = 5)
