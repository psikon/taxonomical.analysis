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
bakteria <- prune_taxa(taxa_sums(bakteria) > 3, bakteria)
bakteria

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


