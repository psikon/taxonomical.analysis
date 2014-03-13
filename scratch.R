library(ProjectTemplate)
load.project()
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
registerDoParallel(cores = 20)

#######################
# init the data       #
#######################

# load all databases and generate a biom file with it
path <- generateBiomFile("bacterial")
# init the phyloseq object
phylo <- generatePhyloseq("bacterial")

#############################################################

#######################
# General separations #
#######################

# subset the phyloseq object based on taxonomy expressions
fungi <- subset_taxa(phylo, superkingdom == "k__Eukaryota")
# only artefacts
vir <- subset_taxa(phylo, superkingdom == "k__Viruses")
bakteria <- subset_taxa(phylo, superkingdom == 'k__Bacteria')
# remove contamination
bakteria <- remove_taxa(bakteria, "2")
bakteria

#############################################################

#############################################
# Generating Overview plots for every level #
#############################################

# Overview Phylum #

# agglomerate the various taxa to map all to phylum level
phylum <- tax_glom(bakteria, taxrank = "phylum")
# remove taxa were all samples have less than 10 reads 
phylum <- prune_taxa(taxa_sums(phylum) > 10, phylum)
# calculate relative abundance
rel_phylum = transform_sample_counts(phylum, function(x) x/sum(x))
# create overview plot
pdf("graphs/overview/phylum.pdf")
    rel_phylum <- remove_Underscore(rel_phylum)
    p <- plot_bar(rel_phylum, x = "SampleName", y = "Abundance", fill = "phylum",
                  title = "Overview at phylum level")
    p$labels$y <- "Relative Abundance (in percent)" 
    p + geom_bar(aes(color = phylum, fill = phylum), 
                 stat = 'identity', position = 'stack') + facet_grid(. ~Environment, 
                                                                     scales = "free", 
                                                                     space = "free")
dev.off()
# save otu_table
otu_phylum <- otu_table(phylum)

# Overview Class #

class <- tax_glom(bakteria, taxrank = "class")
class <- prune_taxa(taxa_sums(class) > 10, class)
rel_class = transform_sample_counts(class, function(x) x/sum(x))
pdf("graphs/overview/class.pdf")
    rel_class <- remove_Underscore(rel_class)
    p <- plot_bar(rel_class, x = "SampleName", y = "Abundance", fill = "class",
                  title = "Overview at class level")
    p$labels$y <- "Relative Abundance (in percent)" 
    p + geom_bar(aes(color = class, fill = class), 
                 stat = 'identity', position = 'stack') + facet_grid(. ~Environment, 
                                                                     scales = "free", 
                                                                     space = "free")
dev.off()
otu_class <- otu_table(rel_class)

# Overview Order #

order <- tax_glom(bakteria, taxrank = "order")
order <- prune_taxa(taxa_sums(order) > 10, order)
rel_order = transform_sample_counts(order, function(x) x/sum(x))
pdf("graphs/overview/order.pdf")
    rel_order <- remove_Underscore(rel_order)
    p <- plot_bar(rel_order, x = "SampleName", y = "Abundance", fill = "order",
                  title = "Overview at order level")
    p$labels$y <- "Relative Abundance (in percent)" 
    p + geom_bar(aes(color = order, fill = order), 
                 stat = 'identity', position = 'stack') + facet_grid(. ~Environment, 
                                                                     scales = "free", 
                                                                     space = "free")
dev.off()
otu_order <- otu_table(rel_order)

# Overview family #

family <- tax_glom(bakteria, taxrank = "family")
family <- prune_taxa(taxa_sums(family) > 10, family)
rel_family = transform_sample_counts(family, function(x) x/sum(x))
pdf("graphs/overview/family.pdf")
    rel_family <- remove_Underscore(rel_family)
    p <- plot_bar(rel_family, x = "SampleName", y = "Abundance", fill = "family",
                  title = "Overview at family level")
    p$labels$y <- "Relative Abundance (in percent)" 
    p + geom_bar(aes(color = family, fill = family), 
                 stat = 'identity', position = 'stack') + facet_grid(. ~Environment, 
                                                                     scales = "free", 
                                                                     space = "free")
dev.off()
otu_family <- otu_table(rel_family)

# Overview genus #

genus <- tax_glom(bakteria, taxrank = "genus")
genus <- prune_taxa(taxa_sums(genus) > 10, genus)
rel_genus = transform_sample_counts(genus, function(x) x/sum(x))
pdf("graphs/overview/genus.pdf")
    rel_genus <- remove_Underscore(rel_genus)
    p <- plot_bar(rel_genus, x = "SampleName", y = "Abundance", fill = "genus",
                  title = "Overview at genus level")
    p$labels$y <- "Relative Abundance (in percent)" 
    p + geom_bar(aes(color = genus, fill = genus), 
                stat = 'identity', position = 'stack') + facet_grid(. ~Environment, 
                                                                    scales = "free", 
                                                                    space = "free")
dev.off()
otu_genus <- otu_table(rel_genus)

# Overview species #

species <- tax_glom(bakteria, taxrank = "species")
species <- prune_taxa(taxa_sums(species) > 10, species)
rel_species = transform_sample_counts(species, function(x) x/sum(x))
pdf("graphs/overview/species.pdf")
    rel_species <- remove_Underscore(rel_species)
    p <- plot_bar(rel_species, x = "SampleName", y = "Abundance", fill = "species",
                  title = "Overview at species level")
    p$labels$y <- "Relative Abundance (in percent)" 
    p + geom_bar(aes(color = species, 
                     fill = species), 
                 stat = 'identity', 
                 position = 'stack') + facet_grid(. ~Environment, 
                                                  scales = "free", 
                                                  space = "free") + theme(legend.text = element_text(size = 9),
                                                                          legend.key.size = unit(0.5,"cm"))
dev.off()
otu_species <- otu_table(rel_species)

########################################################################



