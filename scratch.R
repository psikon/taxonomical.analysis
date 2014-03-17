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
                                              space = "free") + theme(legend.text = element_text(size = 7),
                                                                      legend.key.size = unit(0.25,"cm"))
dev.off()
otu_species <- otu_table(rel_species)

########################################################################
load.project()
######################
# Most Abundant Taxa #
######################

# phylum level

# agglomerate all taxa at phylum level
phylum <- tax_glom(bakteria, taxrank = "phylum")
# calculate relative abundance
rel_phylum <- transform_sample_counts(phylum, function(x) x/sum(x))
# sort the taxa with decreasing relative abundance
top_phylum <-  sort(tapply(taxa_sums(rel_phylum),tax_table(rel_phylum)[, "phylum"], sum), TRUE)
# select only the taxa with an abundance of 10% or more
top_phylum <- top_phylum[which(top_phylum >= 0.2)]
# reconstruct phyloseq object
top_phylum <- transform_sample_counts(subset_taxa(phylum, phylum %in% names(top_phylum)), function(x) x/sum(x))

pdf("graphs/most_abundant/ma_phylum.pdf")
    top_phylum <- remove_Underscore(top_phylum)
    p <- plot_bar(top_phylum, x = "SampleName", y = "Abundance", fill = "phylum",
                title = paste0("Top ",nrow(otu_table(top_phylum))," taxa at phylum level"))
    p$labels$y <- "Relative Abundance (in percent)"
    p + geom_bar(aes(color = phylum, 
                     fill = phylum), 
                stat = 'identity', position = 'stack') + facet_grid(. ~Environment, 
                                                                    scales = "free", 
                                                                    space = "free") 
dev.off()

# class level

class <- tax_glom(bakteria, taxrank = "class")
rel_class <- transform_sample_counts(class, function(x) x/sum(x))
top_class <-  sort(tapply(taxa_sums(rel_class),tax_table(rel_class)[, "class"], sum),TRUE)
top_class <- top_class[which(top_class >= 0.2)]
top_class <- transform_sample_counts(subset_taxa(class, class %in% names(top_class)), 
                                      function(x) x/sum(x))
pdf("graphs/most_abundant/ma_class.pdf")
    top_class <- remove_Underscore(top_class)
    p <- plot_bar(top_class, x = "SampleName", y = "Abundance", fill = "class",
                title = "Top 5 taxa at class level")
    p$labels$y <- "Relative Abundance (in percent)"
    p + geom_bar(aes(color = class, 
                    fill = class), 
                stat = 'identity', position = 'stack') + facet_grid(. ~Environment, 
                                                                    scales = "free", 
                                                                    space = "free")
dev.off()

# order level

order <- tax_glom(bakteria, taxrank = "order")
rel_order <- transform_sample_counts(order, function(x) x/sum(x))
top_order <-  sort(tapply(taxa_sums(rel_order),tax_table(rel_order)[, "order"], sum),TRUE)
top_order <- top_order[which(top_order >= 0.2)]
top_order <- transform_sample_counts(subset_taxa(order, order %in% names(top_order)), 
                                     function(x) x/sum(x))
pdf("graphs/most_abundant/ma_order.pdf")
    top_order <- remove_Underscore(top_order)
    p <- plot_bar(top_order, x = "SampleName", y = "Abundance", fill = "order",
                  title = "Top 9 taxa at order level")
    p$labels$y <- "Relative Abundance (in percent)"
    p + geom_bar(aes(color = order, 
                    fill = order), 
                stat = 'identity', position = 'stack') + facet_grid(. ~Environment, 
                                                                    scales = "free", 
                                                                    space = "free")
dev.off()

# family level

family <- tax_glom(bakteria, taxrank = "family")
rel_family <- transform_sample_counts(family, function(x) x/sum(x))
top_family <-  sort(tapply(taxa_sums(rel_family),tax_table(rel_family)[, "family"], sum),TRUE)
top_family <- top_family[which(top_family >= 0.2)]
top_family <- transform_sample_counts(subset_taxa(family, family %in% names(top_family)), 
                                     function(x) x/sum(x))
pdf("graphs/most_abundant/ma_family.pdf")
    top_family <- remove_Underscore(top_family)
    p <- plot_bar(top_family, x = "SampleName", y = "Abundance", fill = "family",
                    title = "Top 10 taxa at family level")
    p$labels$y <- "Relative Abundance (in percent)"
    p + geom_bar(aes(color = family, 
                    fill = family), 
                 stat = 'identity', position = 'stack') + facet_grid(. ~Environment, 
                                                                     scales = "free", 
                                                                     space = "free")
dev.off()

# genus level

genus <- tax_glom(bakteria, taxrank = "genus")
rel_genus <- transform_sample_counts(genus, function(x) x/sum(x))
top_genus <-  sort(tapply(taxa_sums(rel_genus),tax_table(rel_genus)[, "genus"], sum),TRUE)
top_genus <- top_genus[which(top_genus >= 0.3)]
top_genus <- transform_sample_counts(subset_taxa(genus, genus %in% names(top_genus)), 
                                      function(x) x/sum(x))
pdf("graphs/most_abundant/ma_genus.pdf")
    top_genus <- remove_Underscore(top_genus)
    p <- plot_bar(top_genus, x = "SampleName", y = "Abundance", fill = "genus",
                title = "Top 10 taxa at genus level")
    p$labels$y <- "Relative Abundance (in percent)"
    p + geom_bar(aes(color = genus, 
                 fill = genus), 
                 stat = 'identity', position = 'stack') + facet_grid(. ~Environment, 
                                                                    scales = "free", 
                                                                    space = "free")
dev.off()

# species level

species <- tax_glom(bakteria, taxrank = "species")
rel_species <- transform_sample_counts(species, function(x) x/sum(x))
top_species <-  sort(tapply(taxa_sums(rel_species),tax_table(rel_species)[, "species"], sum),TRUE)
top_species <- top_species[which(top_species >= 0.2)]
top_species <- transform_sample_counts(subset_taxa(species, species %in% names(top_species)), 
                                     function(x) x/sum(x))
pdf("graphs/most_abundant/ma_species.pdf")
    top_species <- remove_Underscore(top_species)
    p <- plot_bar(top_species, x = "SampleName", y = "Abundance", fill = "species",
                title = "Top 10 taxa at species level")
    p$labels$y <- "Relative Abundance (in percent)"
    p + geom_bar(aes(color = species, 
                     fill = species), 
                stat = 'identity', position = 'stack') + facet_grid(. ~Environment, 
                                                                 scales = "free", 
                                                                 space = "free")
dev.off()

######################################################################################
load.project()
###############################
# Rarefication of the samples #
###############################

rare = rarefy_even_depth(bakteria)
plot_bar(rare)
