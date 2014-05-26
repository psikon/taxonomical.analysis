library(ProjectTemplate)
load.project()
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
registerDoParallel(cores=20)

# load all databases and generate a biom file with it
path <- generateBiomFile("bacterial")
# init the phyloseq object
phylo <- generatePhyloseq("bacterial")

sample_data(phylo)


#"genus"

# subset the phyloseq object based on taxonomy expressions
fungi <- subset_taxa(phylo, superkingdom == "k__Eukaryota")
# only artefacts
#vir <- subset_taxa(phylo, superkingdom == "k__Viruses")
bakteria <- subset_taxa(phylo, superkingdom == 'k__Bacteria')
# remove contamination
bakteria <- remove_taxa(bakteria,"2")
bakteria

merged <- merge_samples(bakteria, group="pool_info")
pdf("graphs/test.pdf")
    plot_bar(tax_glom(merged,"class"), fill="class")
dev.off()

# Overview Phylum
phylum <- tax_glom(bakteria,taxrank="phylum")
phylum <- prune_taxa(taxa_sums(phylum) > 10, phylum)
phylum <- merge_samples(phylum, group = "pool_info")
pdf("graphs/test.pdf")
    p = plot_bar(tax_glom(merged,"family"), fill="family")
    p + facet_wrap(~Environment)
dev.off()
plot_overview_bar(phylum,level="phylum",seperator="Environment")
otu_phylum <- otu_table(phylum)
nrow(otu_phylum)

# Overview Class
class <- tax_glom(bakteria, taxrank="class")
otu_class <- otu_table(class)
nrow(otu_class)
class <- prune_taxa(taxa_sums(class) > 10, class)
plot_overview_bar(family,level="class",seperator="Environment")

# Overview Order
order <- tax_glom(bakteria, taxrank="order")
otu_order <- otu_table(order)
nrow(otu_order)
order <- prune_taxa(taxa_sums(order) > 10, order)
plot_overview_bar(family,level="order",seperator="Environment")

# Overview family
family <- tax_glom(bakteria,taxrank="family")
otu_family <- otu_table(family)git pull

nrow(otu_family)
family <- prune_taxa(taxa_sums(family) > 10, family)
plot_overview_bar(family,level="family",seperator="Environment")
plot_bar(phylum, facet_grid=~Environment)
# Overview genus
genus <- tax_glom(bakteria, taxrank="genus")
otu_genus <- otu_table(genus)
nrow(otu_genus)
genus <- prune_taxa(taxa_sums(genus) > 10, genus)
plot_overview_bar(genus,level="genus",seperator="Environment")

# Overview species
species <- tax_glom(bakteria, taxrank="species")
otu_table(species)
nrow(otu_table(species))
species <- prune_taxa(taxa_sums(species) > 5, species)
plot_overview_bar(species,level="species",seperator="Environment")




