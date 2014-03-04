library(ProjectTemplate)
load.project()
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
registerDoParallel(cores=20)
# get all metadata objects
source("src/metadata.R")
# connect to all databses
source("src/connections.R")


# generate biom object from databases and write it to disk
biom <- to_biom(list(sample60, sample62, 
                     sample66, sample68, sample70,
                     sample76, sample78, 
                     sample80, sample82))
write_biom(biom, biom_file = "data/bacterial.biom")

# create phyloseq object
phylo <- import_biom("data/bacterial.biom", parallel = TRUE)
# adjuxt colnames of the tax_table to map metaR specifications
colnames(tax_table(phylo)) = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

# subset the phyloseq object based on taxonomy expressions
euk <- subset_taxa(phylo, superkingdom == "k__Eukaryota")
fun <- subset_taxa(phylo, phylum == "p__Opisthokonta")
vir <- subset_taxa(phylo, superkingdom == "k__Viruses")
bakteria <- subset_taxa(phylo, superkingdom == 'k__Bacteria')
bak <- prune_taxa(taxa_sums(bak) > 2, bak)
sum(count(get_taxa(phylo,"80"))$x*count(get_taxa(phylo,"80"))$freq)
sum(count(get_taxa(phylo,"62"))$x*count(get_taxa(phylo,"62"))$freq)
bak

bakteria <- remove_taxa(bakteria,"2")
bak
otu_table(bakteria)
taxon(1029823)
taxon(964)
taxon(85698)
taxon(672)
taxon(662)
taxon(645)
taxon(346)  
taxon(33917)
taxon(1773)
taxon(231049)
taxon(1236)
taxon(1224)
taxon(1158459)
taxon(1110389)
taxon(38293)
pdf("graphs/bacteria.pdf")
    plot_overview_bar(phyloseq=bakteria, level="family",seperator="Environment")
dev.off()



plot_heatmap(bak)
pdf("graphs/initial_heat.pdf")
plot_heatmap(bak, low = "#000033", high = "#FF3300")
dev.off()
sample60

