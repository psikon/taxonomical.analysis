library(ProjectTemplate)
load.project()
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
registerDoParallel(cores=20)

path <- generateBiomFile("bacterial")
phylo <- generatePhyloseq("bacterial")

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

