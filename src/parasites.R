source("src/initAnalysis.R")
# colnames(mdf)[3] <- "Abundance"
fixInNamespace(psmelt, pos = "package:phyloseq")

## barplots ####
ggtheme_bar <- theme(
    legend.text = element_text(family = "Times", size = 10),
    legend.title = element_text(family = "Times", size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(hjust = 1, angle = 45),
    strip.text = element_text(size = 10),
    plot.margin = unit(c(0.025,0.025,.025,0.025), "npc")
)

###############################
### less abundant eukaryota ###
###############################

# first select only eukaryota
eukaryota <- subset_taxa(generate.phyloseq("data/eukaryotic.biom"), 
                         superkingdom == 'k__Eukaryota')
eukaryota <- rm.taxa(eukaryota, "2759")
# remove underscore
eukaryota <- rm.underscore(eukaryota)

# get rare Eukaryota on phylum level
eukaryota <- tax_glom(eukaryota, "phylum")
# remove OTUs with only one hit for all samples
eukaryota <- prune_taxa(taxa_sums(eukaryota) > 1, eukaryota)
# transform sample counts 
eukaryota <- transform_sample_counts(eukaryota, function(OTU) OTU/sum(OTU)*100)
# select only the rare phyla 
eukaryota
euk.lt1 <- prune_taxa(taxa_sums(euk.lt1)/nsamples(euk.lt1) < 0.1, 
                      euk.lt1)
# create rare phyla overview 
unique_rare <- unique(tax_table(euk.lt1)[, rank_names(euk.lt1)[2]])
unique_rare
rare_c <- round(sort(taxa_sums(eukaryota), decreasing = TRUE)/sum(taxa_sums(eukaryota))*100, 4)[-1]
rare_c
# split rare phyla into groups

# group 1 - parasites

# select all eukaryota and select only up to class level
parasites <- tax_glom(eukaryota, "class")
# show tax table
tax_table(parasites)

# create parasite table for harry palm 
pl <- subset_taxa(parasites, phylum == "Platyhelminthes")
an <- subset_taxa(parasites, phylum == "Annelida")
ac <- subset_taxa(parasites, phylum == "Acanthocephala")
ne <- subset_taxa(parasites, phylum == "Nematoda")
cr <- subset_taxa(parasites, phylum == "Arthropoda")
# combine the results
parasites <- merge_phyloseq(pl, an, ac, ne, cr)
# extract the abundance
abundance <- otu_table(parasites)
# extract the taxonomy ids from the rownames
tax_ids <- rownames(abundance)
# search for linages
taxa <- taxonDB(tax_ids)
# select only interesting linages
selected <- c(1,2,6,7,8,9,14)
selected_linages <- taxa[selected]
selected_abundance <- rowSums(abundance[selected,])
x <- c(selected_linages[[1]]@Lineage@ScientificName[9],
       selected_linages[[2]]@Lineage@ScientificName[8],
       selected_linages[[3]]@Lineage@ScientificName[9],
       selected_linages[[4]]@Lineage@ScientificName[9],
       selected_linages[[5]]@Lineage@ScientificName[9],
       selected_linages[[6]]@Lineage@ScientificName[13],
       selected_linages[[7]]@Lineage@ScientificName[13])
names(selected_abundance) <- x
tbl <- c(selected_abundance[1:3], selected_abundance[4]+ selected_abundance[5],
         selected_abundance[6] + selected_abundance[7])
tbl

# create group 1 for second paper

# select all eukaryota and select only up to class level
parasites <- tax_glom(eukaryota, "phylum")
parasites <- rm.taxa(parasites, "2759")
parasites <- rm.underscore(parasites)
# remove OTUs with only one hit for all samples
parasites <- prune_taxa(taxa_sums(parasites) > 1, parasites)
# make counts relative
parasites <- transform_sample_counts(parasites, function(OTU) OTU/sum(OTU)*100)
# create parasite table for harry palm 
pl <- subset_taxa(parasites, phylum == "Platyhelminthes")
an <- subset_taxa(parasites, phylum == "Annelida")
ac <- subset_taxa(parasites, phylum == "Acanthocephala")
ne <- subset_taxa(parasites, phylum == "Nematoda")
nem <- subset_taxa(parasites, phylum == "Nemertea")
# combine the results
parasites <- merge_phyloseq(pl, an, ac, ne, nem)

# ATTANTION: counts are not right. Only the distributions are correct. 
# Only for graphical interpretation  
# create initial barplot
par_pl <- plot_bar(parasites, x = "SampleName", fill = "phylum")
# tweak optical parameter of plot
par_pl <- par_pl + facet_grid(~ HoldingCondition, 
                                scales="free_x") +
    xlab("Samples") + 
    ylab("Relative abundance [%]") +
    geom_bar(aes(color=phylum, fill=phylum), 
             stat="identity", position="stack") + 
    theme_bw() + ggtheme_bar 
# show plot
par_pl
# and save it
ggsave("graphs/less_abundant/parasites.pdf", width = 6.3, height = 4.6)

# determine mean counts for 2 breeding conditions
round(mean(colSums(otu_table(parasites)[,1:5])),4)
round(mean(colSums(otu_table(parasites)[,6:11])), 4)
