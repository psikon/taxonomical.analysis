source("src/initAnalysis.R")
# create initial bacteria and eukaryota phyloseq object
eukaryota <- subset_taxa(generate.phyloseq("data/eukaryotic.biom"), superkingdom == 'k__Eukaryota')

##########################
### Eukaryota Overview ###
##########################

## Mean and range of remaining samples
mean(colSums(otu_table(eukaryota)))
sum(colSums(otu_table(eukaryota)))

which(apply(tax_table(rm.underscore(eukaryota)), 1, function(x) {
    sum(str_count(x,"not classified")) == length(x)-1
}))
# remove OTUs onlxy mapping to eukaryota
eukaryota <- rm.taxa(eukaryota, "2759")
# remove underscore
eukaryota <- rm.underscore(eukaryota)
# map all to phylum
eukaryota <- tax_glom(eukaryota, "phylum")
# remove OTUs with only one hit for all samples
eukaryota <- prune_taxa(taxa_sums(eukaryota) > 1, eukaryota)
# transform otu_table to percent values
eukaryota <- transform_sample_counts(eukaryota, function(OTU) OTU/sum(OTU)*100)

## Classify phyla in decreasing order of abundance
unique_phyla <- unique(tax_table(eukaryota)[, rank_names(eukaryota)[2]])
phyla <- tax_glom(eukaryota, taxrank = rank_names(eukaryota)[2])
unique_phyla <- as.vector(tax_table(phyla)[order(taxa_sums(phyla), 
                                                 decreasing = TRUE),
                                           rank_names(eukaryota)[2]])
length(unique_phyla)
unique_phyla

eukaryota_c <- round(sort(taxa_sums(eukaryota), decreasing = TRUE)/sum(taxa_sums(eukaryota))*100, 4)
eukaryota_c

# most abundant family level

euk <- transform_sample_counts(eukaryota, function(OTU) OTU/sum(OTU)*100)
euk <- tax_glom(euk, "order")
ma <- prune_taxa(taxa_sums(euk)/nsamples(euk) >= 0.1, euk)
unique <- unique(tax_table(ma)[, rank_names(ma)[]])
unique
euk_c <- round(sort(taxa_sums(euk), 
                   decreasing = TRUE)/sum(taxa_sums(euk))*100, 2)
length(euk_c)
unique[,4]

#########################
### Bacteria Overview ###
#########################

bakteria <- subset_taxa(generate.phyloseq("data/bacterial.biom"), superkingdom == 'k__Bacteria')

# 865 - OTUs
bakteria
# sample 66 has only 6111 reads -> filter
colSums(otu_table(bakteria))
bakteria <- subset_samples(bakteria, SampleName!="ef_free_66")
bakteria
with.2 <- sum(colSums(otu_table(bakteria)))
# remove contamination
bakteria <- rm.taxa(bakteria, "2")
# 930 - OTUs
bakteria
## Total number of valid sequence reads
without.2 <- sum(colSums(otu_table(bakteria)))
# Total annotated reads
with.2
# usable reads
without.2
# lost reads
with.2 - without.2


## Mean and range of remaining samples
summary(colSums(otu_table(bakteria)))

# remove underscore
bakteria <- rm.underscore(bakteria)
## Classify phyla in decreasing order of abundance
unique_phyla <- unique(tax_table(bakteria)[, rank_names(bakteria)[2]])
phyla <- tax_glom(bakteria, taxrank = rank_names(bakteria)[2])
unique_phyla <- as.vector(tax_table(phyla)[order(taxa_sums(phyla), decreasing = TRUE), rank_names(bakteria)[2]])
length(unique_phyla)
paste0(unique_phyla, collapse = ", ")
unique_phyla

## percentages of phyla overall
phyla <- round(sort(taxa_sums(phyla), decreasing = TRUE)/sum(taxa_sums(phyla))*100, 4)
names(phyla) <- unique_phyla
phyla

ma <- transform_sample_counts(bakteria, function(OTU) OTU/sum(OTU)*100)
ma <- tax_glom(ma, "order")
ma <- prune_taxa(taxa_sums(ma)/nsamples(ma) >= 1, ma)
ma_c <- round(sort(taxa_sums(ma), 
                   decreasing = TRUE)/sum(taxa_sums(ma))*100, 4)
round(sort(taxa_sums(ma), decreasing = TRUE)/sum(taxa_sums(ma))*100, 4)



