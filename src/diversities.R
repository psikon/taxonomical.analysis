require(grid)
require(ggplot2)
source("src/initAnalysis.R")
# colnames(mdf)[3] <- "Abundance"
fixInNamespace(psmelt, pos = "package:phyloseq")

#################################
### bacterial alpha diversity ###
#################################

alpha_meas <- c("Observed", "Shannon", "Chao1")
bac.rar <- rarefy_even_depth(bacteria, rngseed = 1234, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
bac.rar <- prune_taxa(taxa_sums(bac.rar) > 0, bac.rar)

# calculate richness as table
table_rich <- estimate_richness(bac.rar, measures = alpha_meas)

sample_data(bac.rar)$HoldingCondition <- c(rep("free-living", 4),
                                           rep("mariculture", 7))
bac.rich <- plot_richness(bac.rar, x = "HoldingCondition", 
                          color = "HoldingCondition", 
                          measures = alpha_meas)
bac.rich <- bac.rich + guides(color = FALSE)
## pop out the original points
bac.rich$layers <- bac.rich$layers[-1]
bac.rich <- bac.rich + geom_point(size = 3, alpha = 0.5) +
    geom_boxplot(aes(x = HoldingCondition, y = value, 
                     fill = HoldingCondition), alpha = 0.1, 
                 outlier.size = 0) +
    ggtheme_alpha + xlab(NULL) +
    labs(fill = "breeding condition") +
    scale_color_manual(values = c("#6e8a3d", "#eeb422")) +
    scale_fill_manual(values = c("#6e8a3d", "#eeb422"))
bac.rich
ggsave("graphs/diversity/alpha_div_rarefy_even.pdf", width = 6.3, height = 4.3)

tiff("graphs/diversity/alpha_diversity.tiff", compression = "lzw",
     width = 6.3, height = 4.3, units = "in",res = 600)
    bac.rich
dev.off()

jpeg("graphs/diversity/alpha_diversity.jpg", width = 6.3,
     height = 4.3, units = "in", res = 300)
    bac.rich
dev.off()

# Observed medians
median(table_rich[1:4,1])
median(table_rich[5:11,1])
# chao1 medians
median(table_rich[1:4,2])
median(table_rich[5:11,2])
# Shannon medians
median(table_rich[1:4,4])
median(table_rich[5:11,4])
# calculate the evenness J = observed/log(shannon)
mean(table_rich[1:4,1])/log(mean(table_rich[1:4,4]))
mean(table_rich[5:11,1])/log(mean(table_rich[5:11,4]))

################################
### bacterial beta diversity ###
################################

bac.beta <- transform_sample_counts(bacteria, function(x) 1e6*x/sum(x))
phylum.sum <- tapply(taxa_sums(bac.beta), tax_table(bac.beta)[, "phylum"], sum, na.rm = TRUE)
topPhyla <- names(sort(phylum.sum, TRUE))[1:8]
p1 <- prune_taxa((tax_table(bac.beta)[, "phylum"] %in% topPhyla), bac.beta)

p1.ord <- ordinate(p1, "MDS", "bray")
pp1 <- plot_ordination(p1, p1.ord,
                       type = "samples", 
                       color = "HoldingCondition") +
    geom_point(size = 5) +
    scale_colour_brewer(type = "qual", palette = "Set1")
pp1

ggsave("doc/fig/nmds_bray_all.pdf", width = 6.3, height = 3.6)


tiff("graphs/final/nmds_bray_ordination.tiff", compression = "lzw", ,
     width = 6.3, height = 3.6, units = "in",res = 600)
ord
dev.off()

##################################
### eukaryotic alpha diversity ###
##################################

alpha_meas <- c("Observed", "Shannon", "Chao1")
euk.rar <- rarefy_even_depth(eukaryotic, rngseed = 1234, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
euk.rar <- prune_taxa(taxa_sums(euk.rar) > 0, euk.rar)

# calculate richness as table
table_rich <- estimate_richness(euk.rar, measures = alpha_meas)

sample_data(euk.rar)$HoldingCondition <- c(rep("free-living", 5),
                                           rep("mariculture", 7))
euk.rich <- plot_richness(euk.rar, x = "HoldingCondition", 
                          color = "HoldingCondition", 
                          measures = alpha_meas)
euk.rich <- euk.rich + guides(color = FALSE)
## pop out the original points
euk.rich$layers <- euk.rich$layers[-1]
euk.rich <- euk.rich + geom_point(size = 3, alpha = 0.5) +
    geom_boxplot(aes(x = HoldingCondition, y = value, 
                     fill = HoldingCondition), alpha = 0.1, 
                 outlier.size = 0) +
    ggtheme_alpha + xlab(NULL) +
    labs(fill = "breeding condition") +
    scale_color_manual(values = c("#6e8a3d", "#eeb422")) +
    scale_fill_manual(values = c("#6e8a3d", "#eeb422"))
euk.rich
ggsave("graphs/diversity/alpha_div_rarefy_even.euk.pdf", width = 6.3, height = 4.3)

tiff("graphs/diversity/alpha_diversity.euk.tiff", compression = "lzw",
     width = 6.3, height = 4.3, units = "in",res = 600)
euk.rich
dev.off()

#################################
### eukaryotic beta diversity ###
#################################



## Ordination plot ####





