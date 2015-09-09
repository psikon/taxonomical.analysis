source("src/initAnalysis.R")
# colnames(mdf)[3] <- "Abundance"
fixInNamespace(psmelt, pos = "package:phyloseq")


########################
### Eukaryotic Reads ###
########################

# make OTU counts relative
ef.euk <- transform_sample_counts(eukaryotic, function(OTU) OTU/sum(OTU)*100)

# Most Abundant >1%
ef.euk.gt1 <- prune_taxa(taxa_sums(ef.euk)/nsamples(ef.euk) >= 0.1, 
                         ef.euk)
ef.euk.gt1 <- tax_glom(ef.euk.gt1, "order")
ef.euk.gt1 <- rm.underscore(ef.euk.gt1)
colnames(tax_table(ef.euk.gt1))[4] <- "Order"
tax_table(ef.euk.gt1)["163134",4] <- "Centrarchiformes"
sample_data(ef.euk.gt1)$HoldingCondition <- c(rep("free-living",5),
                                              rep("mariculture",7))
test <- otu_table(ef.euk.gt1)[c("8113","8030"),]
sum(as.vector(test[1,1:5]))
sum(test[1,6:12])
# create initial plot

pl.ef.euk.gt1 <- plot_bar(ef.euk.gt1, x = "SampleName",
                          fill = "Order")
pl.ef.euk.gt1
pl.ef.euk.gt1 <- pl.ef.euk.gt1 + facet_grid(~ HoldingCondition, 
                                            scales="free_x") +
    xlab("Samples") + 
    ylab("Relative abundance [%]") +
    geom_bar(aes(color=Order, fill=Order), 
             stat="identity", position="stack") + 
    theme_bw() + ggtheme_bar 
# show plot 
pl.ef.euk.gt1

# create a new copy
ef.euk <- prune_taxa(taxa_sums(eukaryotic) > 1, eukaryotic)
ef.euk <- transform_sample_counts(ef.euk, function(OTU) OTU/sum(OTU)*100)

# select rare > 1%
ef.euk.lt1 <- prune_taxa(taxa_sums(ef.euk)/nsamples(ef.euk) < 0.1, 
                         ef.euk)
ef.euk.lt1 <- tax_glom(ef.euk.lt1, "phylum")
ef.euk.lt1 <- rm.underscore(ef.euk.lt1)
# change names for plotting
colnames(tax_table(ef.euk.lt1))[2] <- "Phylum"
sample_data(ef.euk.lt1)$HoldingCondition <- c(rep("free-living", 5),
                                              rep("mariculture", 7))
# remove chordata
ef.euk.lt1 <- rm.taxa(ef.euk.lt1,"119488")
# remove not classified
ef.euk.lt1 <- rm.taxa(ef.euk.lt1,"33208")

# create initial plot
pl.euk.lt1 <- plot_bar(ef.euk.lt1, x = "SampleName",
                          fill = "Phylum")
pl.euk.lt1
pl.euk.lt1 <- pl.euk.lt1 + facet_grid(~ HoldingCondition, 
                                            scales="free_x") +
    xlab("Samples") + 
    ylab("Relative abundance [%]") +
    geom_bar(aes(color=Phylum, fill=Phylum), 
             stat="identity", position="stack") + 
    theme_bw() + ggtheme_bar 
# show plot 
pl.euk.lt1

tiff("graphs/euk.barplot.tiff", compression = "lzw", ,
     width = 8, height = 9, units = "in", res = 600)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1)))
    print(pl.ef.euk.gt1, vp = viewport(layout.pos.row = 1, 
                                    layout.pos.col = 1))
    print(pl.euk.lt1, vp = viewport(layout.pos.row = 2, 
                                    layout.pos.col = 1))
dev.off()

jpeg("graphs/euk.barplot.jpg", width = 7.3, height = 8, 
     units = "in", res = 300)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1)))
    print(pl.ef.euk.gt1, vp = viewport(layout.pos.row = 1, 
                                       layout.pos.col = 1))
    print(pl.euk.lt1, vp = viewport(layout.pos.row = 2, 
                                       layout.pos.col = 1))
dev.off()

# parasites

otu_table(ef.euk.lt1)["935648",]
tax_table(ef.euk.lt1)[,2]

#######################
### Bacterial Reads ###
#######################

# create a copy
ef.bac.gt1 <- bacteria
# reduce to order level
ef.bac.gt1 <- tax_glom(bacteria, "order")
# remove underscore syntax
ef.bac.gt1 <- rm.underscore(ef.bac.gt1)
# change some parameter for plotting correctly
sample_data(ef.bac.gt1)$HoldingCondition <- c(rep("free-living", 4),
                                              rep("mariculture", 7))
colnames(tax_table(ef.bac.gt1))[4] <- "Order"
# transform to percentage
ef.bac.gt1 <- transform_sample_counts(ef.bac.gt1, function(OTU) OTU/sum(OTU)*100)
# extract content for correct values
content <- ef.bac.gt1
# select only dominant taxa 99%
ef.bac.gt1 <- prune_taxa(taxa_sums(ef.bac.gt1)/nsamples(ef.bac.gt1) > 1, 
                         ef.bac.gt1)

# create initial plot
pl.ef.bac.gt1 <- plot_bar(ef.bac.gt1, x = "SampleName",
                          fill = "Order")
pl.ef.bac.gt1
pl.ef.bac.gt1 <- pl.ef.bac.gt1 + facet_grid(~ HoldingCondition, 
                                            scales="free_x") +
    xlab("") + 
    ylab("Relative abundance [%]") +
    geom_bar(aes(color=Order, fill=Order), 
             stat="identity", position="stack") +
    theme_bw() + ggtheme_bar 
# show plot 
pl.ef.bac.gt1

unique(tax_table(content))
otu_table(content)["",c("70","76","78","80","82")]
mean(otu_table(content)["28901",c("70","76","78","80","82")])
sd(otu_table(content)["53413",c("70","76","78","80","82")])
sum(otu_table(content)[c("658062", "107709", "28216", "976",
                     "186801", "1235990", "1131287", "932037",
                     "91061", "28211", "1239", "1236", "1224"),"60"])

# create copy
ef.bac.lt1 <- bacteria
# reduce to class level
ef.bac.lt1 <- tax_glom(bacteria, "phylum")
# remove underscore syntax
ef.bac.lt1 <- rm.underscore(ef.bac.lt1)
sample_data(ef.bac.lt1)$HoldingCondition <- c(rep("free-living", 4),
                                              rep("mariculture", 7))
colnames(tax_table(ef.bac.lt1))[2] <- "Phylum"
# transform to percentage
ef.bac.lt1 <- transform_sample_counts(ef.bac.lt1, function(OTU) OTU/sum(OTU)*100)
# extract content 
content <- ef.bac.lt1
# select only taxa with < 1% abundance
ef.bac.lt1 <- prune_taxa(taxa_sums(ef.bac.lt1)/nsamples(ef.bac.lt1) < 1, 
                         ef.bac.lt1)
ef.bac.lt1
# create initial plot
pl.ef.bac.lt1 <- plot_bar(ef.bac.lt1, x = "SampleName",
                          fill = "Phylum")
pl.ef.bac.lt1
pl.ef.bac.lt1 <- pl.ef.bac.lt1 + facet_grid(~ HoldingCondition, 
                                            scales="free_x") +
    xlab("") + 
    ylab("Relative abundance [%]") +
    geom_bar(aes(color=Phylum, fill=Phylum), 
             stat="identity", position="stack") +
    theme_bw() + ggtheme_bar 
# show plot 
pl.ef.bac.lt1

tiff("graphs/bac.barplot.tiff", compression = "lzw", ,
     width = 7.3, height = 8, units = "in", res = 600)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1)))
    print(pl.ef.bac.gt1, vp = viewport(layout.pos.row = 1, 
                                       layout.pos.col = 1))
    print(pl.ef.bac.lt1, vp = viewport(layout.pos.row = 2, 
                                       layout.pos.col = 1))
dev.off()

jpeg("graphs/bac.barplot.jpg", width = 7.3, height = 8, 
     units = "in", res = 300)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 1)))
    print(pl.ef.bac.gt1, vp = viewport(layout.pos.row = 1, 
                                   layout.pos.col = 1))
    print(pl.ef.bac.lt1, vp = viewport(layout.pos.row = 2, 
                                   layout.pos.col = 1))
dev.off()

tax_table(content)
mean(otu_table(content)["848","60"])
mean(otu_table(content)["1164990","74"])

ef.euk.phy <- tax_glom(ef.euk, "phylum")
ef.euk.phy <- rm.underscore(ef.euk.phy)
colnames(tax_table(ef.euk.phy))[2] <- "Phylum"
pl.ef.euk.phy <- plot_bar(ef.euk.phy, x = "SampleName",
                          fill = "Phylum")
pl.ef.euk.phy
pl.ef.euk.phy <- pl.ef.euk.phy + facet_grid(~ HoldingCondition, 
                                        scales="free_x") +
    xlab("") + 
    ylab("Relative abundance [%]") +
    geom_bar(aes(color=Phylum, fill=Phylum), 
             stat="identity", position="stack") +
    theme_bw() + ggtheme_bar 
# show plot 
pl.ef.euk.phy
jpeg("graphs/euk.phy.jpg", width = 7.3, height = 6, 
     units = "in", res = 300)
    pl.ef.euk.phy
dev.off()

