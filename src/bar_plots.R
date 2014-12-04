source("src/initAnalysis.R")
# colnames(mdf)[3] <- "Abundance"
fixInNamespace(psmelt, pos = "package:phyloseq")

## barplots ####
ggtheme_bar <- theme(
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 8),
    plot.margin = unit(c(0.025,0.025,.025,0.025), "npc")
)


# make OTU counts relative
ef <- transform_sample_counts(bakteria, function(OTU) OTU/sum(OTU)*100)

# Most Abundant >1%
ef.gt1 <- prune_taxa(taxa_sums(ef)/nsamples(ef) >= 0.1, ef)
ef.gt1 <- rm.underscore(ef.gt1)
# create initial plot
pl.ef.gt1 <- plot_bar(ef.gt1, x = "SampleName",fill = "phylum")

pl.ef.gt1 <- pl.ef.gt1 + facet_grid( ~ HoldingCondition, scales="free_x") +
    xlab("") +
    ylab("Relative abundance [%]") +
    geom_bar(aes(color=phylum, fill=phylum), 
             stat="identity", position="stack") +
    theme(strip.text=element_text(size=rel(1.25), face="italic"),
          legend.key.height=unit(1.0, "lines")) + ggtheme_bar + theme_bw()
# show plot 
pl.ef.gt1




ef.lt1 <- prune_taxa(taxa_sums(bakteria)/nsamples(bakteria) <= 0.1, bakteria)
