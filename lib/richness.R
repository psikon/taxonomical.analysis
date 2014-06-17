plot_richness_overview <- function(phyloseq, 
                                   file = "graphs/richness/richness.overview.pdf",
                                   measures = c("Observed", "Shannon", "Simpson","ACE", "Chao1")) {
    # adjust theme options (see utils.R)
    ggtheme_alpha <- get_richness_theme()
    # rarefy to even depth
    prar <- rarefy_even_depth(phyloseq, rngseed = 1234, replace = TRUE, trimOTUs = TRUE)
    # remove taxa with all abundance values == 0
    prar <- prune_taxa(taxa_sums(prar) > 0, prar)
    # create richness graph
    pp_rich <- plot_richness(prar, color = "HoldingCondition", measures = measures)
    # disable legend
    pp_rich <- pp_rich + guides(color = FALSE)
    ## pop out the original points
    pp_rich$layers <- pp_rich$layers[-1]
    # adjust graph theme
    pp_rich <- pp_rich + geom_point(size = 3, alpha = 0.5) +
        geom_boxplot(aes(y = value, fill = HoldingCondition), alpha = 0.1, outlier.size = 0) +
        ggtheme_alpha +
        xlab("Samples") +
        labs(fill = "Holding condition") +
        scale_color_manual(values = c("#6e8a3d", "#eeb422")) +
        scale_fill_manual(values = c("#6e8a3d", "#eeb422"))
    # plot to file
    pdf(file)
        plot(pp_rich) 
    dev.off()
}

get_richness <- function(phyloseq,
                         measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson")) {
    data <- rarefy_even_depth(phyloseq, rngseed = 1234, replace = TRUE, trimOTUs = TRUE)
    data <- estimate_richness(data, measures = measures)
    if (nrow(sample_data(phyloseq))>2) rownames(data) <- as.character(sample_data(phyloseq)$SampleName)
    data <- round(data,2)
    data
}
