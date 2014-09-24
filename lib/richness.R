#' plot.richness
#'
#' wrapper function for phyloseq plot_richness function
#' 
#'@description create a plot containing all desired diversity 
#'             indices with alpha diversity measures on the 
#'             y axis and the samples on the x axis, each 
#'             indices has its own plot which will be sticked 
#'             together in one graphic
#' 
#'@param phyloseq   phyloseq object
#'@param file       output file
#'@param measures   list of diversity indices
#'
#'@return ggplot2 object
#'@export
#'
plot.richness <- function(phyloseq, 
                          file = "graphs/richness/richness.overview.pdf",
                          measures = c("Observed", "Shannon", 
                                       "Simpson","ACE", 
                                       "Chao1", "Fisher")) {
    
    # adjust theme options (see utils.R)
    ggtheme.alpha <- get.richnessTheme()

    # remove taxa with abundance values == 0
    phyloseq <- prune_taxa(taxa_sums(phyloseq) > 0, phyloseq)
    
    # create richness graph
    pp_rich <- plot_richness(phyloseq, 
                             color = "HoldingCondition", 
                             measures = measures)
    # disable legend
    pp_rich <- pp_rich + guides(color = FALSE)
    # pop out the original points
    pp_rich$layers <- pp_rich$layers[-1]
    
    # create the plot again
    pp_rich <- pp_rich + geom_point(size = 3, alpha = 0.5) +
               geom_boxplot(aes(y = value, fill = HoldingCondition), 
                            alpha = 0.1, outlier.size = 0) +
               ggtheme.alpha +
               xlab("Samples") +
               labs(fill = "Holding condition") +
               scale_color_manual(values = c("#6e8a3d", "#eeb422")) +
               scale_fill_manual(values = c("#6e8a3d", "#eeb422"))
    
    # save plot to file
    if(!is.null(file)) ggsave(file)
    return(pp_rich)
}
#' get.richness
#'
#' wrapper for phyloseq estimate_richness function 
#' 
#'@description wrapper to create data.frame with 
#'             all desired diversity indices per 
#'             sample.
#' 
#'@param phyloseq   phyloseq object
#'@param measures   list of diversity indices
#'
#'@return data.frame 
#'@export
#'
get.richness <- function(phyloseq,
                         measures = c("Observed", "Chao1", "ACE", 
                                      "Shannon", "Simpson", "Fisher")) {
    # calculate diversity indices
    data <- estimate_richness(phyloseq, measures = measures)
    # bring up more informations to data
    if (nrow(sample_data(phyloseq)) > 2) {
        rownames(data) <- as.character(sample_data(phyloseq)$SampleName)
    }
    data <- round(data, 2)
    return(data)
}
