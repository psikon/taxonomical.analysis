#' plot.ordintion.OTUs
#'
#' wrapper for the plot_ordination function of phyloseq
#'
#'@description  create a plot containing the ordination of OTUs 
#'              at a specific taxonomical level. the otus will be 
#'              filtered first, using only OTUs with hits occuring 
#'              'num_samples' and then the 'num_best' of them will 
#'              be used for plotting.
#' 
#'@param phyloseq       phyloseq object
#'@param file           output file
#'@param method         ordination method
#'@param distance       ordination distance
#'@param hits           min number of hits in a sample
#'@param num_samples    number of samples, where the hit must occur
#'@param num_best       number, how many of the best will be plotted
#'@param level          taxonomical level
#'@param title          plot title
#'@param facet          seperate the plots
#'@param sep            number of seperations
#'
#'@return ggplot2 object
#'@export
#'
plot.ordination.OTUs <- function(phyloseq, 
                                 file = "graphs/ordination/ordinate.otus.pdf",
                                 method = "NMDS", 
                                 distance = "bray", 
                                 hits = 1, 
                                 num_samples = 0.5, 
                                 num_best = 5, 
                                 level = "phylum",
                                 title = "Ordination of Taxa", 
                                 facet = T, 
                                 sep = 5) {
    
    # Remove OTUs that do not show appear more than x hits in 
    # more than num_samples the samples
    preprocess <- genefilter_sample(phyloseq, 
                                   filterfun_sample(function(x) x >= hits), 
                                   A = num_samples * nsamples(phyloseq))
    phyloseq <- prune_taxa(preprocess, phyloseq)
    
    # Transform to even sampling depth
    phyloseq <- transform_sample_counts(phyloseq, 
                                       function(x) 1e+06 * x/sum(x))
    
    # Keep only the most abundant num_best phyla
    level.sum <- tapply(taxa_sums(phyloseq), 
                       tax_table(phyloseq)[, level], 
                       sum,
                       na.rm = TRUE)
    toplevel <- compactNA(names(sort(level.sum, TRUE))[1:num_best])
    phyloseq <- prune_taxa((tax_table(phyloseq)[, level] %in% toplevel), 
                          phyloseq)
    # calculate ordination
    phyloseq.ord <- ordinate(phyloseq, 
                             level, 
                             method = method, 
                             distance = distance)
    # remove underscore symtax in tax_table
    phyloseq <- rm.underscore(phyloseq)
   
    # setup the theme
    theme_set(theme_bw())
    # create the plot
    p = plot_ordination(phyloseq, 
                        phyloseq.ord, 
                        type = "taxa", 
                        color = level, 
                        title = title)
    if (facet) {
        p <- p + eval(parse(text = paste0("facet_wrap(~", 
                                         as.name(level), ",", 
                                         sep, ")"))) 
    }
    # save plot to file
    if(!is.null(file)) ggsave(file)
    return(p)
}
#' plot.ordintion.samples
#'
#' wrapper for the plot_ordination function of phyloseq
#'
#'@description create a plot containing the ordination of 
#'             samples at a specific taxonomical level. the 
#'             otus of the samples will be filtered first, 
#'             using only OTUs with hits occuring 'num_samples' 
#'             and then the 'num_best' of them will be used
#' 
#'@param phyloseq       phyloseq object
#'@param file           output file
#'@param method         ordination method
#'@param distance       ordination distance
#'@param hits           min number of hits in a sample
#'@param num_samples    number of samples, where the hit must occur
#'@param num_best       number, how many of the best will be plotted
#'@param level          taxonomical level
#'@param title          plot title
#'
#'@return ggplot2 object
#'@export
#'
plot.ordination.samples <- function(phyloseq, 
                                    file = "graphs/ordination/ordinate.samples.pdf",
                                    method ="NMDS", 
                                    distance="bray", 
                                    hits = 5, 
                                    num_samples = 0.5, 
                                    num_best = 5,
                                    level, 
                                    title = "Ordination of Samples") {
    
    # Remove OTUs that do not show appear more than x hits 
    # in more than num_samples the samples
    preprocess <- genefilter_sample(phyloseq, 
                                   filterfun_sample(function(x) x >= hits), 
                                   A = num_samples * nsamples(phyloseq))
    phyloseq <- prune_taxa(preprocess, phyloseq)
    
    # Transform to even sampling depth
    phyloseq <- transform_sample_counts(phyloseq, 
                                       function(x) 1e+06 * x/sum(x))
    # Keep only the most abundant num_best phyla
    level.sum <- tapply(taxa_sums(phyloseq), 
                       tax_table(phyloseq)[, level], 
                       sum, 
                       na.rm = TRUE)
    toplevel <- compactNA(names(sort(level.sum, TRUE))[1:num_best])
    phyloseq <- prune_taxa((tax_table(phyloseq)[, level] %in% toplevel), phyloseq)
    phyloseq <- tax_glom(phyloseq, level)  
   
    # setup the theme
    theme_set(theme_bw())
    
    # claculate ordination
    phyloseq.ord <- ordinate(phyloseq, 
                             method = "MDS", 
                             distance = "bray")
    # create the plot
    p <- plot_ordination(phyloseq, 
                         phyloseq.ord, 
                         type = "samples", 
                         color = "SampleName", 
                         shape = "HoldingCondition")
    p <- p + geom_polygon(aes(fill = SampleName)) + geom_point(size = 3) + ggtitle(title)
    # save plot to file
    if(!is.null(file)) ggsave(file)
    return(p)
}
