
#' plot.mostAbundant.habitat
#'
#' create a plot with the most abundant species on the 
#' x axis seperated into the two habitats
#'             
#'@param phyloseq   phyloseq object
#'@param level      phylogenetic level for plotting
#'@param threshold  threshold for most abundant
#'@param title      title of plot
#'@param file       name of file
#'
#'@return ggplot2 plot
#'@export
plot.mostAbundant.habitat <- function(phyloseq, 
                                      level = "order", 
                                      threshold = 0.01,
                                      title = "Most abundant taxa per Holding Confition", 
                                      file = NULL) {
    # reduce phyloseq object to level
    phyloseq <- tax_glom(phyloseq, level)
    
    # select only the most abundance taxa by threshold
    f1  <- filterfun_sample(topp(threshold))
    wh1 <- genefilter_sample(phyloseq, f1, A = 1)
    phyloseq <- prune_taxa(wh1, phyloseq)
    
    # clean the taxa levels from underscore syntax
    phyloseq <- rm.underscore(phyloseq)
    
    # create a string for colorize the plot
    geom_str <- paste0("geom_bar(aes(color = ",level,", 
                       fill = ",level,"), stat = 'identity', 
                       position = 'stack')") 
    
    # draw plot
    p <- plot_bar(phyloseq, level, fill = level, 
                  facet_grid = ~HoldingCondition) + 
         eval(parse(text = geom_str)) + 
         xlab("Taxa") + 
         ylab("Abundance") + 
         ggtitle(title) +     
         scale_y_continuous(labels = comma, 
                            breaks = pretty_breaks(n = 10)) 
    # save plot in file
    if(!is.null(file)) ggsave(file)
    
    return(p)
}

#' plot.mostAbundant.sample
#'
#' create a plot with the samples on the x axis and the 
#' stacked abundance of the taxa on the y axis seperated 
#' under the two conditions
#' 
#'@param phyloseq   phyloseq object
#'@param level      phylogenetic level for plotting
#'@param threshold  threshold for most abundant
#'@param title      title of plot
#'@param file       name of file
#'
#'@return ggplot2 plot
#'@export
plot.mostAbundant.sample <- function(phyloseq, 
                                     level = "order", 
                                     threshold = 0.01,
                                     title = "Most abundant taxa per Sample", 
                                     file = NULL) {
    # reduce phyloseq object to level
    phyloseq <- tax_glom(phyloseq, level) 
    
    # select only the most abundance taxa by threshold
    f1  <- filterfun_sample(topp(threshold))
    wh1 <- genefilter_sample(phyloseq, f1, A = 1)
    phyloseq <-  prune_taxa(wh1, phyloseq)
    
    # clean the taxa levels from underscore syntax
    phyloseq <- rm.underscore(phyloseq)
    
    # create a string for colorize the plot
    geom_str <- paste0("geom_bar(aes(color = ",level,", fill = ",level,
                       "), stat = 'identity', position = 'stack')") 
    # draw plot
    p <- plot_bar(phyloseq, fill = level) + 
         facet_wrap(~ HoldingCondition, scales="free_x") +
         eval(parse(text = geom_str)) + 
         xlab("Samples") + 
         ylab("Abundance") + 
         ggtitle(title) +     
         scale_y_continuous(labels = comma, 
                            breaks = pretty_breaks(n = 10)) 
    # save to file
    if(!is.null(file)) ggsave(file)
    
    return(p)
}

# # overwrites the standard phyloseq function plot_bar() to fullfit 
# # custom needs of some graphic issues
# plot_bar <- function(phyloseq, 
#                      file = NULL,
#                      level, 
#                      sep = "HoldingCondition",
#                      title = "Overview") {
#     
#     # remove the underscores from tax_table syntax
#     phyloseq <- rm.underscore(phyloseq)
#     # create a geom string for dynamically selection of levels
#     geom <- paste0("geom_bar(aes(color = ", level ,", fill = ", level,
#                    "), stat = 'identity', position = 'stack')")
#     # create a facet_wrap string for dynamically selection of seperators
#     if (!is.null(sep)) {
#         facet <-  paste0("facet_wrap(~", sep, ", scales= 'free_x')")
#     } else {
#         facet <- ""
#     }
#     # create the barplot
#     p <- phyloseq::plot_bar(tax_glom(phyloseq, level), fill = level) 
#     # change some ggplot2 parameters
#     p <- p + eval(parse(text = geom)) + ggtitle(title) + theme_bw() + 
#         scale_y_continuous(labels = comma, breaks = pretty_breaks(n = 15)) +
#         xlab("Abundance") + ylab("Sample") + eval(parse(text = ""))
#     
#     if(!is.null(file)) ggsave(file)  
#     return(p)
# }

#' plot.taxonomy.graph
#'
#' create a graph showing the taxonomy of the annotated items in the
#' phyloseq object
#' 
#'@param phyloseq   phyloseq object
#'@param file       name of file
#'@param level      phylogenetic level for plotting
#'@param type       type of the phylogentic graph
#'
#'@export
plot.taxonomy.graph <- function(phyloseq, 
                                file = NULL, 
                                level = "order", 
                                type = "phylogram") {
    # convert the taxa table to data.frame in clear syntax
    taxa <- as.data.frame(tax_table(rm.underscore(phyloseq)))
    # seperate the desired level from data.frame
    df <- taxa[, 2:which(colnames(taxa) == level)]
    # change items with no clear classification at this level to 
    # a well defined label
    no.class <- apply(df, 1, function(x){ any("not classified" %in% x)}) == FALSE
    df <- df[-no.class, ]
    # build a newick string of the taxonomies
    tree <- read.tree(text = df2newick(df, innerlabel = T))
    # open file for plotting
    if(!is.null(file)) pdf(file)
    # plot the tree
    plot.phylo(tree, 
               type = type, 
               cex = 0.8, 
               root.edge = T,
               node.depth = 1,
               label.offset = 0.5)
    # change nodelabels
    nodelabels(tree$node.label, 
               adj = c(1, 1.5), 
               frame = "none", 
               cex = 0.6)
    # close the file
    if(!is.null(file)) dev.off()
}
