plot_most_abundant_per_habitat <- function(phyloseq, level = "order", threshold = 0.01,
                                           title = "Most abundant taxa per Holding Confition", 
                                           file = NULL) {
    # reduce phyloseq object to level
    phyloseq <- tax_glom(phyloseq, level) 
    # select only the most abundance taxa by threshold
    f1  <- filterfun_sample(topp(threshold))
    wh1 <- genefilter_sample(phyloseq, f1, A = 1)
    phyloseq <- prune_taxa(wh1, phyloseq)
    # clean the taxa levels from underscore syntax
    phyloseq <- remove_Underscore(phyloseq)
    # create a string for colorize the plot
    geom_str <- paste0("geom_bar(aes(color = ",level,", fill = ",level,"), stat = 'identity', position = 'stack')") 
    # draw plot
    p <- plot_bar(phyloseq, level, fill = level, facet_grid = ~HoldingCondition) + 
         eval(parse(text = geom_str)) + xlab("Taxa") + ylab("Abundance") + ggtitle(title) +     
         scale_y_continuous(labels = comma, breaks = pretty_breaks(n = 10)) 
    if(!is.null(file)) {
        ggsave(file)
    }
    return(p)
}

plot_most_abundant_per_sample <- function(phyloseq, level = "order", threshold = 0.01,
                                          title = "Most abundant taxa per Sample", 
                                          file = NULL) {
    # reduce phyloseq object to level
    phyloseq <- tax_glom(phyloseq, level) 
    # select only the most abundance taxa by threshold
    f1  <- filterfun_sample(topp(threshold))
    wh1 <- genefilter_sample(phyloseq, f1, A = 1)
    phyloseq <-  prune_taxa(wh1, phyloseq)
    # clean the taxa levels from underscore syntax
    phyloseq <- remove_Underscore(phyloseq)
    # create a string for colorize the plot
    geom_str <- paste0("geom_bar(aes(color = ",level,", fill = ",level,
                       "), stat = 'identity', position = 'stack')") 
    # draw plot
    p <- plot_bar(phyloseq, fill = level) + facet_wrap(~ HoldingCondition, scales="free_x") +
         eval(parse(text = geom_str)) + xlab("Samples") + ylab("Abundance") + ggtitle(title) +     
         scale_y_continuous(labels = comma, breaks = pretty_breaks(n = 10)) 
    # save to file
    if(!is.null(file)) ggsave(file)
    return(p)
}