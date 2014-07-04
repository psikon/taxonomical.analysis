# create a plot containing the ordination of OTUs at a specific taxonomical level. 
# the otus will be first filtered, using only OTUs with hits occuring 'num_samples' 
# and then the 'num_best' of them will be used
plot.ordination.OTUs <- function(phyloseq, file = "graphs/ordination/ordinate.otus.pdf",
                                 method = "NMDS", distance = "bray", 
                                 hits = 1, num_samples = 0.5, num_best = 5, 
                                 level = "phylum", title = "Ordination of Taxa", facet = T, sep = 5) {
    
    # Remove OTUs that do not show appear more than x hits in more than num_samples the samples
    preprocess = genefilter_sample(phyloseq, 
                                   filterfun_sample(function(x) x >= hits), 
                                   A = num_samples * nsamples(phyloseq))
    phyloseq = prune_taxa(preprocess, phyloseq)
    # Transform to even sampling depth
    phyloseq = transform_sample_counts(phyloseq, function(x) 1e+06 * x/sum(x))
    # Keep only the most abundant num_best phyla
    level.sum = tapply(taxa_sums(phyloseq), tax_table(phyloseq)[, level], sum, na.rm = TRUE)
    toplevel = compactNA(names(sort(level.sum, TRUE))[1:num_best])
    phyloseq = prune_taxa((tax_table(phyloseq)[, level] %in% toplevel), phyloseq)
    
    phyloseq.ord <- ordinate(phyloseq, level, method = method, distance = distance)
    phyloseq <- rm.underscore(phyloseq)
    # create the plot
    theme_set(theme_bw())
    p = plot_ordination(phyloseq, phyloseq.ord, type = "taxa", color = level, title = title)
    if (facet) {
        p = p + eval(parse(text = paste0("facet_wrap(~", as.name(level), ",", sep, ")"))) 
    }
    if(!is.null(file)) ggsave(file)
    return(p)
}
   
plot.ordination.samples <- function(phyloseq, file = "graphs/ordination/ordinate.samples.pdf",
                                    method ="NMDS", distance="bray", 
                                    hits = 5, num_samples = 0.5, num_best = 5,
                                    level, title = "Ordination of Samples") {
    # Remove OTUs that do not show appear more than x hits in more than num_samples the samples
    preprocess = genefilter_sample(phyloseq, 
                                   filterfun_sample(function(x) x >= hits), 
                                   A = num_samples * nsamples(phyloseq))
    phyloseq = prune_taxa(preprocess, phyloseq)
    
    # Transform to even sampling depth
    phyloseq = transform_sample_counts(phyloseq, function(x) 1e+06 * x/sum(x))
    # Keep only the most abundant num_best phyla
    level.sum = tapply(taxa_sums(phyloseq), tax_table(phyloseq)[, level], sum, na.rm = TRUE)
    toplevel = compactNA(names(sort(level.sum, TRUE))[1:num_best])
    phyloseq = prune_taxa((tax_table(phyloseq)[, level] %in% toplevel), phyloseq)
    phyloseq <- tax_glom(phyloseq, level)  
    
    theme_set(theme_bw())
    
    phyloseq.ord <- ordinate(phyloseq, method = "MDS", distance = "bray")
    p = plot_ordination(phyloseq, phyloseq.ord, type = "samples", 
                        color = "SampleName", shape = "HoldingCondition")
    p = p + geom_polygon(aes(fill = SampleName)) + geom_point(size = 3) + ggtitle(title)
    
    if(!is.null(file)) ggsave(file)
    return(p)
}
