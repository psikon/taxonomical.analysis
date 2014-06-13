ordination_preprocess <- function(phyloseq, 
                                  hits = 1, 
                                  num_samples = 0.5,
                                  level = 'class',
                                  num_best = 10) {
    # Remove OTUs that do not show appear more than x hits in more than num_samples the samples
    preprocess = genefilter_sample(phyloseq, 
                                   filterfun_sample(function(x) x >= hits), 
                                   A = num_samples * nsamples(phyloseq))
    phyloseq = prune_taxa(preprocess, phyloseq)
    # Transform to even sampling depth
    phyloseq = transform_sample_counts(phyloseq, function(x) 1e+06 * x/sum(x))
    #phyloseq = transform_sample_counts(phyloseq, threshrankfun(10))
    # Keep only the most abundant num_best phyla
    level.sum = tapply(taxa_sums(phyloseq), tax_table(phyloseq)[, level], sum, na.rm = TRUE)
    toplevel = compactNA(names(sort(level.sum, TRUE))[1:num_best])
    phyloseq = prune_taxa((tax_table(phyloseq)[, level] %in% toplevel), phyloseq)
    phyloseq
}

plot_ordination_OTUs <- function(phyloseq, file = "graphs/ordination/ordinate.otus.pdf",
                                 method = "NMDS", distance = "bray", 
                                 level, title = "Ordination of Taxa", 
                                 facet = T, sep = 5) {
    theme_set(theme_bw())
    phyloseq.ord <- ordinate(phyloseq, method=method, distance=distance)
    p = plot_ordination(phyloseq, phyloseq.ord, type = "taxa", color = level, title=title)
    if (facet) {
        p = p + eval(parse(text=paste0("facet_wrap(~",as.name(level),",",sep,")"))) 
    }
    pdf(file)
        plot(p)
    dev.off()
}
   
plot_ordination_Samples <- function(phyloseq, file = "graphs/ordination/ordinate.samples.pdf",
                                    method ="NMDS", distance="bray", 
                                    level, title = "Ordination of Samples") {
    phyloseq.ord <- ordinate(phyloseq, method = method, distance = distance)
    p = plot_ordination(phyloseq, phyloseq.ord, type = "samples", color = "SampleName")
    p = p + geom_polygon(aes(fill = SampleName)) + geom_point(size = 3) + ggtitle(title)
    
    pdf(file)
        plot(p)
    dev.off()
}
