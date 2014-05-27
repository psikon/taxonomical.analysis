# functions for evaluation of the different steps of the usearch pipeline

create_input_overview <- function(file, name) {
    # create overview of 16S input
    res <- readLines(file)
    res <- unlist(strsplit(res[grep('seqs', res)], ' '))
    res <- as.integer(res[grep('seqs',res)-1])
    output <- paste0("graphs/", name, '.jpeg')
    jpeg(output)
    barplot(res, 
            main = "Distribution of extracted \n16S sequences over samples",
            names.arg=c('60', '62', '64', '66', '68', '70', '74', '76', '78', '80', '82'),
            col=brewer.pal(11,"RdYlGn"), 
            xlab='Samples',
            ylab='Number of Sequences')
    dev.off()
    print(output)
}

plot_taxa_resolution <- function(phyloseq, 
                             filename, 
                             ranks = c('phylum', 'class', 'order', 'family', 'genus'),
                             absolute = FALSE,
                             sep = TRUE,
                             length_group1 = 5, 
                             length_group2 = 6,
                             title = "Resolution of OTUs") {
    # create data.fram with taxonomic ranks per sample
    res <- sapply(rownames(sample_data(phyloseq)), function(x) {
        # generate for every sample a vector with numbers of taxonomic levels
        get_phylo_levels(phyloseq, x, ranks, absolute)
    })
    colnames(res) <- sample_data(phyloseq)$SampleName
    # convert data for ggplot2
    data2 <- melt(res)
    # order the data like given ranks
    data2$X1 <- factor(data2$X1, levels = ranks)
    if (sep) {
        # generate environment seperation for data
        dest <- factor(c(rep("free living", length_group1 * length(ranks)),
                         rep("aqua culture", length_group2 * length(ranks))),
                       levels = c("free living", "aqua culture"))
        df <- cbind(data2, dest)
    } else {
        df <- data2
    }
    
    k <- ggplot(df, aes(x = X2, y = value, fill = X1)) + 
        # type of plot
        geom_bar(stat = "identity", position = "dodge") + 
        # names of the axis
        xlab("\nSample") + ylab("OTUs") +
        # adjust Y-Axis
        theme_bw() + scale_y_continuous(labels = comma, breaks = pretty_breaks(n = 15)) + 
        # adjust legend      
        guides(fill = guide_legend("Taxonomic Level")) +
        # rotate x axis labels      
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        # legend title
        ggtitle(title) 
    if (sep) {
        # seperate into two windows
        k <- k + facet_wrap(~dest, scales = "free_x") 
    } else {
        k <- k
    }
    # generate output path
    output <- paste0("graphs/", filename, ".pdf")
    # create output file
    pdf(output)
    plot(k)
    dev.off()
    # return location of output
    return(output)
}

plot_grouped_abundance <- function(phyloseq, 
                                   filename, 
                                   absolute = FALSE,
                                   sep = TRUE,
                                   length_group1 = 5, 
                                   length_group2 = 6,
                                   title = "Resolution of OTUs") {
}

