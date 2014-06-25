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
                             file = NULL, 
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
                         rep("mariculture", length_group2 * length(ranks))),
                       levels = c("free living", "mariculture"))
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
    # write to file
    if (!is.null(file)) ggsave(file)
    return(k)
}

plot_grouped_abundance <- function(phyloseq, 
                                   file = NULL, 
                                   absolute = FALSE,
                                   sep = TRUE,
                                   length_group1 = 4, 
                                   length_group2 = 5,
                                   title = "Grouped Abundance") {
    
    # get a data.frame with grouped abundance values for every sample
    res <- sapply(colnames(otu_table(phyloseq)),function(x) {
        # convert otu_table to data.frame
        data <- as.data.frame(otu_table(phyloseq))
        #generate counts
        counts <- count(data[, x])
        if (absolute) {
            # for absolute values sum all findings for manuell groups
            df <- c("no hits" = counts$freq[which(counts$x == 0)],
                    "1<x<10" = sum(counts$freq[which(counts$x>0 & counts$x<10)]),
                    "10<x<100" = sum(counts$freq[which(counts$x>10 & counts$x<100)]),
                    "100<x<1000" = sum(counts$freq[which(counts$x>100 & counts$x<1000)]),
                    "x>1000" = sum(counts$freq[which(counts$x>1000)]))
        } else {
            # for absolute values calculate percentage for manuell groups
            df <- c("no hits" = counts$freq[which(counts$x == 0)]*100/sum(counts$freq),
                    "1<x<10" = sum(counts$freq[which(counts$x>0 & counts$x<10)])*100/sum(counts$freq),
                    "10<x<100" = sum(counts$freq[which(counts$x>10 & counts$x<100)])*100/sum(counts$freq),
                    "100<x<1000" = sum(counts$freq[which(counts$x>100 & counts$x<1000)])*100/sum(counts$freq),
                    "x>1000" = sum(counts$freq[which(counts$x>1000)])*100/sum(counts$freq))
        }
       
    })
    # change labels for graph
    colnames(res) <- sample_data(phyloseq)$SampleName
    # convert for ggplot2
    data2 <- melt(res)
    # change names for easier access
    colnames(data2) <- c("group", "name", "value")
    # order the factor for ordered plot
    data2$group <- factor(data2$group, levels = c("no hits", "1<x<10", "10<x<100",
                                                  "100<x<1000", "x>1000"))
    # generate seperation if value sep=TRUE
    if (sep) {
        # generate environment seperation for data
        dest <- factor(c(rep("free living", length_group1 * 5),
                         rep("mariculture", length_group2 * 5)),
                       levels = c("free living", "mariculture"))
        df <- cbind(data2, dest)
    } else {
        df <- data2
    }
    # generate stacked or dodged plot
    k <- ggplot(df, aes(x = name, y = value, fill = group))
         # type of plot
    if (absolute) {
        k <- k + geom_bar(stat = "identity", position = "dodge") 
    } else {
        k <- k + geom_bar(stat = "identity", position = "stack") 
    }
    # generate seperation
    if (sep) {
        # seperate into two windows
        k <- k + facet_wrap(~dest, scales = "free_x") 
    } else {
        k <- k
    }
    # adjust rest for plotting
    k <- k + xlab("\nSamples") + ylab("Number of occurences") +
        theme_bw() + scale_y_continuous(labels = comma, breaks = pretty_breaks(n = 15)) +  
        guides(fill = guide_legend("Groups")) +    
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(title)  
    if (!is.null(file)) ggsave(file)
    k
}

plot_database_count <- function(data, 
                                names, 
                                file = NULL,                                    
                                sep = TRUE,
                                length_group1 = 9, 
                                length_group2 = 12,
                                title = "Number of Taxonomies") {
    # get length of taxonomy table per sample
    res <- sapply(data,function(x){
        db_query(conn(x),"SELECT COUNT(*) from taxonomy",1)
    })
    # name the rows like the databases
    names(res) <- names
    # convert for ggplot2
    data2 <- melt(res)
    # insert names
    data2$name <- rownames(data2)
    if (sep) {
        # generate environment seperation for data
        dest <- factor(c(rep("free living", length_group1),
                         rep("mariculture", length_group2 )),
                       levels = c("free living", "mariculture"))
        df <- cbind(data2, dest)
    } else {
        df <- data2
    }
    # generate  plot
    k <- ggplot(df, aes(x = name, y = value, fill= name))
    # type of plot
    k <- k + geom_bar(stat = "identity") 
    # generate seperation
    if (sep) {
        # seperate into two windows
        k <- k + facet_wrap(~dest, scales = "free_x") 
    } else {
        k <- k
    }
    # adjust rest for plotting
    k <- k + xlab("\nTaxonomy databases") + ylab("Number of taxonomies") +
        theme_bw() + scale_y_continuous(labels = comma, breaks = pretty_breaks(n = 15)) +  
        guides(fill = FALSE) +    
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(title)  
    # generate plot file
    if (!is.null(file)) ggsave(file)
    k
}

