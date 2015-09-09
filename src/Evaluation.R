source("src/initAnalysis.R")
# colnames(mdf)[3] <- "Abundance"
fixInNamespace(psmelt, pos = "package:phyloseq")

#' plot.taxaResolution
#'
#' plot the occurences of taxonomical levels per sample 
#' 
#'@param phyloseq       phyloseq object
#'@param file           output file
#'@param ranks          taxonomic levels
#'@param absolute       plot absolute or percentage values
#'@param sep            seperate plot by 'HoldingCondition'
#'@param length_group1  number of free_living samples
#'@param length_group2  number of mariculture samples
#'@param title          plot title
#'
#'@return ggplot2 object
#'@export
#'
plot.taxaResolution <- function(phyloseq, 
                                file = NULL, 
                                ranks = c('phylum', 'class', 'order', 
                                          'family', 'genus'),
                                absolute = FALSE,
                                sep = TRUE,
                                length_group1 = 5, 
                                length_group2 = 6,
                                title = "Resolution of OTUs") {
    
    # create data.fram with taxonomic ranks per sample
    res <- sapply(rownames(sample_data(phyloseq)), 
                  function(x) {
                      # generate for every sample a vector with numbers 
                      # of taxonomic levels
                      get.phyloLevels(phyloseq, x, ranks, absolute)
                  })
    colnames(res) <- sample_data(phyloseq)$SampleName
    
    # convert data for ggplot2
    data2 <- melt(res)
    
    # order the data like given ranks
    data2$X1 <- factor(data2$X1, levels = ranks)
    
    if (sep) {
        # generate environment seperation for data
        dest <- factor(c(rep("free living", 
                             length_group1 * length(ranks)),
                         rep("mariculture", 
                             length_group2 * length(ranks))),
                       levels = c("free living", "mariculture"))
        df <- cbind(data2, dest)
    } else {
        df <- data2
    }
    
    # create the plot
    k <- ggplot(df, aes(x = X2, y = value, fill = X1)) + 
        # type of plot
        geom_bar(stat = "identity", position = "dodge") + 
        # names of the axis
        xlab("\nSample") + ylab("OTUs") +
        # adjust theme
        theme_bw() + 
        # adjust y axis
        scale_y_continuous(labels = comma, 
                           breaks = pretty_breaks(n = 15)) + 
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

#' plot.groupedAbundance
#'
#' plot distribution of abundance values grouped in
#' defined abundance groups
#' 
#'@param phyloseq       phyloseq object
#'@param file           output file
#'@param absolute       plot absolute or percentage values
#'@param sep            seperate plot by 'HoldingCondition'
#'@param length_group1  number of free_living samples
#'@param length_group2  number of mariculture samples
#'@param title          plot title
#'
#'@return ggplot2 object
#'@export
#'  
plot.groupedAbundance <- function(phyloseq, 
                                  file = NULL, 
                                  absolute = FALSE,
                                  sep = TRUE,
                                  length_group1 = 4, 
                                  length_group2 = 5,
                                  title = "Grouped Abundance") {
    
    # get a data.frame with grouped abundance values for every sample
    res <- sapply(colnames(otu_table(phyloseq)), 
                  function(x) {
                      # convert otu_table to data.frame
                      data <- as.data.frame(otu_table(phyloseq))
                      
                      #generate counts
                      counts <- count(data[, x])
                      
                      if (absolute) {
                          # for absolute values sum all findings for manuell groups
                          df <- c("no hits" = counts$freq[which(counts$x == 0)],
                                  "1<x<10" = sum(counts$freq[which(counts$x > 0 & counts$x < 10)]),
                                  "10<x<100" = sum(counts$freq[which(counts$x > 10 & counts$x < 100)]),
                                  "100<x<1000" = sum(counts$freq[which(counts$x > 100 & counts$x < 1000)]),
                                  "x>1000" = sum(counts$freq[which(counts$x > 1000)]))
                      } else {
                          # for absolute values calculate percentage for manuell groups
                          df <- c("no hits" = counts$freq[which(counts$x == 0)] * 100/sum(counts$freq),
                                  "1<x<10" = sum(counts$freq[which(counts$x > 0 & counts$x < 10)]) * 100/sum(counts$freq),
                                  "10<x<100" = sum(counts$freq[which(counts$x > 10 & counts$x < 100)]) * 100/sum(counts$freq),
                                  "100<x<1000" = sum(counts$freq[which(counts$x > 100 & counts$x < 1000)]) * 100/sum(counts$freq),
                                  "x>1000" = sum(counts$freq[which(counts$x > 1000)]) * 100/sum(counts$freq))
                      }
                      
                  })
    # change labels for graph
    colnames(res) <- sample_data(phyloseq)$SampleName
    # convert for ggplot2
    data2 <- melt(res)
    # change names for easier access
    colnames(data2) <- c("group", "name", "value")
    # order the factor for ordered plot
    data2$group <- factor(data2$group, 
                          levels = c("no hits", "1<x<10", "10<x<100",
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
        theme_bw() + 
        scale_y_continuous(labels = comma, 
                           breaks = pretty_breaks(n = 15)) +  
        guides(fill = guide_legend("Groups")) +    
        theme(axis.text.x = element_text(angle = 45, 
                                         hjust = 1)) +
        ggtitle(title) 
    # save in file and return plot
    if (!is.null(file)) ggsave(file)
    return(k)
}

# generate read overview graphics
require(reshape2)
require(plyr)
require(scales)



RAW <- unique(read.table("data/evaluation/reads-raw.txt", 
                         header = F, as.is = T)[1])/4
classify <- unique(read.table("data/evaluation/reads-merged.txt",
                              header = F, as.is = F)[1])/2
r1.paired <- unique(read.table("data/evaluation/reads-paired.quality.txt",
                               header = F, as.is = T)[1])/4
r1.single <- unique(read.table("data/evaluation/reads-single.quality.txt",
                               header = F, as.is = T)[1])/4

RAW
quality <- r1.paired + r1.single

bacteria <- subset_taxa(generate.phyloseq("data/bacterial.biom"), superkingdom == 'k__Bacteria')
bacteria <- rm.taxa(bacteria, 2)
eukaryotic <- subset_taxa(generate.phyloseq("data/eukaryotic.biom"), superkingdom == 'k__Eukaryota')
bacterial <- colSums(otu_table(bacteria))
eukarotical <- colSums(otu_table(eukaryotic))

mapped <- bacterial + eukarotical
functional <- c(45488, 243651, 403117, 55133, 53924, 78806, 
                51673, 162204, 48962, 47896, 51876, 37484) 

func <- functional
tax <- mapped - func
qual <- quality[-13, ] - tax
raw <- RAW[-13, ] - qual

data <- data.frame(raw, qual, tax, func)


colnames(data) <- c("RAW", "quality-processed", "taxonomical annotated", "functional annotated")
rownames(data) <- c("free1", "free2", "free3", "free4", "free5", 
                    "mari1", "mari2", "mari3", "mari4", "mari5", 
                    "mari6", "mari7")
data[,"SampleName"] <- rownames(data)

## reorder the columns
data <- data[, c("SampleName", "functional annotated",
                 "taxonomical annotated","quality-processed", "RAW")]

## rename columns
#names(data) <- c("SampleID", "classified reads",
#                 "quality-filtered reads", "raw reads")
## add marginal
data <- rbind(data, c(SampleName = "Mean", trunc(colMeans(data[, 2:ncol(data)]))))
data[, 2:ncol(data)] <- lapply(data[, 2:ncol(data)], as.integer)
## melt
mcount <- melt(data, id.vars="SampleName", variable.name="Reads", value.name="value")
colnames(mcount) <- c("SampleName", "Reads", "value")

p <- ggplot(mcount, aes(x = `SampleName`, y = value, fill = Reads)) +
    geom_bar(stat="identity") + 
    scale_y_continuous(labels = comma, breaks = pretty_breaks(n=10)) +
    scale_fill_grey(start = 0.33, end = 0.67) +
    ylab("Number of reads") +
    xlab("Samples") +
    #guides(fill = guide_legend(title = NULL)) +
    ggtheme_reads +
    geom_hline(yintercept = sum(subset(mcount, `SampleName` == "Mean", 
                                       "value", drop = TRUE)),
               linetype = "dashed", colour = "gray33", size = 1) +
    geom_hline(yintercept = subset(mcount, `SampleName` == "Mean" &
                                       Reads == "taxonomical annotated", 
                                   "value", drop = TRUE) + 
                            subset(mcount, `SampleName` == "Mean" &
                              Reads == "functional annotated", 
                          "value", drop = TRUE),
               linetype = "dashed", colour = "gray66", size = 1)
p
ggsave("graphs/evaluation/read_distribution.pdf", width=6.3, height=4.3)

tiff("graphs/evaluation/read_distribution.tiff", compression = "lzw", ,
     width = 6.3, height = 4.3, units = "in",res = 600)
    p
dev.off()

jpeg("graphs/evaluation/read_distribution.jpg", width = 6.3, 
     height = 4.3, units = "in", res = 300)
p
dev.off()
