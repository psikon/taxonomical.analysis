require(grid)
require(ggplot2)
source("src/initAnalysis.R")

# create a copy
shared <- bacteria
# select shared OTUs
core_p_all <- filter_taxa(shared, function(x) all(x > 0), prune=TRUE)
# function for most abundant number

# select only 10 most abundant ones
shared_all_10 <- prune_taxa(most_abundant_taxa(core_p_all, 10), core_p_all)
# show it
tax_table(shared_all_10)[order(taxa_sums(shared_all_10), decreasing=TRUE), ]

core <- bacteria

# Rare Curves

abu <- sort(taxa_sums(core), decreasing = TRUE)
df.abu <- data.frame(Sample = "All OTUs", 
                          Read.abundance = unname(core.abu),
                          Rank.abundance = seq_along(core.abu))
df.abu <- rbind(df.core.abu, list("All OTUs", 0, 
                                       nrow(df.core.abu) + 1))
# all core
abu.core <- sort(attr(core_otus(core, 4), "abundance"), 
                      decreasing = TRUE)
df.abu.core <- data.frame(Sample = "Shared OTUs", 
                        Read.abundance = unname(abu.core),
                        Rank.abundance = seq_along(abu.core))
df.abu.core <- rbind(df.abu.core, list("Shared OTUs", 0, 
                               nrow(df.abu.core) + 1))

# free core
free.abu.core <- sort(attr(core_otus(subset_samples(core, 
                            HoldingCondition == "free living")), 
                            "abundance"), decreasing = TRUE)
df.free.abu.core <- data.frame(Sample = "Shared OTUs\nfree-living", 
                               Read.abundance = unname(free.abu.core),
                              Rank.abundance = seq_along(free.abu.core))
df.free.abu.core <- rbind(df.free.abu.core, list("Shared OTUs\nfree-living", 
                              0, nrow(df.free.abu.core) + 1))
# mari core
mari.abu.core <- sort(attr(core_otus(subset_samples(core, 
                            HoldingCondition == "mariculture")), 
                           "abundance"), decreasing = TRUE)
df.mari.abu.core <- data.frame(Sample = "Shared OTUs\nmariculture", 
                               Read.abundance = unname(mari.abu.core),
                              Rank.abundance = seq_along(mari.abu.core))
df.mari.abu.core <- rbind(df.mari.abu.core, list("Shared OTUs\nmariculture", 
                            0, nrow(df.mari.abu.core) + 1))
###
df <- rbind(df.abu, df.abu.core, df.free.abu.core, df.mari.abu.core)
df$Ind <- 0

snm <- sample_data(core)$SampleName
x <- list()
for (nm in snm) {
    x <- c(x, list(sort(compactZero(taxa_sums(subset_samples(core, SampleName == nm))), decreasing = TRUE)))
}
n.otu <- sapply(x, length)
df.ind <- do.call("rbind", Map(function(x, n, i) data.frame(Sample = "Individual fish", Read.abundance = x,
                                                            Rank.abundance = seq_len(n), Ind = i),
                               x = x, n = n.otu, i = seq_along(n.otu)))
df <- rbind(df, df.ind)
df <- subset(df, Sample == "All OTUs" | 
                 Sample == "Individual fish" | 
                 Sample == "Shared OTUs")
df <- subset(df, Sample == "All OTUs" | 
                 Sample == "Individual fish")
pp <- ggplot(df[df$Ind == 0,], aes(x = Rank.abundance, y = Read.abundance, colour = Sample)) +
    geom_line(size = 2) +
    geom_line(aes(x = Rank.abundance, y = Read.abundance, group = Ind),
              data = df[df$Ind != 0, ], size = 1, color = "gray", alpha = 0.4) +
    scale_x_log10(limits = c(1, 500)) +
    scale_y_log10(limits = c(1, 3e6)) +
    xlab("OTU rank abundance (n)") +
    ylab("OTU read abundance (n)") +
    theme(legend.position = "top") +
    labs(colour = "")
pp  

# define function
rel_abundance <- function(x, n) (x/n)/(sum(x)/n)
# create a new copy without __
core <- rm.underscore(bacteria)

# select free living samples
free <- subset_samples(core, HoldingCondition == "free living")
# control
sum(sort(taxa_sums(free), decreasing = TRUE) > 0)
length(core_otus(free))

# select core otus in free
core_free <- prune_taxa(taxa_names(free) %in% names(core_otus(free)), free)
# count order
tax_count <- table(tax_table(core_free)[, "class"])
# reduce to order
core_free <- tax_glom(core_free, taxrank = "class")
core_free <- rm.taxa(core_free,
                    rownames(tax_table(core_free)[which(tax_table(core_free)[,"class"] == "not classified")]))
# count again
tax_count <- tax_count[as.vector(tax_table(core_free)[, "class"])]
# create data.frame
free.df <- data.frame(
    Condition = "free-living",
    Abundance = rel_abundance(taxa_sums(core_free), 
                              nsamples(core_free)),
    nOTU = tax_count,
    Class = as.vector(tax_table(core_free)[, "class"])
)
# control of
sum(taxa_sums(core_free))/sum(taxa_sums(free))

# select mariculture samples
mari <- subset_samples(core, HoldingCondition == "mariculture")
# control
sum(sort(taxa_sums(mari), decreasing = TRUE) > 0)
length(core_otus(mari))

# select core of mari
core_mari <- prune_taxa(taxa_names(mari) %in% names(core_otus(mari)), mari)
# count on order
tax_count <- table(tax_table(core_mari)[, "class"])
# reduce to order
core_mari <- tax_glom(core_mari, taxrank = "class")
core_mari <- rm.taxa(core_mari,
                    rownames(tax_table(core_mari)[which(tax_table(core_mari)[,"class"] == "not classified")]))
# count again
tax_count <- tax_count[as.vector(tax_table(core_mari)[, "class"])]
# create df
mari.df <- data.frame(
    Condition = "mariculture",
    Abundance = rel_abundance(taxa_sums(core_mari), 
                              nsamples(core_mari)),
    nOTU = tax_count,
    Class = as.vector(tax_table(core_mari)[, "class"])
)
# control
sum(taxa_sums(core_mari))/sum(taxa_sums(mari))

# select shared
core_all <- prune_taxa(taxa_names(core) %in% names(core_otus(core)), core)
# count on order
tax_count <- table(tax_table(core_all)[, "class"])
# reduce to order
core_all <- tax_glom(core_all, taxrank = "class")
core_all <- rm.taxa(core_all,
    rownames(tax_table(core_all)[which(tax_table(core_all)[,"class"] == "not classified")]))
# count again
tax_count <- tax_count[as.vector(tax_table(core_all)[, "class"])]
all.df <- data.frame(
    Condition = "shared",
    Abundance = rel_abundance(taxa_sums(core_all),
                              nsamples(core_all)),
    nOTU = tax_count,
    Class = as.vector(tax_table(core_all)[, "class"])
)
# control
sum(taxa_sums(core_all))/sum(taxa_sums(core))

perc <- all.df$Abundance
perc <- round(perc*100,2)
names(perc) <- all.df$level
perc

# create data.frame for plotting
df <- rbind(free.df, mari.df, all.df)
df <- df[df$Abundance > 0.01, ]
df$Abundance <- df$Abundance * 100
pl <- ggplot(df, aes(x = factor(1), y = Abundance, fill = Class)) + 
    geom_bar(stat = "identity") +
    scale_y_discrete(breaks = c(0,25,50,75)) +
    facet_grid(facets = . ~ Condition) + 
    coord_polar(theta = "y") +
    xlab('') + ylab('') + ggtheme_core
pl
ggsave("graphs/shared_otus.pdf", width = 12, height = 3.5, dpi = 150)

tiff("graphs/shared_otus.tiff", compression = "lzw", ,
     width = 10, height = 2.7, units = "in",res = 600)
pl
dev.off()

jpeg("graphs/shared_otus.jpg", width = 10, height = 2.7, 
     units = "in", res = 300)
pl
dev.off()

