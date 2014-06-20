

plot_variance <- function(phyloseq, x.text = "log10(variance)", title = "Variance of OTUs" ) {
  hist(log10(apply(otu_table(phyloseq),1, var)), xlab = x.text, main = title)
}

deseq_differential_OTUs <- function(phyloseq, variance_threshold, alpha) {
    keepOTUs <- apply(otu_table(phyloseq), 1, var) > variance_threshold
    data <- prune_taxa(keepOTUs, phyloseq)
    #create cds
    cds <- suppressWarnings(phyloseq_to_DESeq(data, 
                                              "HoldingCondition", 
                                              sharingMode = "gene-est-only"))
    res <- nbinomTest(cds, "free living", "mariculture")
    sigtab <- res[which(res$padj < alpha),]
    sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(data)[sigtab$id, ], "matrix"))
    sigtab <- sigtab[order(sigtab$padj),]
    sigtab
}

plot_differences <- function(sigtab) {
    theme_set(theme_bw())
    scale_fill_discrete <- function(palname = "Set1", ...) {
        scale_fill_brewer(palette = palname, ...)
    }
    
    # Phylum order
    x <- tapply(sigtab$log2FoldChange, sigtab$phylum, function(x)
        max(x))
    x <- sort(x,TRUE)
    sigtab$phylum <- factor(as.character(sigtab$phylum), levels = names(x))
    
    # Genus order
    x <- tapply(sigtab$log2FoldChange, sigtab$genus, function(x)
        max(x))
    x <- sort(x,TRUE)
    sigtab$genus <- factor(as.character(sigtab$genus), levels = names(x))
    
    ggplot(sigtab, 
           aes(x = genus, y = log2FoldChange, color = phylum)) + 
        geom_point(size = 3) + 
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
}

deseq2_differntial_OTUs <- function(phyloseq, alpha) {
    deseq2 <- suppressWarnings(phyloseq_to_deseq2(phyloseq, ~HoldingCondition))
    deseq2 <- DESeq(deseq2, test = "Wald", fitType="local")
    res <- results(deseq2)
    sigtab <- res[which(res$padj < alpha),]
    sigtab = cbind(as(sigtab,"data.frame"), 
                   as(tax_table(phyloseq)[rownames(sigtab), ], "matrix")) 
    sigtab
}
deseq2 <- deseq2_differntial_OTUs(bakteria, 0.3)
plot_differences(deseq2)
dge <- phyloseq_to_edgeR(bakteria, "HoldingCondition")
dge
et = exactTest(dge)
et
tt = topTags(et, n = nrow(dge$table),adjust.method= "BH", sort.by = "PValue")
res = tt@.Data[[1]]
res
alpha = 0.003
sigtab = res[(res$FDR < alpha),]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(bakteria)[rownames(sigtab), 
                                                               ], "matrix"))
dim(sigtab)
head(sigtab)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$phylum = factor(as.character(sigtabgen$phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$genus = factor(as.character(sigtabgen$genus), levels = names(x))
ggplot(sigtabgen, aes(x = genus, y = logFC, color = phylum)) + geom_point(size = 3) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))

