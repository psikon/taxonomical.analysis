get_core_microbiom <- function(phyloseq) {
    otu_tbl <- otu_table(phyloseq)
    sample_tbl <- sample_data(phyloseq)
    tax_tbl <- tax_table(phyloseq)
    no_zero = apply(otu_tbl, 1, function(row) all(row !=0 ))
    otu_tbl <- otu_tbl[no_zero,]
    remain_taxa <- rownames(otu_tbl)
    tax_tbl <- tax_tbl[rownames(tax_tbl) %in%  remain_taxa]
    phyloseq <- phyloseq(otu_tbl,sample_tbl,tax_tbl)
    phyloseq
}

plot_core_venn <- function(phyloseq, filename) {
    # seperate the samples by environment
    free <- get_core_microbiom(get_free(bakteria))
    aqua <- get_core_microbiom(get_aqua(bakteria))
    output <- paste0(filename, ".pdf")
    pdf(output)
        grid.draw(venn.diagram(x = list(
                                        free = rownames(otu_table(free)),
                                        aqua = rownames(otu_table(aqua))),
                               filename = NULL, main = "Core Microbiom",
                               main.cex = 2, main.fontface = "bold",
                               col = c("red","blue"), fill = c("red","blue"),
                               lwd = 2, alpha = 0.3, cex = 1.5, 
                               cat.cex = 1.5, cat.fontface = "bold"))
    dev.off()
}

get_core_table <- function(phyloseq, file = NULL, col.names = F, row.names = F) {

    data <- as.data.frame(t(sapply(seq_along(1:nrow(tax_table(phyloseq))), function(x){
                c("tax id" = rownames(otu_table(phyloseq)[x]),
                  "scientific.name" = sub("_"," ",last_taxa(tax_table(phyloseq)[x])),
                  "tax.rank" = last_rank(tax_table(phyloseq)[x]),
                  "count"= as.vector(otu_table(phyloseq)[x]))
    })))
    colnames(data) <- c("tax.id", "scientific.name", "tax.rank", 
                        paste("counts", colnames(otu_table(phyloseq))))
    #nr #tax_id #tax-name #level #count
    
    if(!is.null(file)) {
        write.table(data,file, sep="\t",quote=F,row.names=row.names,col.names=col.names)
    }
    data 
}
