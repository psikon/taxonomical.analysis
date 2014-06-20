get_core_microbiom <- function(phyloseq) {
    # split up the phyloseq object in the 3 sperate parts
    otu_tbl <- otu_table(phyloseq)
    sample_tbl <- sample_data(phyloseq)
    tax_tbl <- tax_table(phyloseq)
    # find all OTUs contained in every sample
    no_zero = apply(otu_tbl, 1, function(row) all(row !=0 ))
    # get desired subset the OTUs 
    otu_tbl <- otu_tbl[no_zero,]
    # extract tax_ids of the remaining taxa
    remain_taxa <- rownames(otu_tbl)
    # subset the desired taxa from tax_table to make both equal
    tax_tbl <- tax_tbl[rownames(tax_tbl) %in%  remain_taxa]
    # build a new phyloseq object
    phyloseq <- phyloseq(otu_tbl, sample_tbl, tax_tbl)
    phyloseq
}

plot_core_venn <- function(phyloseq, filename) {
    # seperate the samples by environment
    free <- get_core_microbiom(get_free(bakteria))
    aqua <- get_core_microbiom(get_aqua(bakteria))
    # generate output path
    output <- paste0(filename, ".pdf")
    pdf(output)
        # build the venn diagram
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
    # build a data frame with desired inormations
    data <- as.data.frame(t(sapply(seq_along(1:nrow(tax_table(phyloseq))), function(x){
                c("tax id" = rownames(otu_table(phyloseq)[x]),
                  "scientific.name" = sub("_"," ",last_taxa(tax_table(phyloseq)[x])),
                  "tax.rank" = last_rank(tax_table(phyloseq)[x]),
                  "count"= as.vector(otu_table(phyloseq)[x]))
    })))
    # adjust column names
    colnames(data) <- c("tax.id", "scientific.name", "tax.rank", 
                        paste("counts", colnames(otu_table(phyloseq))))
    # save the data frame in a tab separeted file
    if(!is.null(file)) {
        write.table(data,file, sep="\t",quote=F,row.names=row.names,col.names=col.names)
    }
    data 
}