extract_singletons_for_sample <- function(phyloseq, num_samples) {
    non_zero <- rownames(otu_table(get_core_microbiom(phyloseq)))
    contain_singletons <- phyloseq(otu_table(phyloseq)[-which(rownames(otu_table(phyloseq)) %in% non_zero)],
                                    sample_data(phyloseq),
                                    tax_table(phyloseq)[-which(rownames(tax_table(phyloseq)) %in% non_zero)])
    contain_singletons
}