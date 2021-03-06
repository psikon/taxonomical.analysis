# Evaluation

# Resolution of Taxa 
plot_taxa_resolution(bakteria, 
                     file = "graphs/evaluation/taxa_resolution.absolute.pdf",
                     absolute = T,
                     sep = T,
                     length_group1 = 3, 
                     length_group2 = 7,
                     title = "Resolution of OTUs (absolute)")

plot_taxa_resolution(bakteria, 
                     file = "graphs/evaluation/taxa_resolution.percent.pdf",
                     absolute = F,
                     sep = T,
                     length_group1 = 3, 
                     length_group2 = 7,
                     title = "Resolution of OTUs (in percent)")

# Abundance Table

plot_grouped_abundance(bakteria,
                       file = "graphs/evaluation/grouped_abundance.absolute.pdf",
                       sep = T,
                       length_group1 = 3,
                       length_group2 = 7, 
                       absolute = T,
                       title = "Grouped abundance of OTUs (absolute)")

plot_grouped_abundance(bakteria,
                       file = "graphs/evaluation/grouped_abundance.percent.pdf",
                       sep = T,
                       length_group1 = 3,
                       length_group2 = 7, 
                       absolute = F,
                       title = "Grouped abundance of OTUs (in percent)")

# Evaluation of TaxonomyReportDBs

plot_database_count(data = c(sample60.old,sample60.new,sample62.old,sample64.new,
                             sample66.old,sample68.old,sample68.new,sample70.old,
                             sample72.new,sample74.new,sample76.old,sample76.new,
                             sample78.old,sample80.old,sample80.new,sample82.old,
                             sample82.new),
                    names = c("sample60 metpipe", "sample60 meta-pipeline", "sample62 metpipe", 
                              "sample64 meta-pipeline","sample66 metpipe", "sample68 metpipe", 
                              "sample68 meta-pipeline", "sample70 metpipe", "sample72 meta-pipeline", 
                              "sample74 meta-pipeline", "sample76 metpipe", "sample76 meta-pipeline",
                              "sample78 metpipe", "sample80 metpipe", "sample80 meta-pipeline", 
                              "sample82 metpipe", "sample82 meta-pipeline"),
                    filename = "TaxonomyReportDB-distribution",
                    sep = T,
                    length_group1 = 7,
                    length_group2 = 10)