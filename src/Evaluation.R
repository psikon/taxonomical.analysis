source("src/initAnalysis.R")
# colnames(mdf)[3] <- "Abundance"
fixInNamespace(psmelt, pos = "package:phyloseq")

#############################################################
######## Evaluation of different pipeline steps #############
#############################################################

# plot absolute distribution of taxonomical levels
plot.taxaResolution(phylo.new,
                    file = "graphs/evaluation/taxa_res.abs.pdf",
                    absolute = TRUE, sep = TRUE, 
                    length_group1 = 5, length_group2 = 7,
                    title = "Taxonomical Resolution per samples \n(abs)")

# plot realtive distribution of taxonomical levels
plot.taxaResolution(phylo.new,
                    file = "graphs/evaluation/taxa_res.perc.pdf",
                    absolute = FALSE, sep = TRUE, 
                    length_group1 = 5, length_group2 = 7,
                    title = "Taxonomical Resolution per samples \n(percent)")

# plot the absolute apperance of specified occurence groups
plot.groupedAbundance(phylo.new, 
                      file = "graphs/evaluation/abundance.abs.pdf",
                      absolute = TRUE, sep = TRUE,
                      length_group1 = 5, length_group2 = 7,
                      title = "Abundance in defined groups per sample\n(absolute)")

# plot the relative apperance of specified occurence groups
plot.groupedAbundance(phylo.new, 
                      file = "graphs/evaluation/abundance.perc.pdf",
                      absolute = FALSE, sep = TRUE,
                      length_group1 = 5, length_group2 = 7,
                      title = "Abundance in defined groups per sample\n(percent)")

# plot the number of database hits per sample
plot.DBCount(data = get.DBConnection.new(get.metadata.list()) , 
             names = c("sample 60", "sample62","sample 64", "sample 66","sample 68", "sample 70", 
                       "sample 72", "sample 74", "sample 76", "sample 78", 
                       "sample 80", "sample 82"), 
             file = "graphs/evaluation/database_counts.pdf",
             sep = TRUE, length_group1 = 5, length_group2 = 7,
             title = "Hits in taxonomyReportDB per sample")
