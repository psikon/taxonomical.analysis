library(ProjectTemplate)
load.project()
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
registerDoParallel(cores = 20)
# in zeile 46 hinzufügen
# colnames(mdf)[3] <- "Abundance"
fixInNamespace(psmelt, pos = "package:phyloseq")


metadataList <- get.metadata.list()
bac <- get.bacterialDB(metadataList)
euk <- get.eukaryotaDB(metadataList)
biom <- generate.biomFile(file = "data/eukaryotic.biom", euk)
biom <- generate.biomFile(file = "data/bacterial.biom", bac)


