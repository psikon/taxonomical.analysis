library(goseq)
library(org.Hs.eg.db)

# getConnections()
# rm(sample60.old,sample62.old,sample64.old,sample64.new,sample68.old,
#   sample68.new,sample70.new,sample70.old,sample72.new,sample74.new,
#   sample76.new,sample76.old,sample66.old,sample78.new,sample78.old,
#   sample80.new,sample80.old,sample82.new,sample82.old)
# 
# hits <- db_query(conn(sample60.new), "Select * from hit where query_id between 1 and 20")
# geneids <- db_query(conn(sample60.new), "Select gene_id from hit where query_id between 1 and 20")
# getgo(geneids, "hg19","ensGene")
# supportedGenomes()
