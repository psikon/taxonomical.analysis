#' get.eukaryotaDB 
#'
#' wrapper for the creation of the taxonomyReportDB connections 
#' of the eukaryotic database runs
#'
#'@param metadataList result of function getMetaDataList()
#'                    
#'@keyword internal, deprecated
#'@return list of db con
#'
get.eukaryotaDB <- function(metadataList) {
    sample60.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample60.tax.db", metadataList[["60"]])
    sample62.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample62.tax.db", metadataList[["62"]]) 
    sample64.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample64.tax.db", metadataList[["64"]]) 
    sample66.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample66.tax.db", metadataList[["66"]])
    sample68.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample68.tax.db", metadataList[["68"]]) 
    sample70.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample70.tax.db", metadataList[["70"]]) 
    sample72.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample72.tax.db", metadataList[["72"]]) 
    sample74.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample74.tax.db", metadataList[["74"]]) 
    sample76.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample76.tax.db", metadataList[["76"]]) 
    sample78.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample78.tax.db", metadataList[["78"]]) 
    sample80.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample80.tax.db", metadataList[["80"]]) 
    sample82.euk <- taxonomyReportDBConnect("data/eukaryotaDB/taxonomy/sample82.tax.db", metadataList[["82"]]) 
    return(list(sample60.euk, sample62.euk, sample64.euk, 
                sample66.euk, sample68.euk, sample70.euk, 
                sample72.euk, sample74.euk, sample76.euk, 
                sample78.euk, sample80.euk, sample82.euk))
}
#' get.bacterialDB
#'
#' wrapper for the creation of the taxonomyReportDB connections 
#' of the new experiment run with meta-pipeline
#'
#'@param metadataList result of function getMetaDataList()
#'                      
#'@keyword internal
#'@return list of db con
#'
get.bacterialDB <- function(metadataList) {
    sample60.bac <- taxonomyReportDBConnect("data/bacterialDB/sample60.new.bac.db", metadataList[["60"]])
    sample62.bac <- taxonomyReportDBConnect("data/bacterialDB/sample62.new.bac.db", metadataList[["62"]])
    sample64.bac <- taxonomyReportDBConnect("data/bacterialDB/sample64.new.bac.db", metadataList[["64"]])
    sample66.bac <- taxonomyReportDBConnect("data/bacterialDB/sample66.new.bac.db", metadataList[["66"]]) 
    sample68.bac <- taxonomyReportDBConnect("data/bacterialDB/sample68.new.bac.db", metadataList[["68"]]) 
    sample70.bac <- taxonomyReportDBConnect("data/bacterialDB/sample70.new.bac.db", metadataList[["70"]]) 
    sample72.bac <- taxonomyReportDBConnect("data/bacterialDB/sample72.new.bac.db", metadataList[["72"]]) 
    sample74.bac <- taxonomyReportDBConnect("data/bacterialDB/sample74.new.bac.db", metadataList[["74"]]) 
    sample76.bac <- taxonomyReportDBConnect("data/bacterialDB/sample76.new.bac.db", metadataList[["76"]]) 
    sample78.bac <- taxonomyReportDBConnect("data/bacterialDB/sample78.new.bac.db", metadataList[["78"]]) 
    sample80.bac <- taxonomyReportDBConnect("data/bacterialDB/sample80.new.bac.db", metadataList[["80"]]) 
    sample82.bac <- taxonomyReportDBConnect("data/bacterialDB/sample82.new.bac.db", metadataList[["82"]])
    return(list(sample60.bac, sample62.bac, sample64.bac, 
                sample66.bac, sample68.bac, sample70.bac, 
                sample72.bac, sample74.bac, sample76.bac, 
                sample78.bac, sample80.bac, sample82.bac))
}