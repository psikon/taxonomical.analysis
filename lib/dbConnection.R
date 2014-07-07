# wrapper for the creation of the taxonomyReportDB connections of the old experiment run with metpipe
get.DBConnection.old <- function(metadataList) {
    sample60.old <- taxonomyReportDBConnect("data/taxonomyDB/sample60_bac.db", metadataList[["60"]])
    sample62.old <- taxonomyReportDBConnect("data/taxonomyDB/sample62_bac.db", metadataList[["62"]]) 
    sample66.old <- taxonomyReportDBConnect("data/taxonomyDB/sample66_bac.db", metadataList[["66"]])
    sample68.old <- taxonomyReportDBConnect("data/taxonomyDB/sample68_bac.db", metadataList[["68"]]) 
    sample70.old <- taxonomyReportDBConnect("data/taxonomyDB/sample70_bac.db", metadataList[["70"]]) 
    sample76.old <- taxonomyReportDBConnect("data/taxonomyDB/sample76_bac.db", metadataList[["76"]]) 
    sample78.old <- taxonomyReportDBConnect("data/taxonomyDB/sample78_bac.db", metadataList[["78"]]) 
    sample80.old <- taxonomyReportDBConnect("data/taxonomyDB/sample80_bac.db", metadataList[["80"]]) 
    sample82.old <- taxonomyReportDBConnect("data/taxonomyDB/sample82_bac.db", metadataList[["82"]]) 
    return(list(sample60.old, sample62.old, sample66.old, sample68.old,sample70.old,
                sample76.old, sample78.old, sample80.old, sample82.old))
}
# wrapper for the creation of the taxonomyReportDB connections of the new 
# experiment run with meta-pipeline
get.DBConnection.new <- function(metadataList) {
    sample60.new <- taxonomyReportDBConnect("data/taxonomyDB/sample60.new.bac.db", metadataList[["60"]])
    #sample62.new <- taxonomyReportDBConnect("data/taxonomyDB/sample62.new.bac.db", metadataList[["62"]])
    sample64.new <- taxonomyReportDBConnect("data/taxonomyDB/sample64.new.bac.db", metadataList[["64"]])
    #sample66.new <- taxonomyReportDBConnect("data/taxonomyDB/sample66.new.bac.db", metadataList[["66"]]) 
    sample68.new <- taxonomyReportDBConnect("data/taxonomyDB/sample68.new.bac.db", metadataList[["68"]]) 
    sample70.new <- taxonomyReportDBConnect("data/taxonomyDB/sample70.new.bac.db", metadataList[["70"]]) 
    sample72.new <- taxonomyReportDBConnect("data/taxonomyDB/sample72.new.bac.db", metadataList[["72"]]) 
    sample74.new <- taxonomyReportDBConnect("data/taxonomyDB/sample74.new.bac.db", metadataList[["74"]]) 
    sample76.new <- taxonomyReportDBConnect("data/taxonomyDB/sample76.new.bac.db", metadataList[["76"]]) 
    sample78.new <- taxonomyReportDBConnect("data/taxonomyDB/sample78.new.bac.db", metadataList[["78"]]) 
    sample80.new <- taxonomyReportDBConnect("data/taxonomyDB/sample80.new.bac.db", metadataList[["80"]]) 
    sample82.new <- taxonomyReportDBConnect("data/taxonomyDB/sample82.new.bac.db", metadataList[["82"]])
    return(list(sample60.new, sample64.new, sample68.new, sample70.new, sample72.new,
                sample74.new, sample76.new, sample78.new, sample80.new, sample82.new))
}