#' get.DBConnection.old 
#'
#' wrapper for the creation of the taxonomyReportDB connections 
#' of the old experiment run with metpipe
#'
#'@param metadataList result of function getMetaDataList()
#'                    
#'@keyword internal, deprecated
#'@return list of db con
#'
get.eukaryotaDB <- function(metadataList) {
    sample60.euk <- taxonomyReportDBConnect("data/eukaryotaDB/sample60.complete.db", metadataList[["60"]])
    sample62.euk <- taxonomyReportDBConnect("data/", metadataList[["62"]]) 
    sample66.euk <- taxonomyReportDBConnect("data/", metadataList[["66"]])
    sample68.euk <- taxonomyReportDBConnect("data/", metadataList[["68"]]) 
    sample70.euk <- taxonomyReportDBConnect("data/", metadataList[["70"]]) 
    sample72.euk <- taxonomyReportDBConnect("data/", metadataList[["72"]]) 
    sample75.euk <- taxonomyReportDBConnect("data/", metadataList[["74"]]) 
    sample76.euk <- taxonomyReportDBConnect("data/", metadataList[["76"]]) 
    sample78.euk <- taxonomyReportDBConnect("", metadataList[["78"]]) 
    sample80.euk <- taxonomyReportDBConnect("", metadataList[["80"]]) 
    sample82.euk <- taxonomyReportDBConnect("", metadataList[["82"]]) 
    return(list(sample60.old, sample62.old, sample66.old, sample68.old,sample70.old,
                sample76.old, sample78.old, sample80.old, sample82.old))
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
    return(list(sample60.bac, sample62.bac, sample64.bac, sample66.bac, sample68.bac, 
                sample70.bac, sample72.bac, sample74.bac, sample76.bac, sample78.bac, 
                sample80.bac, sample82.bac))
}
#' create.QueryIDFile
#'
#' extract from taxonomyReportDB object the query_ids with 
#' successfull annotation and create file of it
#'
#'@param path       path for output file
#'@param position   position of samples in db structure
#'                      
#'@keyword internal
#'@return list of files
#'
create.QueryIDFile <- function(path = "data/functional", 
                               position = NULL) {
    
    # get all taxonomyReportDB's
    db <- get.DBConnection.new(get.metadata.list())
    
    # reduce to desired dbs
    if (!is.null(position)) {
        db <- db[position]
    }
    
    files <- lapply(db, function(x) {
        # create filename for index file
        file <- paste0(path, 
                       "/", 
                       as.character(x$.metadata["SampleName"]), ".txt")
        # extract the query_ids
        query_ids <- db_query(conn(x),
                              "Select query_def from query", 1)
        # write query_ids to file line by line
        write.table(query_ids, 
                    file = file, 
                    quote = FALSE, 
                    row.names = F, 
                    col.names = F)
        # return filepath
        file
    })
    return(files)
}

#' create.FastaFromTaxonomyReportDB
#'
#' extract from fasta input file the reads with annotation 
#' result in database 
#'
#'@param path           path for output file
#'@param fasta.files    list of input fasta files
#'@param position       position of samples in db structure
#'                      
#'@keyword internal
#'
create.FastaFromTaxonomyReportDB <- function(path = "data/functional", 
                                             fasta.files, 
                                             position = NULL) {
    # create the index files
    index.files <- create.QueryIDFile(path = path, 
                                      position = position)
    for (i in seq_along(index.files)) {
        # resolve filename
        name <- strsplit(strsplit(index.files[[i]], 
                                  split = "/")[[1]][3], '.', 1)[[1]][1]
        message(paste0("create fasta file for: ", name))
        # construct output path
        output <- paste0(path, "/", name, ".fasta")
        # generate commandline string
        cmd <- paste0("python src/fastaExtractor.py --header ", 
                      index.files[[i]], " --fasta ", 
                      fasta.files[i], " --output ",
                      output)
        # run src/fastaExtractor.py
        system(cmd, wait = T)
    }
}


