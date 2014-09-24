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
#' get.DBConnection.new
#'
#' wrapper for the creation of the taxonomyReportDB connections 
#' of the new experiment run with meta-pipeline
#'
#'@param metadataList result of function getMetaDataList()
#'                      
#'@keyword internal
#'@return list of db con
#'
get.DBConnection.new <- function(metadataList) {
    sample60.new <- taxonomyReportDBConnect("data/taxonomyDB/sample60.new.bac.db", metadataList[["60"]])
    sample62.new <- taxonomyReportDBConnect("data/taxonomyDB/sample62.new.bac.db", metadataList[["62"]])
    sample64.new <- taxonomyReportDBConnect("data/taxonomyDB/sample64.new.bac.db", metadataList[["64"]])
    sample66.new <- taxonomyReportDBConnect("data/taxonomyDB/sample66.new.bac.db", metadataList[["66"]]) 
    sample68.new <- taxonomyReportDBConnect("data/taxonomyDB/sample68.new.bac.db", metadataList[["68"]]) 
    sample70.new <- taxonomyReportDBConnect("data/taxonomyDB/sample70.new.bac.db", metadataList[["70"]]) 
    sample72.new <- taxonomyReportDBConnect("data/taxonomyDB/sample72.new.bac.db", metadataList[["72"]]) 
    sample74.new <- taxonomyReportDBConnect("data/taxonomyDB/sample74.new.bac.db", metadataList[["74"]]) 
    sample76.new <- taxonomyReportDBConnect("data/taxonomyDB/sample76.new.bac.db", metadataList[["76"]]) 
    sample78.new <- taxonomyReportDBConnect("data/taxonomyDB/sample78.new.bac.db", metadataList[["78"]]) 
    sample80.new <- taxonomyReportDBConnect("data/taxonomyDB/sample80.new.bac.db", metadataList[["80"]]) 
    sample82.new <- taxonomyReportDBConnect("data/taxonomyDB/sample82.new.bac.db", metadataList[["82"]])
    return(list(sample60.new, sample62.new, sample64.new, sample66.new, sample68.new, 
                sample70.new, sample72.new, sample74.new, sample76.new, sample78.new, 
                sample80.new, sample82.new))
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


