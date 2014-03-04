library(ShortRead)
library(plyr)
library(blastr)

metacv60 <- scan("../epinephelus/sample60/metacv_result/metpipe.res",
                 what = list(header = character(0), 
                             score = integer(0),
                             geneid = character(0),
                             KEGGid = character(0),
                             COGid = character(0),
                             taxid = character(0),
                             taxname = character(0)),
                 sep = "\t")

metacv60 <- as.data.frame(metacv60, stringsAsFactors=FALSE)
metacv60 <- metacv60[which(metacv60$score >= 20),]
metacv60 <- unique(metacv60)
fastq60 <- readFastq(dirPath="../epinephelus/sequences/",pattern="60_S1_L001_R1_001.fastq")

subset60 <- metacv60[1:50,]
idx <- match(subset60$header,strsplitN(as.vector(id(fastq60)),"\\s",1))
seqs <- fastq60[idx]
writeFasta(seqs,file="seqs.fasta")
seqs <- readFasta("seqs.fasta")
blast <- blastn(sread(seqs[1]),db="nt", outfmt="table")
options(blastr.blastdb.path)
compare <- data.frame(header=subset60$header,score=subset60$score,taxname=subset60$taxname)
compare

blast <- blastReport("metacv.test.xml")
