# function for holding the metadata of the experiment for every sample
get.metadata.list <- function() {
    
    metadata60 <- list(
        SampleId = 60, SampleName = "ef_free_60", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition = "free living", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 300, 
        TotalNumReads = 448479 , ProkReads = 307267, EukReads = 64, VirReads = 0, 
        Artifacts = 141150, OTUs = 174
    )
    
    metadata62 <- list(
        SampleId = 62, SampleName = "ef_free_62", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition = "free living", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 330,
        TotalNumReads = 326143 , ProkReads = 19096, EukReads = 155, VirReads = 0, 
        Artifacts = 451670, OTUs = 174
    )

    metadata64 <- list(
        SampleId = 64, SampleName = "ef_free_64", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition = "free living", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 187,
        TotalNumReads = 326142 , ProkReads = 292536, EukReads = 188, VirReads = 0, 
        Artifacts =  33413, OTUs = 450
    )
    
    metadata66 <- list(
        SampleId = 66, SampleName = "ef_free_66", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition = "free living", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 330,
        TotalNumReads = 5260 , ProkReads = 1565, EukReads = 74, VirReads = 0, 
        Artifacts = 451670, OTUs = 450
    )
    
    metadata68 <- list(
        SampleId = 68, SampleName = "ef_free_68", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition = "free living", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 328,
        TotalNumReads = 590435 , ProkReads = 394503, EukReads = 90, VirReads = 0, 
        Artifacts =  195840, OTUs = 107
    )
    
    metadata70 <- list(
        SampleId = 70, SampleName = "ef_cage_70", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition= "mariculture", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 200,
        TotalNumReads = 830720 , ProkReads = 151707, EukReads = 60, VirReads = 1, 
        Artifacts = 678950, OTUs = 104
    )
    
    metadata72 <- list(
        SampleId = 72, SampleName = "ef_cage_72", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition= "mariculture", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 300,
        TotalNumReads = 621426 , ProkReads = 377240, EukReads = 45, VirReads = 0, 
        Artifacts = 244140, OTUs = 85
    )
    
    metadata74 <- list(
        SampleId = 74, SampleName = "ef_cage_74", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition= "mariculture", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 75,
        TotalNumReads = 86449 , ProkReads = 44107, EukReads = 208, VirReads = 1, 
        Artifacts = 42124, OTUs = 443
    )
    
    metadata76 <- list(
        SampleId = 76, SampleName = "ef_cage_76", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition= "mariculture", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 68,
        TotalNumReads = 810522 , ProkReads = 134490, EukReads = 51, VirReads = 0, 
        Artifacts = 675980, OTUs = 142
    )
    
    metadata78 <- list(
        SampleId = 78, SampleName = "ef_cage_78", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition= "mariculture", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 140,
        TotalNumReads = 950972 , ProkReads = 148478, EukReads = 64, VirReads = 0, 
        Artifacts = 802430, OTUs = 87
    )
    
    metadata80 <- list(
        SampleId = 80, SampleName = "ef_cage_80", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition= "mariculture", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 190,
        TotalNumReads = 745240 , ProkReads = 121627, EukReads = 62, VirReads = 0, 
        Artifacts = 623550, OTUs = 96
    )
    
    metadata82 <- list(
        SampleId = 82, SampleName = "ef_cage_82", Location = "Jakarta", 
        Host = "Epinephelus fuscoguttatus", HoldingCondition= "mariculture", 
        Description = "fecal sample", DNAQuality = "good", DNAQuantity = 138,
        TotalNumReads = 546183 , ProkReads = 94464, EukReads = 49, VirReads = 0, 
        Artifacts = 451670, OTUs = 106
    )
    metadata <- list(metadata60, metadata62, metadata64, metadata66, metadata68,
                     metadata70, metadata72, metadata74, metadata76, metadata78,
                     metadata80, metadata82)
    names(metadata) <- c("60", "62", "64", "66", "68", "70", "72", "74", "76", "78", "80", "82")
    return(metadata)
}
