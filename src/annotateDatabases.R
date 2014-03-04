## Big DataSets
source(file="src/metadata.R")

## Big DataSets
generate.TaxonomyReport(blast_db_path = "../sample64/analysis/blastn.db",
                        metadata = metadata64,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/sample64.db",
                        bitscore_tolerance = 0.90)
 
generate.TaxonomyReport(blast_db_path = "../sample74/analysis/blastn.db",
                        metadata = metadata74,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/sample74.db",
                        bitscore_tolerance = 0.90)

## MiSeq-Rest

# free living # 

generate.TaxonomyReport(blast_db_path = "data/blastDB/blast60_bac.db",
                         metadata = metadata60,
                         taxon_db_path = "data/taxonomyDB/sample60_bac.db",
                         coverage_threshold = 0.30,
                         min_match = 50,
                         bitscore_tolerance = 1, log = NULL)

generate.TaxonomyReport(blast_db_path = "data/blastDB/blast62_bac.db",
                        metadata = metadata62,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample62.db",
                        bitscore_tolerance = 0.90)

generate.TaxonomyReport(blast_db_path = "data/blastDB/blast66_bac.db",
                        metadata = metadata66,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample66_bac.db",
                        bitscore_tolerance = 0.90)

generate.TaxonomyReport(blast_db_path = "data/blastDB/blast68_bac.db",
                        metadata = metadata68,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample68_bac.db",
                        bitscore_tolerance = 0.90)

################################################################################

# aquaculture # 

generate.TaxonomyReport(blast_db_path = "data/blastDB/blast70_bac.db",
                        metadata = metadata70,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample70_bac.db",
                        bitscore_tolerance = 0.90)

generate.TaxonomyReport(blast_db_path = "data/blastDB/blast76_bac.db",
                        metadata = metadata76,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample76_bac.db",
                        bitscore_tolerance = 0.90)

generate.TaxonomyReport(blast_db_path = "data/blastDB/blast78_bac.db",
                        metadata = metadata78,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample78_bac.db",
                        bitscore_tolerance = 0.90)

generate.TaxonomyReport(blast_db_path = "data/blastDB/blast80_bac.db",
                        metadata = metadata80,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample80_bac.db",
                        bitscore_tolerance = 0.90)

generate.TaxonomyReport(blast_db_path = "data/blastDB/blast82_bac.db",
                        metadata = metadata82,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample82_bac.db",
                        bitscore_tolerance = 0.90)
