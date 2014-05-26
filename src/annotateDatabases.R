## Big DataSets
source(file="src/metadata.R")

## Big DataSets
generate.TaxonomyReport(blast_db_path = "data/blastDB/sample64.new.bac.db",
                        metadata = metadata64,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample64.new.bac.db",
                        bitscore_tolerance = 0.90)
 
generate.TaxonomyReport(blast_db_path = "data/blastDB/sample74.new.bac.db",
                        metadata = metadata74,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample74.new.bac.db",
                        bitscore_tolerance = 0.90)

## MiSeq-Rest

# free living # 

generate.TaxonomyReport(blast_db_path = "data/blastDB/sample60.new.blast.db",
                         metadata = metadata60,
                         taxon_db_path = "data/taxonomyDB/sample60.new.bac.db",
                         coverage_threshold = 0.30,
                         min_match = 50,
                         bitscore_tolerance = 0.90)

generate.TaxonomyReport(blast_db_path = "data/blastDB/s",
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

generate.TaxonomyReport(blast_db_path = "data/blastDB/sample68.new.bac.db",
                        metadata = metadata68,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample68.new.bac.db",
                        bitscore_tolerance = 0.90)

################################################################################

# aquaculture # 

generate.TaxonomyReport(blast_db_path = "data/blastDB/sample70.new.bac.db",
                        metadata = metadata70,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample70.new.bac.db",
                        bitscore_tolerance = 0.90)

generate.TaxonomyReport(blast_db_path = "data/blastDB/blast76_bac.db",
                        metadata = metadata76,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample76.new.bac.db",
                        bitscore_tolerance = 0.90)

generate.TaxonomyReport(blast_db_path = "data/blastDB/blast78_bac.db",
                        metadata = metadata78,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample78_bac.db",
                        bitscore_tolerance = 0.90)

generate.TaxonomyReport(blast_db_path = "data/blastDB/sample80.new.bac.db",
                        metadata = metadata80,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample80.new.bac.db",
                        bitscore_tolerance = 0.90)

generate.TaxonomyReport(blast_db_path = "data/blastDB/sample82.new.bac.db",
                        metadata = metadata82,
                        coverage_threshold = 0.30,
                        min_match = 50,
                        taxon_db_path = "data/taxonomyDB/sample82.new.bac.db",
                        bitscore_tolerance = 0.90)
