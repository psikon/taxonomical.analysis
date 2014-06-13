plot_richness_overview <- function(phyloseq, 
                                   file = "graphs/richness/richness.overview.pdf",
                                   measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson")) {
   
    p = plot_richness(phyloseq, 
                      color= "samples", 
                      measures = measures,
                      title = "Diversity Overview")
    p = p + geom_point(size = 2, alpha = 0.7) +
        guides(fill = guide_legend("Samples")) + xlab("\nSample")
    
    pdf(file)
        plot(p) 
    dev.off()
}

plot_richness_shannon <- function(phyloseq, 
                                  file = "graphs/richness/richness.shannon.pdf",
                                  measures = "Shannon") {
    p = plot_richness(phyloseq, 
                      color= "samples", 
                      measures = measures,
                      title = "Diversity Shannon index")
    p = p + geom_point(size = 2, alpha = 0.7) +
        guides(fill = guide_legend("Samples")) + xlab("\nSample")
    
    pdf(file)
    plot(p) 
    dev.off()
}

plot_richness_simpson <- function(phyloseq, 
                                  file = "graphs/richness/richness.simpson.pdf",
                                  measures = "Simpson") {
    p = plot_richness(phyloseq, 
                      color= "samples", 
                      measures = measures,
                      title = "Diversity Simpson index")
    p = p + geom_point(size = 2, alpha = 0.7) +
        guides(fill = guide_legend("Samples")) + xlab("\nSample")
    
    pdf(file)
    plot(p) 
    dev.off()
}

plot_richness_Chao <- function(phyloseq, 
                                  file = "graphs/richness/richness.chao.pdf",
                                  measures = "Chao1") {
    p = plot_richness(phyloseq, 
                      color= "samples", 
                      measures = measures,
                      title = "Diversity Chao1 index")
    p = p + geom_point(size = 2, alpha = 0.7) +
        guides(fill = guide_legend("Samples")) + xlab("\nSample")
    
    pdf(file)
    plot(p) 
    dev.off()
}

plot_richness_ACE <- function(phyloseq, 
                               file = "graphs/richness/richness.ACE.pdf",
                               measures = "ACE") {
    p = plot_richness(phyloseq, 
                      color= "samples", 
                      measures = measures,
                      title = "Diversity ACE index")
    p = p + geom_point(size = 2, alpha = 0.7) +
        guides(fill = guide_legend("Samples")) + xlab("\nSample")
    
    pdf(file)
    plot(p) 
    dev.off()
}

get_richness <- function(phyloseq,
                         measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson")) {
    data <- estimate_richness(phyloseq, measures = measures)
    if (nrow(sample_data(phyloseq))>2) rownames(data) <- as.character(sample_data(phyloseq)$SampleName)
    data <- round(data,2)
    data
}
#estimate_richness(bak,measures=c("Observed","Chao1","ACE","Shannon","Simpson"))