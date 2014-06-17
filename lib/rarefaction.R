# # rarified <- rarefy_even_depth(bakteria)
# # 
#rarefy(otu_table(bakteria), sample=3, se=T, MARGIN =1)
#  
# 
#  otus <- otu_table(bakteria)
#  otus <- as.data.frame(otus)
# quantile(rowSums(otus))
#  # rarecurve(otus, step = 1, sample, xlab = "Sample Size", ylab = "Species",label = TRUE)
#  raremax <- min(rowSums(otus))
#  Srare <- rarefy(otus, raremax)
#  pdf("graphs/rare.pdf")
#  rarecurve(otus,step = 1, sample = raremax, col = "blue", cex = 0.6)
#  dev.off()
