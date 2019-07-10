library(dplyr)

#required data
geneloc <- which(eqtls.cis.cvrt$snps == "rs34055023")
genename <- which(gene.name$ensembl_gene_id_version == eqtls.cis.cvrt[45, 2])
genex <- which(gene_vst$geneid == )
genotype <- 
snpID <- 
snpPval <- 

p1 <- ggplot2:::qplot(genotype, as.numeric(genex),
                      xlab = snpID,
                      colour = gt, geom = c("boxplot", "jitter"),
                      ylab = "Gene expression (log2)") + theme(
                        panel.background = element_blank(),
                        legend.position = "none") +
  theme(text = element_text(size=14), axis.title.x = element_text(vjust = -2)) +
  theme(text = element_text(size=14), axis.title.y = element_text(vjust = -1)) +
  theme(plot.margin = unit(c(0.5,0.5,1,2), "cm")) +
  theme(axis.text.x=element_text(size=12, color = "black")) +
  theme(axis.text.y=element_text(size=12, color = "black")) +
  ylim(c(min(.genex) - 0.1, max(.genex))) +
  labs(title = " ") #labs(title = " \n \n ")
