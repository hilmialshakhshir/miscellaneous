
library(biomaRt)
#Create the snpList used for the query
snpList <- SNP_full_list$variantID

#Mart used to map SNPs to Ensembl Gene IDs
grch38.snp = useMart(biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")

#Mart used to map Ensembl Gene IDs to Gene name
grch38 = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

#Mapping SNPs to ENsembl GeneIDs
snp2gene <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id", "chr_name", "associated_gene"), 
                filters = "refsnp_id", 
                values = can_genes$`Ensembl ID`, 
                mart = grch38.snp)

#Mapping ENsembl Gene IDs to Gene names
genename <- getBM(attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "external_gene_name","start_position","end_position"),
                filters = "ensembl_gene_id_version", 
                values =  genes, 
                mart = grch38,
                verbose = F)

snpgenename <- merge(snp2gene,gene2name, by.x = "ensembl_gene_stable_id", by.y="ensembl_gene_id", all.x=T)


#Mapping SNPs to ENsembl GeneIDs
snp2gene <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id", "chr_name", "associated_gene"), 
                  filters = "ensembl_gene", 
                  values = can_genes$`Ensembl ID`, 
                  mart = grch38.snp)


#Alternative code
grch38.snp <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
#snps <- c("rs1071646", "rs10494366")
snp2gene <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id", "chr_name", "associated_gene"), 
                  filters = "snp_filter", 
                  values = altRatio_matrix$variantID, 
                  mart = grch38.snp)

#To obtain the SNPs that are associated with the geneID
gene.snps <- filter(snp2gene, snp2gene$ensembl_gene_stable_id %in% can_genes$`Ensembl ID`)
#To get the SNPs that are associated with the genes
can.gene.ratio <- filter(altRatio_matrix, altRatio_matrix$variantID %in% can.snps$refsnp_id)

#Count number of NAs per patient
na_count <-sapply(can.gene.ratio, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)

#Set all values equal to {0, 1} to NA
cg.ratio.filter <- can.gene.ratio
cg.ratio.filter[cg.ratio.filter == 1] <- NA
cg.ratio.filter[cg.ratio.filter == 0] <- NA

na_count_filter <-sapply(cg.ratio.filter, function(y) sum(length(which(is.na(y)))))
na_count_filter <- data.frame(na_count_filter)

sample.matrix <- Dose_matrix
for (j in 2:ncol(cg.ratio.filter)) {
  sample.matrix[1,j] <- mad(cg.ratio.filter[,j], center = 0.5, na.rm = T)
}
sample.matrix <- sample.matrix[1,2:38]
#Trial without changing {0,1} to NA
cg.ratio <- can.gene.ratio
na.count <- sapply(cg.ratio, function(y) sum(length(which(is.na(y)))))
na.count <- data.frame(na.count)
imb.tot <- Dose_matrix[1,]
for (j in 2:ncol(cg.ratio)) {
  imb.tot[1,j] <- mad(cg.ratio[,j], center = 0.5, na.rm = T)
}


#-------------------------------------------------------------------------
#Mapping Ensembl Gene IDs to Gene names

can_genes$geneidv <- NA
geneidv <- getBM(attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "external_gene_name","start_position","end_position"),
                    filters = "ensembl_gene_id", 
                    values = can_genes$`Ensembl ID`, 
                    mart = grch38,
                    verbose = F)

gx.can <- filter(gene_vst, gene_vst$geneid %in% geneidv$ensembl_gene_id_version)

#----------------------------------
#Get the ensembl gene name for the ID

grch38 = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
gene.name <- getBM(attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "external_gene_name","start_position","end_position"),
                 filters = "ensembl_gene_id_version", 
                 values = eqtls.cis.cvrt$gene, 
                 mart = grch38,
                 verbose = F)

gene.snps <- filter(snp2gene, snp2gene$ensembl_gene_stable_id %in% genename$ensembl_gene_id)
