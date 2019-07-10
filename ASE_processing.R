
#----------------------------------------------------------------------
# This will modify all the Sample ASE files and add the genotype, and allele dose
  
  for (j in 1:length(samples)) {
    ASE <- samples[[j]]
    ASE$gene1 <- NA
    ASE$gene2 <- NA
    ASE$Genotype <- NA
    ASE$alleledose <- 0
    for (i in 1:nrow(ASE)) {
      if (ASE[i, 6] == 0) {
        ASE[i, 14:15] <- ASE[i, 5]
      }
      else if (ASE[i, 7] == 0) {
        ASE[i, 14:15] <- ASE[i, 4]
      }
      else
        ASE[i, 14:15] <- ASE[i, 4:5]
    }
    samples[[j]] <- ASE
  }

#----------------------------------------------------------------------
#This is the loop that calculates the allele dose per SNP  
  for (j in 1:length(samples)) {
    ASE <- samples[[j]]
    ASE$Genotype <- paste(ASE$gene1, ASE$gene2)
    for (i in 1:nrow(ASE)) {
      if (ASE[i, 14] == ASE[i, 5] & ASE[i, 15] == ASE[i, 5]) {
        ASE[i, 17] = ASE[i, 17]+2
      }
      else if (ASE[i, 15] == ASE[i, 5]) {
        ASE[i, 17] = ASE[i, 17]+1
      }
      else if (ASE[i, 14] == ASE[i, 5]) {
        ASE[i, 17] = ASE[i, 17]+1
      }
      else ASE[i, 17] = ASE[i, 17]+0
    }
    samples[[j]] <- ASE
  }


#This is to calculate the refRatio/altRatio
for (j in 1:length(samples)) {
  ASE <- samples[[j]]
  ASE$refRatio <- 0
  ASE$altRatio <- 0
  ASE$refRatio <- ASE$refCount/ASE$totalCount
  ASE$altRatio <- ASE$altCount/ASE$totalCount
  samples[[j]] <- ASE
}

#----------------------------------------------------------------------
# This part will create the master list of all the SNPs amongst all samples, and then will create the dose matrix.
#Initialize the master list
snp2 <- samples[[1]]
snp3 <- samples[[2]]

snp_merge <- merge(snp2[,c(3,17)], snp3[,c(3,17)], by = "variantID", all = TRUE)
snp_merge <- snp_merge[,1:2]

#Loop over all the samples, with merge() we make sure every unique SNP id is present
for (j in 1:(length(samples))) {
  snp2 <- samples[[j]]
  snp_merge <- merge(snp_merge, snp2[,c(3,17)], by = "variantID", all = TRUE)
}
Dose_matrix <- snp_merge
Dose_matrix <- Dose_matrix[,-2]
colnames(Dose_matrix) <- samples2

#Make sure all SNP IDs are unique
length(unique(snp_merge$variantID))

#----------------------------------------------------------------------
#Create a list of all SNPs that are present within all samples, which will be used for MatrixEQTL analysis
snp.list <- merge(samples[[1]][,1:4], samples[[2]][,1:4], by = c("contig", "variantID", "position"), all = T)
for (j in 1:length(samples)) {
  sample <- samples[[j]]
  snp.list <- merge(snp.list, sample[, c(1,2,3,4)], by = c("contig", "variantID", "position"), all = T)
}
snp.list <- snp.list[,1:3]

#----------------------------------------------------------------------
#This is where the SNPs and GeneIDs are annotated with the gene names

library(biomaRt)
# Set Mart used to map Ensembl Gene IDs to Gene name
grch38 = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

#Mapping Ensembl Gene IDs obtained from the eQTL analysis to Gene names
eqtls.gene <- getBM(attributes = c("ensembl_gene_id_version", "chromosome_name","external_gene_name","start_position","end_position","strand"),
                    filters = "ensembl_gene_id_version", 
                    values =  eqtls.cvrt$gene, 
                    mart = grch38,
                    verbose = F)


#This is the filtering step for the Dose_matrix, NA<4, MAF > 0.05
dose.maf <- merge(Dose_matrix, maf[,c(1,7)], by.x = "variantID", by.y = "X", all = T) #merge the dataframes
dose.maf <- dose.maf[rowSums(is.na(dose.maf)) <4, ] #filter out NA>3
dose.maf <- dose.maf[dose.maf$maf >0.05, ] #MAF filter
dose.maf <- dose.maf[,-39]
write.csv(dose.maf, "snp_maf.csv", row.names = F, col.names = T)

#----------------------------------------------------------------------
#This is a quick fix for running the allele dose twice.
samples.test <- samples
for (j in 1:length(samples.test)) {
  sample <- samples[[j]]
  sample$alleledose <- sample$alleledose/2
  samples.test[[j]]<-sample
}

#----------------------------------------------------------------------
#This will create the altAllele.ratio

#
snp2 <- samples[[1]]
snp3 <- samples[[2]]

#Initialize the master list of SNPs by merging the first 2, this will start a template
SNP_full_list <- merge(snp2[,c(3, 19)], snp3[,c(3, 19)], by = "variantID", all = TRUE)
SNP_full_list[,c(2, 3)] <- NULL

#This is where the template is merged with the remainder of the of SNP column of the samples
for (j in 1:(length(samples))) {
  snp2 <- samples[[j]]
  SNP_full_list <- merge(SNP_full_list, snp2[,c(3, 19)], by = "variantID", all = TRUE)
}
altRatio_matrix <- SNP_full_list
colnames(altRatio_matrix) <- samples2

#--------------------------------------------------------------------
#This is to calculate the beta statistic
eqtls.cis.cvrt$BetaStat <- eqtls.cis.cvrt$beta/eqtls.cis.cvrt$statistic

#
snp2 <- samples[[1]]
snp3 <- samples[[2]]

#Initialize the master list of SNPs by merging the first 2, this will start a template
SNP_full_list <- merge(snp2[,c(3, 19)], snp3[,c(3, 19)], by = "variantID", all = TRUE)
SNP_full_list[,c(2, 3)] <- NULL

#This is where the template is merged with the remainder of the of SNP column of the samples
for (j in 1:(length(samples))) {
  snp2 <- samples[[j]]
  SNP_full_list <- merge(SNP_full_list, snp2[,c(3, 19)], by = "variantID", all = TRUE)
}
altRatio_matrix <- SNP_full_list
colnames(altRatio_matrix) <- samples2

#----------------------------------------------------------------------
#This is to obtain a matrix of the total counts per snp per patient
snp2 <- samples[[1]]
snp3 <- samples[[2]]

#Initialize the master list of SNPs by merging the first 2, this will start a template
SNP_full_list <- merge(snp2[,c(3, 8)], snp3[,c(3, 8)], by = "variantID", all = TRUE)
SNP_full_list[,c(2, 3)] <- NULL

#This is where the template is merged with the remainder of the of SNP column of the samples
for (j in 1:(length(samples))) {
  snp2 <- samples[[j]]
  SNP_full_list <- merge(SNP_full_list, snp2[,c(3, 8)], by = "variantID", all = TRUE)
}
totalCount <- SNP_full_list
colnames(totalCount) <- samples2

#------------------------------------------------------------------------
#This will make a genotype matrix
# This part will create the master list of all the SNPs amongst all samples, and then will create the dose matrix.
#Initialize the master list
snp2 <- samples[[1]]
snp3 <- samples[[2]]

snp_merge <- merge(snp2[,c(3,16)], snp3[,c(3,16)], by = "variantID", all = TRUE)
snp_merge <- snp_merge[,1:2]

#Loop over all the samples, with merge() we make sure every unique SNP id is present
for (j in 1:(length(samples))) {
  snp2 <- samples[[j]]
  snp_merge <- merge(snp_merge, snp2[,c(3,16)], by = "variantID", all = TRUE)
}
gt_matrix <- snp_merge
gt_matrix <- gt_matrix[,-2]
colnames(gt_matrix) <- samples2

#--------------------------------------------------
# This is to impute missing values of the dose matrix used for the eQTL analysis

library(pcaMethods)
dose.matrix.maf <- read.csv("/Users/hilmi/DCMUM/snp_maf.csv", header = T, row.names = 1)
dose.matrix <- as.matrix(dose.matrix.maf)
dose.matrix <- prep(t(dose.matrix), scale = "none", center = T)
resSVDI <- pca(dose.matrix, method = "svdImpute", center = F, nPcs = 5)
biplot(resSVDI, choices = 1:2, scale= 1, pc.biplot = F, main = "Dose Matrix Impute PCA")
plot(resSVDI, main = "Dose Matrix Imupte PCA")
plotPcs(resSVDI, pcs = 1:5, sl = NULL, type = "scores", hotelling = 0.95, main = "Dose matrix PCA 1:5")
dose.matrix.impute <- svdImpute(dose.matrix, nPcs = 5, threshold = 0.01, maxSteps = 100)



