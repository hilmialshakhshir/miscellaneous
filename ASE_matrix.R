#This is the loop that will add new columns to determine the Genotype
for (j in 1:length(samples)) {
  ASE <- samples[[j]]
  ASE$gene1 <- NA
  ASE$gene2 <- NA
  ASE$Genotype <- NA
  ASE$alleledose <- 0
  #This is the loop that determines the genotype for every sample
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
  #This will create a column that has the combined alleles
  ASE$Genotype <- paste(ASE$gene1, ASE$gene2)
  samples[[j]] <- ASE
}

***********************************************************************************************
  #This is the loop that calculates the allele dose per SNP
  for (j in 1:length(samples)) {
    ASE <- samples[[j]]
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

***********************************************************************************************
  
  
  for (j in 1:length(samples)) {
    ASE <- samples[[j]]
    ASE$refRatio
    ASE$altRatio
    ASE$refRatio <- ASE$refCount/ASE$totalCount
    ASE$altRatio <- ASE$altCount/ASE$totalCount
    samples[[j]] <- ASE
  }

#This is the configuration to get the matrix with the allele frequency count
snp2 <- samples[[1]]
snp3 <- samples[[2]]

#Initialize the master list of SNPs by merging the first 2, this will start a template
snp_merge <- merge(snp2[,c(3, 18)], snp3[,c(3, 18)], by = "variantID", all = TRUE)
snp_merge[,c(2,3)] <- NULL

#This is where the template is merged with the remainder of the of SNP column of the samples
for (j in 1:(length(samples))) {
  snp2 <- samples[[j]]
  snp_merge <- merge(snp_merge, snp2[,c(3, 18)], by = "variantID", all = TRUE)
}

#This point is to check the previous outcome only has unique variables
length(unique(snp_merge$variantID))

#This is where the matrix is formed, in the first element of the merge
for (j in 1:length(samples)) {
  ase1 <- samples[[j]]
  ase1 <- merge(snp_merge, ase1[,c(3, 17)], by = "variantID", all = TRUE)
  samples_dose[[j]] <- ase1
}

***********************************************************************************************
  #This is the configuration to get the matrix with the allele frequency count
  snp2 <- samples[[1]]
snp3 <- samples[[2]]

#Initialize the master list of SNPs by merging the first 2, this will start a template
SNP_full_list <- merge(snp2[,c(3, 17)], snp3[,c(3, 17)], by = "variantID", all = TRUE)
SNP_full_list[,c(2, 3)] <- NULL

#This is where the template is merged with the remainder of the of SNP column of the samples
for (j in 1:(length(samples))) {
  snp2 <- samples[[j]]
  SNP_full_list <- merge(SNP_full_list, snp2[,c(3, 17)], by = "variantID", all = TRUE)
}
#The variables are deleted so that in the next loop each variable is added with each new loop
SNP_full_list_2 <- SNP_full_list
SNP_full_list_2[,c(2:36)] <- NULL

#This point is to check the previous outcome only has unique variables
length(unique(snp_merge$variantID))
Dose_matrix[,c(2:37)] <- NULL

#This is where the matrix is formed, in the first element of the merge
for (j in 1:length(samples)) {
  ase1 <- samples[[j]]
  Dose_matrix <- merge(Dose_matrix, ase1[,c(3, 17)], by = "variantID", all = TRUE)
}

colnames(Dose_matrix) <- samples1


#This is the configuration to get the matrix with the allele frequency count
snp2 <- samples[[1]]
snp3 <- samples[[2]]

#Initialize the master list of SNPs by merging the first 2, this will start a template
SNP_full_list <- merge(snp2[,c(3, 18)], snp3[,c(3, 18)], by = "variantID", all = TRUE)
SNP_full_list[,c(2, 3)] <- NULL

#This is where the template is merged with the remainder of the of SNP column of the samples
for (j in 1:(length(samples))) {
  snp2 <- samples[[j]]
  SNP_full_list <- merge(SNP_full_list, snp2[,c(3, 18)], by = "variantID", all = TRUE)
}
refRatio_matrix <- SNP_full_list
colnames(refRatio_matrix) <- samples1


#This is the configuration to get the matrix with the allele frequency count
snp2 <- samples[[1]]
snp3 <- samples[[2]]

#Initialize the master list of SNPs by merging the first 2, this will start a template
SNP_full_list <- merge(snp2[,c(3, 19)], snp3[,c(3, 19)], by = "variantID", all = TRUE)
SNP_full_list[,c(2, 3)] <- NULL

#ALTRATIO
#This is where the template is merged with the remainder of the of SNP column of the samples
for (j in 1:(length(samples))) {
  snp2 <- samples[[j]]
  SNP_full_list <- merge(SNP_full_list, snp2[,c(3, 19)], by = "variantID", all = TRUE)
}
altRatio_matrix <- SNP_full_list
colnames(altRatio_matrix) <- samples1

#This is the configuration to get the matrix with the allele frequency count
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
otherbases_matrix <- SNP_full_list
colnames(totCount_matrix) <- samples1
