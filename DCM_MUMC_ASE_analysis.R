#=============================================================================#
# DCM_MUMC_ASE_analysis.R           										  #
#																			  #
# Version: 1.0   															  #
# Date: Apr 24, 2018											       		  #
# Author: Michiel Adriaens, PhD; MaCSBio, Maastricht University               #
# History:																	  #
#  1.0: Creation															  #
#                                                                             #
#=============================================================================#

# Work on prioritization first: how to find interesting genes. A network based
#  approach seems like a nice approach for this (both from an analytical as 
#  well as a visual point of view). Then add to this what it means, i.e. 
#  how does genetic regulation, epigenetic regulation, splicing play into
#  imbalance. 

# Thought: If we treat the concept of allelic imbalance as symmetrical, as we 
#       do with the ratioMAD, then anything that is not imbalanced is zero.
#       May help with batteling missing values and ease PCA/clustering/network
#       based analyses. 

# Idea:  NMD can be verified/hypothesized by looking up (known?) tvs in WES

# Issue: How to define ?interesting? imbalanced SNPs/genes: 
#        is it SNPS/genes with many imbalanced individuals? Or with imbalance 
#        in only a few cases? Or with opposite imbalance across individuals? 
#        and how it could contribute to stratification of the patient group.
# Idea:	 Define ASE classes? E.g.: all imbalanced in same direction (~NMD), 
#        imbalanced in opposite directions, imb. in many, imb. in few;
#        Combined with functional annotation (location? gene function? heart
#        muscle specificness?)
# Idea:  Network approach on all genes with at least one imbalanced individual 
#        to find interesting genes and paths between genes: 
#        - Edge = number of shared imbalanced individuals relative to unique
#           number of heterozygotes, i.e. a percentage -> weighted edges
#        - Node color is average degree of imbalance (absolute? Or relative 
#           to alt allele? (I.e. blue and orange) split node color into two?)

# Idea: Test:
#        eQTL enrichment especially for imbalance within 0.5 +/-
#         median ratioMAD; same for epigenetic control/environmental 
#         influence. 
#        NMD more pronounced (likely) outside 0.5 +/- median ratioMAD; 
#         same with sQTL. (checkable by looking for tvs (in WES?))

# Idea: I'm wondering whether we shouldn't also perform an analysis only on 
#       candidate genes? E.g. one of the strikingly imbalanced SNPs, rs2296612,
#       is in NEBL, which is a cardiac muscle gene. 

# Idea: Limit PCA/clustering to smaller set of SNPs/genes instead of all?
#       E.g. focusing on SNPs with high number of hets, or on SNPs within
#       cardiac muscle genes, or DCM genes, or all SNPs with at least 1
#       significantly imbalanced individual, or SNPs with at most 4 sign.
#       imbalanced individuals? Still worried that PCA/clustering driven
#       by missing/blinded values, i.e. essentially looking at genotypic
#       differences, in which case:
# Idea: why not keep in the alt ratios for the homozygous SNPs as well for a 
#       PCA analysis? Or would that be introducing false info? 

# Idea: since allelic imbalance often found to be specific per individual
#       (e.g. NMD due to rare tv), part of the analysis should focus on those
#       SNPs/genes with only very few (1 or 2?) individuals

# FIXME: In binomial test: previously used read depth set to average
#        coverage (then I need the other matrix from Hilmi), perhaps should
#        do that again. Alternative: switch to average read depth PER SNP,
#        or alternatively individual read depth (i.e. per SNP and per 
#        individual) and combine that with the MAD-ratio to identify genes of
#        interest per individual. 
# From 'Allelic_imbalance_analysis_BinomNmeanCount_strict.R' (Berlin)
# "Briefly: large differences in number of counts are creating large differences
#          in power. I.e. highly expressed genes are significantly imbalanced
#          even if there is only a slight deviation from ~0.5, while for
#          moderately expressed genes, this difference needs to be much larger
#          for it to be called as significantly imbalanced. On the one hand 
#          logical, on the other hand something that needs to be adjusted for.
#          One option is to do the filtering afterwards (i.e. demand for
#          differential ASEs that the alt-ref difference is at least 0.20),
#          another is to set the number of trials in the binomial test to
#          the mean number of trials, i.e. the mean of counts of all 
#          potentially imbalanced genes, ALWAYS (that's the difference with
#          the other script)."

#-----------------------------------------------------------------------------#
# 1. Initialization
#-----------------------------------------------------------------------------#

# Set directories
DATA.DIR <- "~/macsbio_research/DCM_RNA_WES/"
RES.DIR <- "~/macsbio_research/DCM_RNA_WES/Results"
		
# Load packages
require(ggplot2)
require(qvalue)
require(pcaMethods)

# Load data
setwd(DATA.DIR)
load("dose_matrix.RData") 		# allele dosage matrix
load("altRatio_matrix.RData")   # alternative allele ratio matrix
load("totalCount.RData")		# read count matrix

# Some checks
all(Dose_matrix[, 1] == altRatio_matrix[, 1]) # TRUE
all(colnames(Dose_matrix) == colnames(altRatio_matrix)) # TRUE
all(Dose_matrix[, 1] == totalCount_matrix[, 1]) # TRUE
all(colnames(Dose_matrix) == colnames(totalCount_matrix)) # TRUE

# Filter data: require genotyping coverage of at least 90% across samples
#  and leave only SNPs with at least 1 heterozygote
nMissing  <- apply(Dose_matrix, 1, function(x) sum(is.na(x[-1])))
nHet <- apply(Dose_matrix, 1, function(x) {
			sum(as.numeric(x[-1]) == 1, na.rm = T)
		})
selection <- which(nMissing < 4 & nHet > 0)
dosageData   <- Dose_matrix[selection, -1]
altRatioData <- altRatio_matrix[selection, -1]
totalCountData <- totalCount_matrix[selection, -1]
altRatioData[dosageData == 0 | dosageData == 2] <- NA 
rownames(altRatioData) <- rownames(dosageData) <- 
 rownames(totalCountData) <- Dose_matrix[selection, 1]

# Global median altRatio:
altRatioMedians <- apply(altRatioData, 2, median, na.rm = T) # all < 0.5
median(altRatioMedians, na.rm = T) # 0.444

# Basic plotting function
plotAse <- function(snp) {
	png(paste("AllelicImbalancePlot_", snp, ".png", sep = ""), 
			height = 1200, width = 1200,
			res = 188)
	par(las = 2)
	par(oma = c(2,0,0,0))
	.temp <- altRatioData[snp, ][, !is.na(altRatioData[snp, ])]
	asePlotData <- rbind(sort(.temp), 1 - sort(.temp))
	barplot(as.matrix(asePlotData), main = paste("Allelic imbalance", snp),
			xlab = "", col = c("darkblue","orange"), ylab = 'Ratio',
			ylim = c(0, 1.2))
	legend("topleft", c('Alternative ratio', 'Reference ratio'), 
			fill = c("darkblue","orange"))
	lines(c(-1E9, 1E9), c(0.5, 0.5), lty = 3, col = "darkgrey", lwd = 2)
	lines(c(-1E9, 1E9), c(0.5 - IS.MEDIAN, 0.5 - IS.MEDIAN), 
			lty = 3, col = "lightgrey", lwd = 2)
	lines(c(-1E9, 1E9), c(0.5 + IS.MEDIAN, 0.5 + IS.MEDIAN), 
			lty = 3, col = "lightgrey", lwd = 2)
	dev.off()
}

#-----------------------------------------------------------------------------#
# 2. Global prioritization metrics
#-----------------------------------------------------------------------------#
# Number of heterozygotes
#------------------------
numberOfHets <- apply(altRatioData, 1, function(x) sum(!is.na(x)))

# Median absolute deviation from 0.5
#-----------------------------------
ratioMad <- function(x) {
      median(abs(x[!is.na(x)] - 0.5)) # should 'expected' be changed to 0.444?
}
imbalanceScore  <- apply(altRatioData, 1, ratioMad)
IS.MEDIAN <- median(imbalanceScore)

# Seems to be yielding high scores for SNPs with low MAF (i.e. where there's 
#  a high likelihood of observing only very few heterozygous individuals).
#---------------------------------------------------------------------------

# Plot 2.1: Imbalance score density split into high/low number of hets 
#  I'm using the mean of the number of heterozygotes (across all SNPs/genes) 
#   instead of median due to the integer nature of that score:
setwd(RES.DIR)
png("ImbalanceScore_DensityPlot_nHet.png", height = 1200, width = 1200,
		res = 188)
d1 <- density(imbalanceScore,  na.rm = T)
d2 <- density(imbalanceScore[numberOfHets > mean(numberOfHets)],  na.rm = T) 
d3 <- density(imbalanceScore[numberOfHets < mean(numberOfHets)], na.rm = T)
plot(d1, ylim = c(0, max(c(d1$y, d2$y, d3$y))), 
		main = 'Density of imbalance score (ratioMAD)',
		col = "lightgrey", xlab = "Imbalance score (ratioMAD)")
lines(d2, col = "darkblue")
lines(d3, col = "orange")
legend("top", c('All', 'N_het > mean', 'N_het <= mean'), 
		fill = c('lightgrey', 'darkblue', 'orange'))
dev.off()

# Plot 2.2: Number of hets density split into high/low imbalance score
setwd(RES.DIR)
png("nHets_DensityPlot_ImbalanceScore.png", height = 1200, width = 1200,
		res = 188)
d1 <- density(numberOfHets,  na.rm = T)
d2 <- density(numberOfHets[imbalanceScore < median(imbalanceScore, na.rm = T)],
		na.rm = T)          
d3 <- density(numberOfHets[imbalanceScore > median(imbalanceScore, na.rm = T)], 
		na.rm = T)
plot(d1, ylim = c(0, max(c(d1$y, d2$y, d3$y))), 
		main = 'Density of number of heterozygotes', col = "lightgrey", 
		xlab = "Number of heterozygotes")
lines(d2, col = "darkblue")
lines(d3, col = "orange")
legend("top", c('All', 'ImbalanceScore < median', 'ImbalanceScore > median'), 
		fill = c('lightgrey', 'darkblue', 'orange'))
dev.off()

# Also possible to combine number of heterozygotes with the imbalance score
#  to identify interesting genes/SNPs:
#--------------------------------------------------------------------------
id <- rownames(altRatioData[numberOfHets > 5 & 
						    imbalanceScore > IS.MEDIAN, ][1, ])
setwd(RES.DIR)
plotAse(id)
# SNP "rs1041939" implicated in hypoxia in SCD (https://goo.gl/KxsL2p)

# How about least varying SNPs (e.g. using quantiles of imbalanceScore)
#----------------------------------------------------------------------
for (prob in c(0.5, 0.33, 0.25, 0.10)) {
print(median(as.vector(as.matrix(altRatioData[which(imbalanceScore < 
	   quantile(imbalanceScore, prob = prob, na.rm = T)),])),
       na.rm = T))
}
# The median alternative allele ratio get's closer to 1:1 ratio for 
#  higher imbalance score cut-offs.

# IMPORTANT: 
# High heterozygosity gene/snp with only a few, but still highly imbalanced 
# individuals are lost with these global scores, e.g.:
#--------------------------------------------------------------------------
altRatioData['rs1049862',]
median(as.numeric(altRatioData['rs1049862',]), na.rm = T)
imbalanceScore['rs1049862']
imbalanceScore[1:10] # c.f.
setwd(RES.DIR)
plotAse('rs1049862')

# I think there is a difference between using the global scores to find
#  globally interesting genes, and using patient specific scores (i.e.
#  statistically significant imbalance, next section) to identify
#  interesting genes per individual patient. 

# A additional score would be simply counting the number of significantly
#  imbalanced individuals: this can be used as a global score as well
#  (next section)

#-----------------------------------------------------------------------------#
# 3. Individual based prioritization metrics
#-----------------------------------------------------------------------------#

asePvalueMatrix <- c() 
for (i in 1:37) {
	nr  <- floor(rowMeans(totalCountData, na.rm = T)) #* (see comment below)
#	nr  <- totalCountData[, i] #* (see comment below)
	ratio <- altRatioData[, i]
	p <- apply(cbind(ratio, nr), 1, function(r) {
				if (is.na(r[1])) {
					NA
				} else {
					binom.test(round(r[1] * r[2]), r[2], 
						 median(ratio, na.rm = T))$p.value
				}
			})
	asePvalueMatrix <- cbind(asePvalueMatrix, p)
	message(i)
}
colnames(asePvalueMatrix) <- colnames(altRatioData)
rownames(asePvalueMatrix) <- rownames(altRatioData)

# Correct p-values
aseQvalueMatrix <- asePvalueMatrix
aseQvalueMatrix[!is.na(aseQvalueMatrix)] <- 
  qvalue(as.vector(as.matrix(aseQvalueMatrix[!is.na(aseQvalueMatrix)])))$qvalues

# Check some stats
nImbalancedSnps <- c() 
for (i in 1:37) {
	nImbalancedSnps <- c(nImbalancedSnps, 
	  sum(aseQvalueMatrix[, i] < 0.05, na.rm = T))
}
cbind(nImbalancedSnps, colMeans(totalCountData, na.rm = T))
cor(nImbalancedSnps, colMeans(totalCountData, na.rm = T), method = "spearman")

#* Higher average coverage means higher power to detect significant imbalance.
#   when leaving nr defined per individual (and per SNP) this skews the result,
#   then the above correlation is actually ~0.614, quite high!
#  So I now leave the number of reads to vary per SNP/gene but not per 
#   individual, c.f. transcriptomics data analysis (i.e. higher expressed = 
#   more power to detect change)

nImbalancedIndividuals <- c() 
for (i in 1:nrow(aseQvalueMatrix)) {
	nImbalancedIndividuals <- c(nImbalancedIndividuals, 
	  sum(aseQvalueMatrix[i, ] < 0.05, na.rm = T))
}
names(nImbalancedIndividuals) <- rownames(aseQvalueMatrix)

# Let's plot some results for strikingly imbalanced SNPs/genes:
setwd(RES.DIR)
setwd('Striking ASE examples')
for (id in names(which(nImbalancedIndividuals > 19))) {
	plotAse(id)	
}

#-----------------------------------------------------------------------------#
# 4. PCA plots  
#-----------------------------------------------------------------------------#

# FIXME: unsure what this means

## based on imbalance q-values
##----------------------------
#pcaInput <- aseQvalueMatrix
##pcaInput[is.na(pcaInput)] <- 1
#pcaInput <- -log10(pcaInput)
#pcaRes <- pca(t(pcaInput), nPcs = 10, method = 'svdImpute')
#setwd(RES.DIR)
#png("PCA_plot_AllelicImbalance_Qvalue.png", height = 3 * 1200, width = 2 * 1200,
#		res = 144)
#par(mfrow = c(3, 2))
#plot(pcaRes@scores[,1], pcaRes@scores[,2], typ = "n")
#text(pcaRes@scores[,1], pcaRes@scores[,2], labels = colnames(altRatioData))
#plot(pcaRes@scores[,1], pcaRes@scores[,3], typ = "n")
#text(pcaRes@scores[,1], pcaRes@scores[,3], labels = colnames(altRatioData))
#plot(pcaRes@scores[,1], pcaRes@scores[,4], typ = "n")
#text(pcaRes@scores[,1], pcaRes@scores[,4], labels = colnames(altRatioData))
#plot(pcaRes@scores[,2], pcaRes@scores[,3], typ = "n")
#text(pcaRes@scores[,2], pcaRes@scores[,3], labels = colnames(altRatioData))
#plot(pcaRes@scores[,2], pcaRes@scores[,4], typ = "n")
#text(pcaRes@scores[,2], pcaRes@scores[,4], labels = colnames(altRatioData))
#plot(pcaRes@scores[,3], pcaRes@scores[,4], typ = "n")
#text(pcaRes@scores[,3], pcaRes@scores[,4], labels = colnames(altRatioData))
#dev.off()
#
## First compontents ~5% range of variance explained
#
## based on altRatio data
##-----------------------
#pcaRes2 <- pca(t(altRatioData), nPcs = 10, method = 'svdImpute')
#setwd(RES.DIR)
#png("PCA_plot_AllelicImbalance_altRatio.png", height = 3 * 1200, width = 2 * 1200,
#		res = 144)
#par(mfrow = c(3, 2))
#plot(pcaRes2@scores[,1], pcaRes2@scores[,2], typ = "n")
#text(pcaRes2@scores[,1], pcaRes2@scores[,2], labels = colnames(altRatioData))
#plot(pcaRes2@scores[,1], pcaRes2@scores[,3], typ = "n")
#text(pcaRes2@scores[,1], pcaRes2@scores[,3], labels = colnames(altRatioData))
#plot(pcaRes2@scores[,1], pcaRes2@scores[,4], typ = "n")
#text(pcaRes2@scores[,1], pcaRes2@scores[,4], labels = colnames(altRatioData))
#plot(pcaRes2@scores[,2], pcaRes2@scores[,3], typ = "n")
#text(pcaRes2@scores[,2], pcaRes2@scores[,3], labels = colnames(altRatioData))
#plot(pcaRes2@scores[,2], pcaRes2@scores[,4], typ = "n")
#text(pcaRes2@scores[,2], pcaRes2@scores[,4], labels = colnames(altRatioData))
#plot(pcaRes2@scores[,3], pcaRes2@scores[,4], typ = "n")
#text(pcaRes2@scores[,3], pcaRes2@scores[,4], labels = colnames(altRatioData))
#dev.off()
#
## 40% variance explained by 1st, 25% by second component?
#
## based on altRatio deviation from 0.5
##-------------------------------------
#pcaInput3 <- altRatioData - 0.5
#pcaInput3[is.na(pcaInput3)] <- 0
#pcaRes3 <- pca(t(pcaInput3), nPcs = 10)
#setwd(RES.DIR)
#png("PCA_plot_AllelicImbalance_altRatioDeviation.png", height = 3 * 1200, 
#		width = 2 * 1200, res = 144)
#par(mfrow = c(3, 2))
#plot(pcaRes3@scores[,1], pcaRes3@scores[,2], typ = "n")
#text(pcaRes3@scores[,1], pcaRes3@scores[,2], labels = colnames(altRatioData))
#plot(pcaRes3@scores[,1], pcaRes3@scores[,3], typ = "n")
#text(pcaRes3@scores[,1], pcaRes3@scores[,3], labels = colnames(altRatioData))
#plot(pcaRes3@scores[,1], pcaRes3@scores[,4], typ = "n")
#text(pcaRes3@scores[,1], pcaRes3@scores[,4], labels = colnames(altRatioData))
#plot(pcaRes3@scores[,2], pcaRes3@scores[,3], typ = "n")
#text(pcaRes3@scores[,2], pcaRes3@scores[,3], labels = colnames(altRatioData))
#plot(pcaRes3@scores[,2], pcaRes3@scores[,4], typ = "n")
#text(pcaRes3@scores[,2], pcaRes3@scores[,4], labels = colnames(altRatioData))
#plot(pcaRes3@scores[,3], pcaRes3@scores[,4], typ = "n")
#text(pcaRes3@scores[,3], pcaRes3@scores[,4], labels = colnames(altRatioData))
#dev.off()

# END
