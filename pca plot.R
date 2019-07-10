library(DESeq2)
countmatrix <- read.csv("/Users/hilmi/Documents/Internship/countmatrix_htseq2.csv", header = T, row.names = 1)
counts <- read.csv("/Users/hilmi/Documents/Internship/countmatrix_htseq2.csv", header = T)
cm <- countmatrix[apply(countmatrix == 0, 1, sum) <= 20, ]
cm.vst <- varianceStabilizingTransformation(as.matrix(cm), blind = T, fitType = "parametric")
km.res <- kmeans(t(cm.vst), 5, nstart = 25)
gx.t <- t(cm.vst)
library(factoextra)
library(cluster)
library(NbClust)
#gx.can.t is the transposed matrix of gene expression for the candidate genes
#gene_vst_2 <- read.csv("/Users/hilmi/DCMUM/gene_vst_2.csv", header = T, row.names = 1)
#km.res <- kmeans(t(gene_vst_2), 4, nstart = 25)
#gx.t <- t(gene_vst_2)
# Visualize k-means clusters
fviz_cluster(km.res, data = gx.t, geom = c("point","text"), repel = T,
             stand = FALSE, frame.type = "norm", ellipse = F, main = "K-means Clustering of Gene expression Data")


library(ggplot2)
library(factoextra)

#plot and compute PCA
res.pca <- prcomp(gx.t, scale = TRUE)
fviz_eig(res.pca, main = "Variance per dimension")
#graph of individuals
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
#graph of variables
#fviz_pca_var(res.pca,
 #            col.var = "contrib", # Color by contributions to the PC
  #           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
   #          repel = TRUE     # Avoid text overlapping
#)


#stats on pca 
# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

gx.pca <- prcomp(gx.t, center = T, scale = F)
scores <- as.data.frame(gx.pca$x)
ggplot(data = scores, aes(x = PC2, y = PC3, label = rownames(scores))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(colour = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot All Genes")


gx.full.t <- t(gx.full)
gx.full.pca <- prcomp(gx.full.t, center = T, scale = F)
scores.full <- as.data.frame(gx.full.pca$x)
ggplot(data = scores.full, aes(x = PC1, y = PC2, label = rownames(scores.full))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(colour = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot All Genes")

ggplot(data = svd.scores, aes(x = PC1, y = PC2, label = rownames(svd.scores))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(colour = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot Dose Matrix")


library(factoextra)
library(cluster)
library(NbClust)

cm.scaled <- scale(gx.t)

# K-means clustering
set.seed(123)
km.res <- kmeans(cm.scaled, 5, nstart = 25)
# k-means group number of each observation
km.res$cluster

# Visualize k-means clusters
fviz_cluster(km.res, data = cm.scaled, geom = "point",
             stand = FALSE, frame.type = "norm")

# PAM clustering
library("cluster")
pam.res <- pam(cm.scaled, 5)
pam.res$cluster

# Visualize pam clusters
fviz_cluster(pam.res, stand = FALSE, geom = "point",
             frame.type = "norm")
# Compute pairewise distance matrices
dist.res <- dist(cm.scaled, method = "euclidean")
# Hierarchical clustering results
hc <- hclust(dist.res, method = "complete")
# Visualization of hclust
plot(hc, labels = FALSE, hang = -1)
# Add rectangle around 3 groups
rect.hclust(hc, k = 5, border = 2:4) 


library(cluster)
k.max <- 15
data <- na.omit(cm.scaled)
sil <- rep(0, k.max)

# Compute the average silhouette width for 
# k = 2 to k = 15
for(i in 2:k.max){
  km.res <- kmeans(data, centers = i, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(data))
  sil[i] <- mean(ss[, 3])
}

# Plot the  average silhouette width
plot(1:k.max, sil, type = "b", pch = 19, 
     frame = FALSE, xlab = "Number of clusters k", main = "Optimal # of clusters - Silhouette")
abline(v = which.max(sil), lty = 2)

# Compute gap statistic
library(cluster)
set.seed(123)
gap_stat <- clusGap(cm.scaled, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
# Base plot of gap statistic
plot(gap_stat, frame = FALSE, xlab = "Number of clusters k", main = "Optimal # of Clusters - GAP")
abline(v = 5, lty = 2)
fviz_gap_stat(gap_stat)


nb <- NbClust(cm.scaled, distance = "euclidean", min.nc = 2,
              max.nc = 8, method = "complete", index ="all")