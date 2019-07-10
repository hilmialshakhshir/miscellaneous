#Calculate the optimal number of clusters

library(factoextra)
library(cluster)
library(NbClust)

cm <- read.csv("/Users/hilmi/Documents/Internship/countmatrix_p.csv", header = T, row.names = 1)
cm <- as.matrix(cm)
cm <- varianceStabilizingTransformation(cm, blind = TRUE, fitType = "parametric")
cm.scaled <- scale(cm)
cm.scaled <- t(cm.scaled)
#na_dat <- do.call(data.frame,lapply(cm, function(x) replace(x, is.infinite(x),NA)))
#prepare the matrix for analysis
cm.scaled <- scale(cm)
cm.scaled <- t(cm.scaled)
set.seed(123)
# Compute and plot wss for k = 2 to k = 15
k.max <- 15 # Maximal number of clusters
data <- na.omit(cm.scaled)
wss <- sapply(1:k.max, function(k){kmeans(data, k, nstart=8 )$tot.withinss})

plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
abline(v = 3, lty =2)

fviz_nbclust(x, FUNcluster, method = c("silhouette", "wss"))

#k-means
fviz_nbclust(cm.scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 5, linetype = 2)

#pam
fviz_nbclust(cm.scaled, pam, method = "wss") +
  geom_vline(xintercept = 5, linetype = 2)

#h-clustering
fviz_nbclust(cm.scaled, hcut, method = "wss") +
  geom_vline(xintercept = 5, linetype = 2)

library(cluster)
k.max <- 15
data <- cm.scaled
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
     frame = FALSE, xlab = "Number of clusters k")
abline(v = which.max(sil), lty = 2)

fviz_nbclust(cm.scaled, kmeans, method = "silhouette")

fviz_nbclust(cm.scaled, hcut, method = "silhouette",
             hc_method = "complete")


#Calculate the optimal number of clusters using 26 different indices
#nb <- NbClust(cm.scaled, distance = "euclidean", min.nc = 2,
#              max.nc = 10, method = "complete", index ="all")

# K-means clustering
set.seed(1)
cm.scaled <- na.omit(cm.scaled)
km.res <- kmeans(cm.scaled, 5, nstart = 25)
# k-means group number of each observation
km.res$cluster

# Visualize k-means clusters
fviz_cluster(km.res, data = cm.scaled, geom = c("point", "text"), repel = T, ellipse = F,
             stand = FALSE, frame.type = "norm", main = "K-means Clustering of Gene expression Data")

# PAM clustering
library("cluster")
pam.res <- pam(cm.scaled,k= 5)
pam.res$cluster

# Visualize pam clusters
fviz_cluster(pam.res, stand = FALSE, geom = c("point","text"), repel = T, ellipse = F,
             frame.type = "norm")

# Compute pairewise distance matrices
dist.res <- dist(cm.scaled, method = "euclidean")
# Hierarchical clustering results
hc <- hclust(dist.res, method = "complete")
# Visualization of hclust
plot(hc, labels = T, hang = -1)
# Add rectangle around 3 groups
rect.hclust(hc, k = 5, border = 2:4) 


#GAP statisitc
set.seed(123)
gap_stat <- clusGap(cm.scaled, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)

