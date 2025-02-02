
## Install bioconductor
## Uncomment the line below to install/ Follow the instructions provided in bioconductor.org
## to install the bioconductor package
## source("https://bioconductor.org/biocLite.R")

library(GEOquery)

## Preprocessing and normalization of Affymetr  ix expression data

## 1.set directory which has .CEL files (raw data)
## Change to your working directory
setwd("C:/Users/vemma/Documents/Master in Bioinformatics/Spring_2024_Semester/Introduction to Bioinformatics II/Lab 10")

## 2. Load the "library" that contains the Affymetrix microarray
## Uncomment the below three lines to install the affy package (if you have not installed earlier)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("affy")

library(affy)

## read the CEL files (first command below) and then 
## summarize and normalize with "rma" introduced in the lecture class
## The provided input files (brain and liver cel) present in your directory will be loaded

## 1. Normalization with rma (Robust Multichip Average (RMA) normalization).
affy.data = ReadAffy() 
eset.rma = rma(affy.data)

## 2. Normalization with mas5
## https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-273 
## preprocessing with MAS5
eset.mas5 = mas5(affy.data)

## The variable "eset.rma" contains normalized expression
## values for all probesets, along with other information.
## Let's continue by getting the expression matrix(probesets/genes in rows,
## chips in columns)
exprSet.nologsmas5 = exprs(eset.mas5)

exprSet.nologsrma = exprs(eset.rma)
## List the column (chip) names
colnames(exprSet.nologs)

## Rename the column names if we want
colnames(exprSet.nologsmas5) = c("brain.1","brain.2","fetal.brain.1",
                             "fetal.brain.2",
                             "fetal.liver.1","fetal.liver.2",
                             "liver.1","liver.2")

colnames(exprSet.nologsrma) = c("brain.1","brain.2","fetal.brain.1",
                                 "fetal.brain.2",
                                 "fetal.liver.1","fetal.liver.2",
                                 "liver.1","liver.2")
# check column names
colnames(exprSet.nologsmas5)
colnames(exprSet.nologsrma)

# Now, we'll apply a log transformation to expression values to achieve a 
#normal distribution. This step is crucial for ratio calculations. 
#Using base 2 for logarithms simplifies ratio transformations, 
#turning a 2-fold increase or decrease into +1 or -1, respectively. 
#Thus, we employ base 2 logarithms for straightforwardness.
#You can experiment with others to find the best fit



# 3. log2 transformation
exprSetmas5 = log(exprSet.nologsmas5, 2)
exprSetmas5

exprSetrma = log(exprSet.nologsrma, 2)
exprSetrma

############ End of getting expresion data ####################

## Clustering
brain.datmas5 = exprSetmas5[,c(1:4)]
liver.datmas5 = exprSetmas5[,c(5:8)]

brain.datrma = exprSetrma[,c(1:4)]
liver.datrma = exprSetrma[,c(5:8)]

## 4. Hierarchical Clustering
#d = dist(as.matrix(t(exprSet)))
#hc = hclust(d,method="average")
#plot(hc)

?hclust

d = dist(as.matrix(t(exprSetmas5)))
hc = hclust(d,method="centroid")
plot(hc)

d = dist(as.matrix(t(exprSetmas5)))
hc = hclust(d,method="average")
plot(hc)

d = dist(as.matrix(t(exprSetmas5)))
hc = hclust(d,method="complete")
plot(hc)

d = dist(as.matrix(t(exprSetmas5)))
hc = hclust(d,method="single")
plot(hc)
##########################################################
d = dist(as.matrix(t(exprSetrma)))
hc = hclust(d,method="centroid")
plot(hc)

d = dist(as.matrix(t(exprSetrma)))
hc = hclust(d,method="average")
plot(hc)

d = dist(as.matrix(t(exprSetrma)))
hc = hclust(d,method="complete")
plot(hc)

d = dist(as.matrix(t(exprSetrma)))
hc = hclust(d,method="single")
plot(hc)

## know other options for hclust by running the below command.
?hclust

## 5. K means clustering
#k = 5
k = 4
fitmas5 <- kmeans(t(exprSetmas5), k) # 2 cluster solution
fitmas5$cluster

k = 2
fitmas5 <- kmeans(t(exprSetmas5), k) # 2 cluster solution
fitmas5$cluster

k = 3
fitrma <- kmeans(t(exprSetrma), k) # 2 cluster solution
fitrma$cluster

k = 2
fitrma <- kmeans(t(exprSetrma), k) # 2 cluster solution
fitrma$cluster


## 6. For documentation
?kmeans
#kmeans(exprSet, centers = 3, nstart = 20)
#k1 = kmeans(exprSet, centers = 4, nstart = 20)
#print(k1)
#plot(exprSet, col = k1$cluster, main = "K-means Clustering Results")

k1 = kmeans(exprSetmas5, centers = 2, nstart = 30)
print(k1)
plot(exprSetmas5, col = k1$cluster, main = "K-means Clustering Results")

k1 = kmeans(exprSetmas5, centers = 3, nstart = 30)
print(k1)
plot(exprSetmas5, col = k1$cluster, main = "K-means Clustering Results")

k1 = kmeans(exprSetmas5, centers = 4, nstart = 30)
print(k1)
plot(exprSetmas5, col = k1$cluster, main = "K-means Clustering Results")

k1 = kmeans(exprSetmas5, centers = 8, nstart = 30)
print(k1)
plot(exprSetmas5, col = k1$cluster, main = "K-means Clustering Results")

###################################################################

k1 = kmeans(exprSetrma, centers = 2, nstart = 30)
print(k1)
plot(exprSetrma, col = k1$cluster, main = "K-means Clustering Results")

k1 = kmeans(exprSetrma, centers = 3, nstart = 30)
print(k1)
plot(exprSetrma, col = k1$cluster, main = "K-means Clustering Results")

k1 = kmeans(exprSetrma, centers = 4, nstart = 30)
print(k1)
plot(exprSetrma, col = k1$cluster, main = "K-means Clustering Results")

k1 = kmeans(exprSetrma, centers = 8, nstart = 30)
print(k1)
plot(exprSetrma, col = k1$cluster, main = "K-means Clustering Results")

## 7. PCA

## Principal Components Analysis
## entering raw data and extracting PCs
## from the correlation matrix
fit = princomp((exprSet),cor=TRUE)
summary(fit) # print variance accounted for
loadings(fit) # pc loadings
plot(fit,type="lines") # scree plot
fit$scores # the principal components
biplot(fit)

?princomp

fit = princomp((exprSet),cor=FALSE, scores = FALSE)
summary(fit) # print variance accounted for
loadings(fit) # pc loadings
plot(fit,type="b") # scree plot
fit$scores # the principal components
biplot(fit, pch=8, scale = 0, col = c("blue", "black"))





