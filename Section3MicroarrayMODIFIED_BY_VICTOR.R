## Install bioconductor
## Uncomment the line below to install/ Follow the instructions provided in bioconductor.org
## to install the bioconductor package
## source("https://bioconductor.org/biocLite.R")

library(GEOquery)

## Preprocessing and normalization of Affymetr  ix expression data

## 1.set directory which has .CEL files (raw data)
## Change to your working directory
#-------------------CHANGING THE WORKING DIRECTORY------------------------
setwd("C:/Users/vemma/Documents/Master in Bioinformatics/Spring_2024_Semester/Introduction to Bioinformatics II/Lab 10/DFU")
#-------------------------------------------------------------------------
## 2. Load the "library" that contains the Affymetrix microarray
## Uncomment the below three lines to install the affy package (if you have not installed earlier)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("affy")

library(affy)

## read the CEL files (first command below) and then 
## summarize and normalize with "rma" introduced in the lecture class
## The provided input files (brain and liver cel) present in your directory will be loaded

#--------------USING THE OLIGO PACKAGE---------------------------------------
#Here, instead of using Affy I used the Oligo package to read the CEL files
# The affy package kept giving an error but the read.celfiles function
# of the oligo package was able to do the job.
BiocManager::install("oligo")
library(oligo)
## 1. Normalization with rma (Robust Multichip Average (RMA) normalization).
#affy.data = ReadAffy() 
celfiles <- list.celfiles(path = "C:/Users/vemma/Documents/Master in Bioinformatics/Spring_2024_Semester/Introduction to Bioinformatics II/Lab 10/DFU")
affy.data <- read.celfiles(celfiles)
#ffy.data = ReadAffy(path = "C:/Users/vemma/Documents/Master in Bioinformatics/Spring_2024_Semester/Introduction to Bioinformatics II/Lab 10/DFU")
#---------------------------------------------------------------------------

#---------------NORMALIZING THE DATA USING THE OLIGO PACKAGE---------------
eset.rma <- oligo::rma(affy.data)
#--------------------------------------------------------------------------
#eset.rma = rma(affy.data)

## 2. Normalization with mas5
## https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-273 
## preprocessing with MAS5
#eset.mas5 = mas5(affy.data)

## The variable "eset.rma" contains normalized expression
## values for all probesets, along with other information.
## Let's continue by getting the expression matrix(probesets/genes in rows,
## chips in columns)
exprSet.nologs = exprs(eset.mas5)

exprSet.nologs = exprs(eset.rma)
## List the column (chip) names
colnames(exprSet.nologs)

# Now, we'll apply a log transformation to expression values to achieve a 
#normal distribution. This step is crucial for ratio calculations. 
#Using base 2 for logarithms simplifies ratio transformations, 
#turning a 2-fold increase or decrease into +1 or -1, respectively. 
#Thus, we employ base 2 logarithms for straightforwardness.
#You can experiment with others to find the best fit



# 3. log2 transformation
exprSet = log(exprSet.nologs, 2)
exprSet

############ End of getting expresion data ####################

## 4. Hierarchical Clustering
#---------USING DIFFERENT METHODS FOR GETTING THE DENDOGRAMS-----------------
#Here we use different methods of Hierachical Clustering to get different
# dendograms and compare them between each other.
d = dist(as.matrix(t(exprSet)))
hc = hclust(d,method="average")
plot(hc)

?hclust

d = dist(as.matrix(t(exprSet)))
hc = hclust(d,method="centroid")
plot(hc)

d = dist(as.matrix(t(exprSet)))
hc = hclust(d,method="complete")
plot(hc)

d = dist(as.matrix(t(exprSet)))
hc = hclust(d,method="single")
plot(hc)
#----------------------------------------------------------------------------

## know other options for hclust by running the below command.
?hclust

## 5. K means clustering
#k = 5
k = 3
fit <- kmeans(t(exprSet), k) # 2 cluster solution
fit$cluster

k = 2
fit <- kmeans(t(exprSet), k) # 2 cluster solution
fit$cluster


## 6. For documentation
?kmeans
#kmeans(exprSet, centers = 3, nstart = 20)
#k1 = kmeans(exprSet, centers = 4, nstart = 20)
#print(k1)
#plot(exprSet, col = k1$cluster, main = "K-means Clustering Results")
#----MODIFYING THE NUMBER OF CLUSTER AND STARTING POINTS---------
#Here we modify the centers parameter to generate different plots with 
# different number of clusters and starting points.
k1 = kmeans(exprSet, centers = 2, nstart = 30)
print(k1)
plot(exprSet, col = k1$cluster, main = "K-means Clustering Results")

k1 = kmeans(exprSet, centers = 3, nstart = 30)
print(k1)
plot(exprSet, col = k1$cluster, main = "K-means Clustering Results")

k1 = kmeans(exprSet, centers = 4, nstart = 30)
print(k1)
plot(exprSet, col = k1$cluster, main = "K-means Clustering Results")

k1 = kmeans(exprSet, centers = 6, nstart = 30)
print(k1)
plot(exprSet, col = k1$cluster, main = "K-means Clustering Results")
#--------------------------------------------------------------------------
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

rownames(gset)
rownames(tT)
df <- as.data.frame(tT)
rownames(df)