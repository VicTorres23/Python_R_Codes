## Introduction to clustering - clustering of expression data

##  Go through this script line by line to read and understand the
##  code. 
## Execute code by typing <cmd><enter> or <ctrl><enter>. When nothing is
##  selected, that will execute the current line and move the cursor to
##  the next line. You can also select more than one line, e.g. to
##  execute a block of code, or less than one line, e.g. to execute
##  only the core of a nested expression.
##
##  Edit code, as required, experiment with options, or just play.
## Especially play around withthe options and parameters after you get an understanding.

##  DO NOT simply source() this whole file!
##  If there are portions you don't understand, use R's help system first or if the question is unresolved, ask me. 

# ====================================================================
# ==================================================
# Data
# ==================================================
#
## Let's find some cell-cycle data in GEO, for clustering.
## The goal is to identify coregulated genes, but we don't
## know what their response to time in the cell-cycle will be. Going up or going down?

## The first part of code is slightly adapted from
## performing a standard GEO2R analysis on the NCBI website for
## "Cell cycle expression profiles in HeLa cells" (GSE26922)
## see: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26922
## The dataset contains triplicate measurements for t0 (blocked) and
## t= 2,4,6,8 and 12h post block-release.

## First, you need to install some analysis packages from bioconductor if you have not installed already
## The following commands do for the bioconductor
# install.packages() does for CRAN.

#3 Then we load the libraries....
library(Biobase)
library(GEOquery)
library(limma)

## load the data using the below command
## Make sure you have saved the .RData file in the same folder
#setwd("C:/Users/sgayathrina/OneDrive - The University of Texas at El Paso/1_Sharan_PhD/Coursework/2024_SEM4_Spring/TA2024/Lab11_2024/MyRuns/Section 1") 
setwd("C:/Users/vemma/Documents/Master in Bioinformatics/Spring_2024_Semester/Introduction to Bioinformatics II/Lab 11")
##set this to your own working directory (use your own path in your computer!!!!)

load("GSE26922.RData")


# Check what we have

head(gset)
str(gset)


# The code below is pretty much verbatim GEO2R ...
# ===== (without detailed explanation) ===========
# Make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Group names for all samples
sml <- c("G0","G0","G0","G1","G1","G1",
         "G2","G2","G2","G3","G3","G3",
         "G4","G4","G4","G5","G5","G5");

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex)
}

# Set up the data and proceed with analysis
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G5-G0, G1-G0, G2-G1,
                             G3-G2, G4-G3, G5-G4,
                             levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

# load NCBI platform annotation
gpl <- annotation(gset)
platf <- getGEO(gpl, AnnotGPL=TRUE)
ncbifd <- data.frame(attr(dataTable(platf), "table"))

# replace original platform annotation
tT <- tT[setdiff(colnames(tT), setdiff(fvarLabels(gset), "ID"))]
tT <- merge(tT, ncbifd, by="ID")
tT <- tT[order(tT$P.Value), ]  # restore correct order

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value",
                          "F","Gene.symbol","Gene.title"))

# =========================================
# so far, the GEO2R code ...

# It has returned to us the 250 top-differentially expressed
# genes across the groups we have defined. Note though, that
# the statistics here implicitly treat the data points as
# independent, not as time-series. However for our purpose
# of demonstrating clustering methods this is fine.

# ==================================================
# Exploring the data
# ==================================================
#
# Let's familiarize ourselves a bit with the structure
# of the data.

head(tT)

# The top gene has the ID 8117594: what are the original values?

#exprs(gset)["8117594",]
#barplot(exprs(gset)["8117594",])

#exprs(gset)["7919642", "GSM662899"]
#exprs(gset)["7919642", 5]

exprs(gset)["7900167",]
barplot(exprs(gset)["7900167",])

exprs(gset)["7900167", "GSM662895"]
exprs(gset)["7900167", 7]


# Note how I use a string constant to get a data row from the table
# of expression values. We can also use a vector - the expression below
# returns the data rows for the top three differentially expressed genes:
# exprs(gset)[c("8117594","7900167", "8151871"),]

exprs(gset)[c("8062571","7919642", "7969374"),]

# ==================================================
# Processing the data for cluster analysis
# ==================================================
#
# For cluster analysis, it's useful to make a table from
# these data that contains only numbers, and just a single
# value (mean) for the biological replicates.

gSym <- tT$Gene.symbol

dat <- c()
for (i in 1:nrow(tT)) {
  v <- c()
  v  <- c(v,  mean(exprs(gset)[tT$ID[i], 1:3]))  # t = 0
  v  <- c(v,  mean(exprs(gset)[tT$ID[i], 4:6]))  # t = 2 hours
  v  <- c(v,  mean(exprs(gset)[tT$ID[i], 7:9]))  # etc...
  v  <- c(v,  mean(exprs(gset)[tT$ID[i], 10:12]))
  v  <- c(v,  mean(exprs(gset)[tT$ID[i], 13:15]))
  v  <- c(v,  mean(exprs(gset)[tT$ID[i], 16:18]))
  dat <- rbind(dat, v)
}
colnames(dat) <- c("t0", "t2", "t4", "t6", "t8", "t12")

# We could use the IDs as rownames, like so ...
rownames(dat) <- tT$ID
# ... or the gene symbols, since the IDs don't really
# tell us anything useful. But: are the gene symbols unique?
rownames(dat) <- tT$Gene.symbol

# If they are not unique, we'll have all sorts of trouble later
# on when we select by rowname and find multiple instances....
# R has the function duplicated() to find repeated values
# in a vector...
as.character(gSym[duplicated(gSym)])

# See the vount on the left side below! 
# There are eleven symbols that reappear. Some of them  are
# formatted like "FAM72A///FAM72D///FAM72B" which may mean that
# a spot on the microarray doesn't distinguish between three
# isoforms ... and some are simply the empty string "".
# Since duplicated() gives us a convenient logical vector to
# identify them, we can simply remove them. This is good enough
# for our clustering exercise, for "real" work we should go back
# to the platform information, find out why there are duplicated
# gene symbols, and address this issue.
dat <- dat[!duplicated(gSym), ]
rownames(dat) <- gSym[!duplicated(gSym)]

# We'll also remove all rows that have spots for isoforms.
dat <- dat[-(grep("/", rownames(dat))), ]

# This completes the creation of our expression dataset for clustering.

# You could store the data in a local file ...

write.csv(dat, file="GSE26922.dat")

# and then read it back like so if you want to analyse further....
dat2 <- as.matrix(read.csv(file="GSE26922.dat",
                           row.names = 1,
                           header=TRUE))
# rm(dat)
## ... or, you could save the object as a binary object
## using the function saveRDS(). Then you can read it
## back in with readRDS(). Note that you can change the
## name when you read it back.

saveRDS(dat,  file="GSE26922.rds")

dat <- readRDS("GSE26922.rds")
# identical(dat, dat2)  # has to be TRUE !

# ==================================================
# First explorations: Heatmap
# ==================================================
#
# Heatmaps are a staple of gene expression analysis.
# You can tweak many of the parameters, but for a first look
# we'll just heatmap the data with default parameters.

# This is a standard view that can be applied to all manners
# of multidimensional data, not just genes.
heatmap(dat, Colv = NA)

# Just for illustration and readability let's map only
# every fifth gene
heatmap(dat[seq(1, nrow(dat), by=5), ], Colv = NA, cexRow = 0.3)

#### I'm changing the code to map every 3th gene.

heatmap(dat[seq(1, nrow(dat), by=3), ], Colv = NA, cexRow = 0.3)

# And error happened and I had to run this function dev.off() to see the plots again.

# What's the actual range of values?
range(dat[,1])

# Study the heatmap, and consider what it tells you.
# For example, there seem to be genes that are low at t4 and t6
# but high at t0, t2 ...
#set1 <- c("FNIP1", "MED13L", "NRIP1", "MSI2", "ZNFX1")

#### Here I have modified the vector with genes that are low at t0 and t2
#### and high at t6 and t8.
set1 <- c("CMIP", "PLCL2", "DDX58", "LDLRAD3", "POU2F1")
# ... and there are genes for which the inverse is true:
#set2 <- c("FBXL20", "CCNE1", "ZBTB14", "HIST1H2AH")

#### Genes that are high at t0 and t2 but low at t6 and t8.

set2 <- c("NRG4", "RFC4", "PEX11B", "FAM111B")

# We can use a "parallel coordinates" plot - matplot()
# to look at the actual expression levels. Note that
# matplot expects the values column-wise ordered, thus
# we have to transpose - t() - the data!
matplot(t(dat[set1,]),
        type="l", lwd=2, col="skyblue", lty=1,
        ylim=c(8,14), xlab="time", ylab="log expression value")

# Then we can use lines() to superimpose the genes for set2.
# No transpose here
for (i in 1:length(set2)) {
  lines(dat[set2[i], ], type="l", lwd=2, col="firebrick")
}

# Indeed, these genes - visibly different in the heatmap
# have similar expression profiles within sets and different
# profiles between sets.

# ==================================================
# Hierarchical clustering
# ==================================================
#
# Hierarchical clustering is probably the most basic technique.
# The dendograms on the rows and columns of the heatmap
# were created by hierarchical clustering.

# For hierarchical clustering, first we need to produce
# a distance table. There are many ways to define distances
# let's just go with the default: "Euclidean distance".
#distDat <- dist(dat)

#### I will use the Canberra method to build the dendogram:
distDatCanberra <- dist(dat, method = "canberra")

# Then we use the clustering distance matrix to produce a
# dendrogram in which the most similar genes are connected, and then
# similar genes or connected groups are added. There are
# several ways to define "most-similar", lets just go with the
# default for now: "complete linkage" hierarchical clustering
hcCanberra <- hclust(distDatCanberra)

plot(hcCanberra)

# But do note that both distance as well as clustering
# method matter, and there is not really a "best" way that
# works for all data. You'll need to explore: what you are looking for
# is a distance metric that gives the clearest block structure.

#dEu <- function(x) dist(x, method="euclidian")
#heatmap(dat, distfun = dEu, Colv = NA)

#dCan <- function(x) dist(x, method="canberra")
#heatmap(dat, distfun = dCan, Colv = NA)

#dMax <- function(x) dist(x, method="maximum")
#heatmap(dat, distfun = dMax, Colv = NA)

#dMink <- function(x) dist(x, method="minkowski")
#heatmap(dat, distfun = dMink, Colv = NA)

dEu <- function(x) dist(x, method="euclidian")
heatmap(dat, distfun = dEu, Colv = NA, main = "Euclidian Distance")

dCan <- function(x) dist(x, method="canberra")
heatmap(dat, distfun = dCan, Colv = NA, main = "Canberra Distance")

dMax <- function(x) dist(x, method="maximum")
heatmap(dat, distfun = dMax, Colv = NA, main = "Maximum Distance")

dMan <- function(x) dist(x, method="manhattan")
heatmap(dat, distfun = dMan, Colv = NA, main = "Manhattan Distance")

dMink <- function(x) dist(x, method="minkowski", p=1)
heatmap(dat, distfun = dMink, Colv = NA, main = "Minkowski Distance P=1")

dMink <- function(x) dist(x, method="minkowski", p=2)
heatmap(dat, distfun = dMink, Colv = NA, main = "Minkowski Distance P=2")

dMink <- function(x) dist(x, method="minkowski", p=3)
heatmap(dat, distfun = dMink, Colv = NA, main = "Minkowski Distance P=3")

# You are not confined to the default distance functions, it
# is quite straightforward to define your own, for example
# using correlation properties. Here is a distance function
# defined as 1 - abs(pearson correlation)...

dCor <- function(x) as.dist(1 - abs(cor(t(x))))
heatmap(dat, distfun = dCor)


# ==================================================
# Partitioning clustering
# ==================================================

# === K-means ======================================

# K-means clusters by assigning elements to a fixed
# number of cluster centres, so that similarity
# within a cluster is maximized.

?kmeans


## Uncomment the below line if the RColorBrewer package is not already installed

install.packages("RColorBrewer")

library(RColorBrewer)

#k <- 3
#cl <- kmeans(dat, k)

#niceCols <- brewer.pal(k, "Spectral")

#plot(dat[,"t2"], dat[,"t6"], col = niceCols[cl$cluster])
#points(cl$centers, col = niceCols[1:k], pch = 8, cex=2)

k <- 6
cl <- kmeans(dat, k)

niceCols <- brewer.pal(k, "Spectral")

plot(dat[,"t2"], dat[,"t12"], col = niceCols[cl$cluster])
points(cl$centers, col = niceCols[1:k], pch = 8, cex=2)

# === K-medoids ======================================

# load library "cluster" for K-medoid partitioning
if (!require(cluster, quietly=TRUE)) {
  install.packages("cluster")
  library(cluster)
}

set.seed(112358)

#k <- 6
#cl<-pam(dat, 4)
#plot(dat[,"t2"],dat[,"t12"], col=niceCols[cl$cluster])

#plot(cl)

k <- 6
cl<-pam(dat, 4)
plot(dat[,"t2"],dat[,"t12"], col=niceCols[cl$cluster])

plot(cl) # shows boundary and silhouette plots

## Have fun coding!
## Sharan Gayathrinathan wishes the best :)








