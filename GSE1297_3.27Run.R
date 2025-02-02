# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
#install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("umap")


#Install the required Libraries
library(GEOquery)
library(limma)
library(umap)
library(dplyr)


# load series and platform data from GEO
gset <- getGEO("GSE80178", GSEMatrix =TRUE, AnnotGPL=TRUE)
gset
if (length(gset) > 1) idx <- grep("GPL16686", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
head(gset)
feature_labels <- fvarLabels(gset)
feature_labels
# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "0000000000001011111000010000110"
sml <- strsplit(gsms, split="")[[1]]

ex <- exprs(gset)
ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) # log2 transform

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("HIPP_AD","HIP_CTRL"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

nall <- nrow(gset)
gset <- gset[complete.cases(exprs(gset)), ]

# calculate precision weights and show plot of mean-variance trend
v <- vooma(gset, design, plot=T)

# OR weights by group
# v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
v$genes <- fData(gset) # attach gene annotations

# fit linear model
fit  <- lmFit(v)

# set up contrasts of interest and recalculate model coefficients

cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))

tT

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file="GSE1297 Hippocampus AD v NonAD", row.names=F, sep="\t")


################## Modify from here  ######################

## Assign missing values as NA and then delete them
tT[tT==""] <- NA
tT <- na.omit(tT)
tT$Gene.symbol <- gsub("\\///.*","",tT$Gene.symbol)
## Write the output to table
write.table(tT, file="GSE1297_HIP_Final.tsv", row.names=F, sep="\t")

################################################################
## Check the dimension and insert in table
dimension<- dim(gset)
table(gset$group)
################################################################


## Add here filtering (logFC and p-value)
#################################################################################################

##sgayathrina's Change; get Up Regulated and Down Regulated
# Relaxing the criteria for testing different filkters
#GSE1297_up_test <- tT %>%
#  filter(logFC > 0.5 & adj.P.Val < 0.2)

#GSE1297_down_test <- tT %>%
#  filter(logFC < -0.5 & adj.P.Val < 0.2)

GSE1297_up_test <- tT %>%
  filter(logFC > 1 & adj.P.Val < 0.3)

GSE1297_down_test <- tT %>%
  filter(logFC < -1 & adj.P.Val < 0.3)

# Check if the relaxed criteria return any rows
print(nrow(GSE1297_up_test))
print(nrow(GSE1297_down_test))
# shows the number or rows in the console window, is satisfied proceed



#Proceed to write them in the files that will be saved into the folder.
write.table(GSE1297_up_test, file="GSE1297_UR_1", row.names=F, sep="\t")
write.table(GSE1297_down_test, file="GSE1297_DR_1", row.names=F, sep="\t")


#Adjusting for p.adj and FDR - proceed

