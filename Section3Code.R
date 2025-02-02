library(GEOquery)

library(limma)

library(umap)


#----------CHANGING THE DATASET-------------------------------------------
gset <- getGEO("GSE80178", GSEMatrix =TRUE, AnnotGPL=FALSE)
gset
#-------------------------------------------------------------------------

#-----------CHANGING THE ANNOTATION CORRESPONDING TO THE DATASET-----------
if (length(gset) > 1) idx <- grep("GPL16686", attr(gset, "names")) else idx <- 1


gset <- gset[[idx]]
#--------------------------------------------------------------------------

#----------CHECKING THE COLUMNS THAT THE FILE HAS--------------------------
feature_labels <- fvarLabels(gset)
feature_labels
#-------------------------------------------------------------------------

fvarLabels(gset) <- make.names(fvarLabels(gset))
fvarLabels(gset)

#-------------------CATEGORIZING CONTROL AND DISEASE SAMPLES AND 
#  ASSIGNING AN "X" TO UNWANTED SAMPLES------------------------------------
gsms <- "111111XXX000"

sml <- strsplit(gsms, split="")[[1]]
#--------------------------------------------------------------------------

#-------FILTERING OUT THE UNWANTED SAMPLES---------------------------------

exs <- which(sml != "X")

sml <- sml[exs]

gset <- gset[ ,exs]
#--------------------------------------------------------------------------

#------------------LOG 2 TRANSFORMATION------------------------------------

ex <- exprs(gset)

qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }
#---------------------------------------------------------------------------

#-------------GROUPING THE SAMPLES AND GROUPING CONTROL AND DISEASE-------
gs <- factor(sml)
groups <- make.names(c("CONTROL","DFU"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

head(gset)
#------------------------------------------------------------------------

# ---------CALCULATING PRECISION WEIGHTS AND SHOWING PLOT OF 
# MEAN-VARIANCE TREND----------------------------------------------------
v <- vooma(gset, design, plot=T)

par(mar=c(7,4,2,1))
title <- paste ("GSE80178", "/", annotation(gset), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

title <- paste ("GSE80178", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)
#------------------------------------------------------------------------


# OR WEIGHTS BY GROUP
# v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
v$genes <- fData(gset) # attach gene annotations

#----------------FIT LINEAR MODEL----------------------------------------
fit  <- lmFit(v)

#fit <- lmFit(gset, design)
#------------------------------------------------------------------------

#----SETTING UP CONTRASTS OF INTEREST AND RECALCULATING MODEL
#        COEFFICIENTS----------------------------------------------------

cts <- paste(groups[1], groups[2], sep="-")

cont.matrix <- makeContrasts(contrasts=cts, levels=design)
cont.matrix

fit2 <- contrasts.fit(fit, cont.matrix)
fit2
#------------------------------------------------------------------------

#-----COMPUTE STATISTICS AND GENERATE TABLE OF TOP SIGNIFICANT GENES-----

fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT["16765648",]

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","RANGE_STRAND","RANGE_START","RANGE_END","GB_ACC","SPOT_ID","RANGE_GB"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

head(tT, 20)
tT[order(-tT$logFC), ]
tT[order(tT$logFC), ]

exprs(gset)["16693409",]
barplot(exprs(gset)["16693409",])
#------------------------------------------------------------------------

#--------GENERATING VOLCANO PLOT OF TOP 250 GENES------------------------
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(gset, x = "logFC", y = "adj.P.Val", lab = tT$ID)
#------------------------------------------------------------------------

#--------FILTERING NON-SIGNIFICANT GENES---------------------------------
GSE80178_up_test <- tT %>%
  filter(logFC > 1 & adj.P.Val < 0.3)

GSE80178_down_test <- tT %>%
  filter(logFC < -1 & adj.P.Val < 0.3)
#------------------------------------------------------------------------

#-----CHECKING IF THE RELAXED CRITERIA RETURN ANY ROWS-------------------
print(nrow(GSE80178_up_test))
print(nrow(GSE80178_down_test))
#------------------------------------------------------------------------

#----------PROCEEDING TO WRITE A FILE WITH THEM AND STORE IT IN THE 
# WORKING DIRECTORY------------------------------------------------------
#Proceed to write them in the files that will be saved into the folder.
write.table(GSE80178_up_test, file="GSE80178_UR_1", row.names=F, sep="\t")
write.table(GSE80178_down_test, file="GSE80178_DR_1", row.names=F, sep="\t")
#------------------------------------------------------------------------