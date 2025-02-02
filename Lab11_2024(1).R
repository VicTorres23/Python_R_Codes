##############################################################

### Don't forget to use the right directory when running this codes

#setwd("C:/Users/sgayathrina/OneDrive - The University of Texas at El Paso/1_Sharan_PhD/Coursework/2024_SEM4_Spring/TA2024/lab13_2024/")
BiocManager::install("minfi")
library(minfi)
library(GEOquery)

#### Extracting important data from the GEODataSet
 gset <- getGEO("GSE68777", GSEMatrix =TRUE, AnnotGPL=FALSE)
 if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
 gset <- gset[[idx]]
 pData(gset)
 head(sampleNames(gset))
 
####################################
gset <- getGEO("GSE68777", GSEMatrix =TRUE, AnnotGPL=FALSE)

## if you're unable to download the files using getGEOSuppFiles, i'll provide and alternative
## I'll go through the alternative in class 
getGEOSuppFiles("GSE68777")
untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
head(list.files("GSE68777/idat", pattern = "idat"))

### This code below unzip the files in the folder GSE68777/idat
file_list <- list.files("GSE68777/idat")

for (file in file_list) {
  # check if the file is a gz archive file
  if (grepl(".gz$", file)) {
    # extract the contents of the tar archive file
    gunzip(paste0("GSE68777/idat", "/", file))
  }
}


rgSet <-  read.metharray.exp("GSE68777/idat")
rgSet

pData(rgSet) ## This returns the phenotype data for the gene expression data GSE68777
head(sampleNames(rgSet))

#####
#options(timeout=600)  # Set timeout to 600 seconds (10 minutes) if reaching time out limitation
#geoMat <- getGEO("GSE68777")

geoMat <- getGEO(filepath="C:/Users/vemma/Documents/Master in Bioinformatics/Spring_2024_Semester/Introduction to Bioinformatics II/Lab 13/GSE68777/GSE68777_RAW")


pD.all <- pData(geoMat[[1]])
View(pD.all)
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
head(pD)
names(pD)[c(3,4)] <- c("group", "sex")
pD$group <- sub("^diagnosis:", "", pD$group)
pD$sex <- sub("^Sex:", "", pD$sex)

sampleNames(rgSet) <- sub(".*_5", "5", sampleNames(rgSet))

# rownames(pD) <- pD$title
# pD <- pD[sampleNames(rgSet),]
# pData(rgSet) <- pD
rgSet

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

grSet <- preprocessQuantile(rgSet) ##performs quantile normalization on the input rgSet, think why?
grSet
granges(grSet) 

## used to retrieve genomic locations of CpG sites on the genome

##Beta values are a common way to represent DNA methylation data 
##from Illumina methylation arrays, and are defined as the ratio of methylated 
##probe intensity to the sum of methylated and unmethylated probe intensities.

beta = getBeta(grSet) ## to retrieve the beta values
beta[1:3, 1:3]
head(getIslandStatus(grSet)) ## to determine the CpG island status of the probes on the methylation array
grSet
pheno = pD$group

dmp <- dmpFinder(beta, pheno = pheno, type = "categorical") ##to identify differentially methylated positions (DMPs) between groups of samples based on their methylation beta values
# sum(dmp$qval<.05)
#dmp[c(1,6,285:289),]
idx <- which(dmp$qval<0.03)
sig.beta = beta[idx,]

d <- dist(as.matrix((sig.beta)))
hc <- hclust(d)
par(mar = rep(2, 4))
plot(hc)
dmp[idx,]


idx = which(dmp$qval<0.03)

new.data = data.frame(sig.beta[,idx], 
                      sig.beta[,-idx])


heatmap(as.matrix(new.data))


### End here
##############################################################
