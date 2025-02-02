#############STEP 1##############

#############STEP 1.3##############

## list all the vignettes in the RforProteomics package
#vignette(package = "RforProteomics")
## Open the vignette called RforProteomics
vignette("RforProteomics", package = "RforProteomics")
## or just
vignette("RforProteomics")

#############STEP 1.4##############
## only first time you install Bioconductor packages
BiocManager::install("RforProteomics", force = TRUE)

if (!requireNamespace("BiocManager", quietly=TRUE, force = TRUE))
  install.packages("BiocManager")
## else
library("BiocManager")
BiocManager::install("RforProteomics", dependencies = TRUE, update = TRUE)

#############STEP 1.5##############
BiocManager::install("RforProteomics", dependencies = TRUE)

library("RforProteomics")

#############STEP 1.6##############
## gets the vignette source
rfile <- system.file("doc/RforProteomics.R",
                     package = "RforProteomics")
rfile

#############STEP 1.7##############
install.packages("reshape2")

library("RColorBrewer") ## Color palettes
library("ggplot2")  ## Convenient and nice plotting
library("reshape2") ## Flexibly reshape data

#############STEP 2##############

#############STEP 2.1.1##############
library("mzR") ## the software package
library("msdata") ## the data package
BiocManager::install("msdata")

## below, we extract the releavant example file
## from the local 'msdata' installation
filepath <- system.file("microtofq", package = "msdata")
file <- list.files(filepath, pattern="MM14.mzML",
                   full.names=TRUE, recursive = TRUE)
## creates a commection to the mzML file
mz <- openMSfile(file)
## demonstraction of data access
basename(fileName(mz))

runInfo(mz)

BiocManager::install("affy")

instrumentInfo(mz)

## once finished, it is good to explicitely
## close the connection
close(mz)

#############STEP 2.1.2##############
file <- system.file("mzid", "Tandem.mzid.gz", package="msdata")
mzid <- openIDfile(file)
mzid

softwareInfo(mzid)

enzymes(mzid)

names(psms(mzid))

head(psms(mzid))[, 1:13]

#############STEP 2.2##############
library("mzID")
mzids <- list.files(system.file('extdata', package = 'mzID'),
                    pattern = '*.mzid', full.names = TRUE)
mzids

id <- mzID(mzids[1])

id

ids <- mzID(mzids[1:2])
ids

fid <- flatten(id)
names(fid)

dim(fid)

#############STEP 3##############
library("affy")
library("MSnbase")
## uses a simple dummy test included in the package
mzXML <- dir(system.file(package="MSnbase",dir="extdata"),
             full.name=TRUE,
             pattern="mzXML$")
basename(mzXML)

## reads the raw data into and MSnExp instance
raw <- readMSData(mzXML, verbose = FALSE, centroided = TRUE)
raw

## Extract a single spectrum
raw[[3]]

plot(raw, full = TRUE)
plot(raw[[3]], full = TRUE, reporters = iTRAQ4)

#############STEP 4##############

#############STEP 4.1##############
BiocManager::install("rpx")

library("rpx")

px1 <- PXDataset("PXD000001")
px1

pxfiles(px1)

## Downloading the mzTab data
mztab <- pxget(px1, "F063721.dat-mztab.txt")
mztab

## Load mzTab peptide data
qnt <- readMzTabData(mztab, what = "PEP", version = "0.9")

sampleNames(qnt) <- reporterNames(TMT6)
head(exprs(qnt))

## remove missing values
qnt <- filterNA(qnt)
processingData(qnt)

## combine into proteins
## - using the 'accession' feature meta data
## - sum the peptide intensities
protqnt <- combineFeatures(qnt,
                           groupBy = fData(qnt)$accession,
                           method = sum)

cls <- brewer.pal(5, "Set1")
matplot(t(tail(exprs(protqnt), n = 5)), type = "b",
        lty = 1, col = cls,
        ylab = "Protein intensity (summed peptides)",
        xlab = "TMT reporters")
legend("topright", tail(featureNames(protqnt), n=5),
       lty = 1, bty = "n", cex = .8, col = cls)

qntS <- normalise(qnt, "sum")
qntV <- normalise(qntS, "vsn")

qntV2 <- normalise(qnt, "vsn")

acc <- c("P00489", "P00924",
         "P02769", "P62894",
         "ECA")

idx <- sapply(acc, grep, fData(qnt)$accession)
idx2 <- sapply(idx, head, 3)
small <- qntS[unlist(idx2), ]

idx3 <- sapply(idx, head, 10)
medium <- qntV[unlist(idx3), ]

m <- exprs(medium)
colnames(m) <- c("126", "127", "128",
                 "129", "130", "131")
rownames(m) <- fData(medium)$accession
rownames(m)[grep("CYC", rownames(m))] <- "CYT"
rownames(m)[grep("ENO", rownames(m))] <- "ENO"
rownames(m)[grep("ALB", rownames(m))] <- "BSA"
rownames(m)[grep("PYGM", rownames(m))] <- "PHO"
rownames(m)[grep("ECA", rownames(m))] <- "Background"

cls <- c(brewer.pal(length(unique(rownames(m)))-1, "Set1"),
         "grey")
names(cls) <- unique(rownames(m))
wbcol <- colorRampPalette(c("white", "darkblue"))(256)  

heatmap(m, col = wbcol, RowSideColors=cls[rownames(m)])

dfr <- data.frame(exprs(small),
                  Protein = as.character(fData(small)$accession),
                  Feature = featureNames(small),
                  stringsAsFactors = FALSE)
colnames(dfr) <- c("126", "127", "128", "129", "130", "131",
                   "Protein", "Feature")
dfr$Protein[dfr$Protein == "sp|P00924|ENO1_YEAST"] <- "ENO"
dfr$Protein[dfr$Protein == "sp|P62894|CYC_BOVIN"]  <- "CYT"
dfr$Protein[dfr$Protein == "sp|P02769|ALBU_BOVIN"] <- "BSA"
dfr$Protein[dfr$Protein == "sp|P00489|PYGM_RABIT"] <- "PHO"
dfr$Protein[grep("ECA", dfr$Protein)] <- "Background"
dfr2 <- melt(dfr)

ggplot(aes(x = variable, y = value, colour = Protein),
       data = dfr2) +
  geom_point() +
  geom_line(aes(group=as.factor(Feature)), alpha = 0.5) +
  facet_grid(. ~ Protein) + theme(legend.position="none") +
  labs(x = "Reporters", y = "Normalised intensity")

#############STEP 4.3##############
mzxml <- pxget(px1, "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzXML")

rawms <- readMSData(mzxml, centroided = TRUE, verbose = FALSE)
qntms <- quantify(rawms, reporters = TMT7, method = "max")
qntms

d <- data.frame(Signal = rowSums(exprs(qntms)[, 1:6]),
                Incomplete = exprs(qntms)[, 7])
d <- log(d)
cls <- rep("#00000050", nrow(qnt))
pch <- rep(1, nrow(qnt))
cls[grep("P02769", fData(qnt)$accession)] <- "gold4" ## BSA
cls[grep("P00924", fData(qnt)$accession)] <- "dodgerblue" ## ENO
cls[grep("P62894", fData(qnt)$accession)] <- "springgreen4" ## CYT
cls[grep("P00489", fData(qnt)$accession)] <- "darkorchid2" ## PHO
pch[grep("P02769", fData(qnt)$accession)] <- 19
pch[grep("P00924", fData(qnt)$accession)] <- 19
pch[grep("P62894", fData(qnt)$accession)] <- 19
pch[grep("P00489", fData(qnt)$accession)] <- 19

mzp <- plotMzDelta(rawms, reporters = TMT6, verbose = FALSE) + ggtitle("")
mzp

plot(Signal ~ Incomplete, data = d,
     xlab = expression(Incomplete~dissociation),
     ylab = expression(Sum~of~reporters~intensities),
     pch = 19,
     col = "#4582B380")
grid()

abline(0, 1, lty = "dotted")
abline(lm(Signal ~ Incomplete, data = d), col = "darkblue")

MAplot(qnt[, c(4, 2)], cex = .9, col = cls, pch = pch, show.statistics = FALSE)

#############STEP 4.4##############
## load packages
library("MALDIquant")
install.packages("MALDIquantForeign")

library("MALDIquantForeign")
## getting test data
datapath <-
  file.path(system.file("Examples",
                        package = "readBrukerFlexData"),
            "2010_05_19_Gibb_C8_A1")
dir(datapath)

sA1 <- importBrukerFlex(datapath, verbose=FALSE)
# in the following we use only the first spectrum
s <- sA1[[1]]

summary(mass(s))
summary(intensity(s))

head(as.matrix(s))

plot(s)

###PreProcessing#####

## sqrt transform (for variance stabilization)
s2 <- transformIntensity(s, method="sqrt")
s2

s3 <- smoothIntensity(s2, method="MovingAverage", halfWindowSize=2)
s3

## baseline subtraction
s4 <- removeBaseline(s3, method="SNIP")
s4

#Peak Picking

## peak picking
p <- detectPeaks(s4)
length(p) # 181

peak.data <- as.matrix(p) # extract peak information

par(mfrow=c(2,3))
xl <- range(mass(s))
# use same xlim on all plots for better comparison
plot(s, sub="", main="1: raw", xlim=xl)
plot(s2, sub="", main="2: variance stabilisation", xlim=xl)

plot(s3, sub="", main="3: smoothing", xlim=xl)
plot(s4, sub="", main="4: base line correction", xlim=xl)

plot(s4, sub="", main="5: peak detection", xlim=xl)
points(p)

top20 <- intensity(p) %in% sort(intensity(p), decreasing=TRUE)[1:20]
labelPeaks(p, index=top20, underline=TRUE)
plot(p, sub="", main="6: peak plot", xlim=xl)

labelPeaks(p, index=top20, underline=TRUE)

#############STEP 4.5##############
BiocManager::install("BRAIN")

library(BRAIN)

atoms <- getAtomsFromSeq("SIVPSGASTGVHEALEMR")
unlist(atoms)

BiocManager::install("Rdisop")

library(Rdisop)
pepmol <- getMolecule(paste0(names(atoms),
                             unlist(atoms),
                             collapse = ""))
pepmol

BiocManager::install("OrgMassSpecR")

library(OrgMassSpecR)


data(itraqdata)

simplottest <-
  itraqdata[featureNames(itraqdata) %in% paste0("X", 46:47)]
sim <- SpectrumSimilarity(as(simplottest[[1]], "data.frame"),
                          as(simplottest[[2]], "data.frame"),
                          top.lab = "itraqdata[['X46']]",
                          bottom.lab = "itraqdata[['X47']]",
                          b = 25)

if (!file.exists("P00924.fasta"))
  eno <- download.file("http://www.uniprot.org/uniprot/P00924.fasta",
                       destfile = "P00924.fasta")

eno <- paste(readLines("P00924.fasta")[-1], collapse = "")
enopep <- Digest(eno, missed = 1)
nrow(enopep) ## 103

# sum(nchar(enopep$peptide) >= minlength) ## 68

pepcnt <- enopep[enopep[, 1] %in% exppep, ]
nrow(pepcnt) ## 13

BiocManager::install("cleaver")

library(cleaver)
cleave("LAAGKVEDSD", enzym = "trypsin")

## miss one cleavage position
cleave("LAAGKVEDSD", enzym = "trypsin", missedCleavages = 1)

## miss zero or one cleavage positions
cleave("LAAGKVEDSD", enzym = "trypsin", missedCleavages = 0:1)

## 15N incorporation rates from 0, 0.1, ..., 0.9, 0.95, 1
incrate <- c(seq(0, 0.9, 0.1), 0.95, 1)
inc <- lapply(incrate, function(inc)
  IsotopicDistributionN("YEVQGEVFTKPQLWP", inc))
#par(mfrow = c(4,3))
par(mfrow = c(1,1))
for (i in 1:length(inc))
  plot(inc[[i]][, c(1, 3)], xlim = c(1823, 1848), type = "h",
       main = paste0("15N incorporation at ", incrate[i]*100, "%"))

######################STEP 4.6#################################
BiocManager::install("isobar")

library("isobar")

## Welcome to isobar (v 1.40.0)
##    'openVignette("isobar")' and '?isobar' provide help on usage.
## 
## Attaching package: 'isobar'
## The following object is masked from 'package:xtable':
## 
##     sanitize
## The following object is masked from 'package:MSnbase':
## 
##     normalize
## The following object is masked from 'package:ProtGenerics':
## 
##     peptides
## The following object is masked from 'package:BiocGenerics':
## 
##     normalize
## The following object is masked from 'package:base':
## 
##     paste0
## Prepare the PXD000001 data for isobar analysis
.ions <- exprs(qnt)
.mass <- matrix(TMT6@mz, nrow(qnt), byrow=TRUE, ncol = 6)
colnames(.ions) <- colnames(.mass) <-
  reporterTagNames(new("TMT6plexSpectra"))
rownames(.ions) <- rownames(.mass) <-
  paste(fData(qnt)$accession, fData(qnt)$sequence, sep = ".")
pgtbl <- data.frame(spectrum = rownames(.ions),
                    peptide = fData(qnt)$sequence,
                    modif = ":",
                    start.pos = 1,
                    protein = fData(qnt)$accession,
                    accession = fData(qnt)$accession)
x <- new("TMT6plexSpectra", pgtbl, .ions, .mass)

##  data.frame columns OK
## Creating ProteinGroup ... done
featureData(x)$proteins <- as.character(fData(qnt)$accession)

x <- correctIsotopeImpurities(x) ## using identity matrix here

## LOG: isotopeImpurities.corrected: TRUE
x <- isobar::normalize(x, per.file = FALSE)

## LOG: is.normalized: TRUE
## LOG: normalization.multiplicative.factor channel 126: 0.8905
## LOG: normalization.multiplicative.factor channel 127: 0.9288
## LOG: normalization.multiplicative.factor channel 128: 1
## LOG: normalization.multiplicative.factor channel 129: 0.949
## LOG: normalization.multiplicative.factor channel 130: 0.8677
## LOG: normalization.multiplicative.factor channel 131: 0.8965
## spikes
spks <- c(protein.g(proteinGroup(x), "P00489"),
          protein.g(proteinGroup(x), "P00924"),
          protein.g(proteinGroup(x), "P02769"),
          protein.g(proteinGroup(x), "P62894"))

cls2 <- rep("#00000040", nrow(x))
pch2 <- rep(1, nrow(x))
cls2[grep("P02769", featureNames(x))] <- "gold4" ## BSA
cls2[grep("P00924", featureNames(x))] <- "dodgerblue" ## ENO
cls2[grep("P62894", featureNames(x))] <- "springgreen4" ## CYT
cls2[grep("P00489", featureNames(x))] <- "darkorchid2" ## PHO
pch2[grep("P02769", featureNames(x))] <- 19
pch2[grep("P00924", featureNames(x))] <- 19
pch2[grep("P62894", featureNames(x))] <- 19
pch2[grep("P00489", featureNames(x))] <- 19

nm <- NoiseModel(x)

## [1]   0.0731332 644.0157028   2.6730885
ib.background <- subsetIBSpectra(x, protein=spks,
                                 direction = "exclude")

## Creating ProteinGroup ... done
nm.background <- NoiseModel(ib.background)

## [1] 0.01346028 2.85118430 0.84630754
ib.spks <- subsetIBSpectra(x, protein = spks,
                           direction="include",
                           specificity="reporter-specific")

## Creating ProteinGroup ... done
nm.spks <- NoiseModel(ib.spks, one.to.one=FALSE, pool=TRUE)

## 4 proteins with more than 10 spectra, taking top 50.
## [1] 0.0000000001 5.8291408884 0.6609997984
ratios <- 10^estimateRatio(x, nm,
                           channel1="127", channel2="129",
                           protein = spks,
                           combine = FALSE)[, "lratio"]

res <- estimateRatio(x, nm,
                     channel1="127", channel2="129",
                     protein = unique(fData(x)$proteins),
                     combine = FALSE,
                     sign.level = 0.01)[, c(1, 2, 6, 8)]
res <- as.data.frame(res)
res$lratio <- -(res$lratio)

cls3 <- rep("#00000050", nrow(res))
pch3 <- rep(1, nrow(res))
cls3[grep("P02769", rownames(res))] <- "gold4" ## BSA
cls3[grep("P00924", rownames(res))] <- "dodgerblue" ## ENO
cls3[grep("P62894", rownames(res))] <- "springgreen4" ## CYT
cls3[grep("P00489", rownames(res))] <- "darkorchid2" ## PHO
pch3[grep("P02769", rownames(res))] <- 19
pch3[grep("P00924", rownames(res))] <- 19
pch3[grep("P62894", rownames(res))] <- 19
pch3[grep("P00489", rownames(res))] <- 19

rat.exp <- c(PHO = 2/2,
             ENO = 5/1,
             BSA = 2.5/10,
             CYT = 1/1)
maplot(x,
       noise.model = c(nm.background, nm.spks, nm),
       channel1="127", channel2="129",
       pch = 19, col = cls2,
       main = "Spectra MA plot")
abline(h = 1, lty = "dashed", col = "grey")

legend("topright",
       c("BSA", "ENO", "CYT", "PHO"),
       pch = 19, col = c("gold4", "dodgerblue",
                         "springgreen4", "darkorchid2"),
       bty = "n", cex = .7)

###########################STEP 4.7###############################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")

library(SummarizedExperiment)

data(msnset)

se <- as(msnset, "SummarizedExperiment")
se

ms <- as(se, "MSnSet")
ms

###########################STEP 4.8###############################
## open the synapter vignette
BiocManager::install("synapter")
library("synapter")

synapterGuide()

###########################STEP 5###############################
###########################STEP 5.3###############################
BiocManager::install("MSnID")
library("MSnID")
msnid <- MSnID(".")

PSMresults <- read.delim(system.file("extdata", "human_brain.txt",
                                     package="MSnID"),
                         stringsAsFactors=FALSE)
psms(msnid) <- PSMresults
show(msnid)

mzids <- system.file("extdata", "c_elegans.mzid.gz", package="MSnID")
msnid <- read_mzIDs(msnid, mzids)

show(msnid)

msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
prop.table(table(msnid$numIrregCleavages))

pepCleav <- unique(psms(msnid)[,c("numMissCleavages", "isDecoy", "peptide")])
pepCleav <- as.data.frame(table(pepCleav[,c("numMissCleavages", "isDecoy")]))
library("ggplot2")
ggplot(pepCleav, aes(x=numMissCleavages, y=Freq, fill=isDecoy)) +
  geom_bar(stat='identity', position='dodge') +
  ggtitle("Number of Missed Cleavages")

msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))

filtObj <- MSnIDFilter(msnid)
filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=5.0)
filtObj$msmsScore <- list(comparison=">", threshold=8.0)
show(filtObj)

evaluate_filter(msnid, filtObj, level="PSM")

filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01,
                                method="Grid", level="peptide",
                                n.iter=500)
show(filtObj.grid)

filtObj.nm <- optimize_filter(filtObj.grid, msnid, fdr.max=0.01,
                              method="Nelder-Mead", level="peptide",
                              n.iter=500)
show(filtObj.nm)

evaluate_filter(msnid, filtObj, level="peptide")

evaluate_filter(msnid, filtObj.grid, level="peptide")

evaluate_filter(msnid, filtObj.nm, level="peptide")

msnid <- apply_filter(msnid, filtObj.nm)
show(msnid)

msnid <- apply_filter(msnid, "!grepl('Contaminant',accession)")
show(msnid)

psm.df <- psms(msnid)
psm.dt <- as(msnid, "data.table")



peps <- MSnID::peptides(msnid)
head(peps)

prots <- accessions(msnid)
head(prots)

prots <- proteins(msnid) # may be more intuitive then accessions
head(prots)

msnset <- as(msnid, "MSnSet")
library("MSnbase")
head(fData(msnset))

head(exprs(msnset))

msnset <- combineFeatures(msnset,
                          fData(msnset)$accession,
                          redundancy.handler="unique",
                          fun="sum",
                          cv=FALSE)

head(fData(msnset))

head(exprs(msnset))

#########################STEP 7##################################
getHpa(id, "hpaSubcellularLoc")

################################STEP 8################################## #other packages
BiocManager::install("RforProteomics")
library("RforProteomics")
pp <- proteomicsPackages()
View(pp)