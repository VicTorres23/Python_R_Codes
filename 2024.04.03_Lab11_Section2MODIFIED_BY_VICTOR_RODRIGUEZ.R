
## install bioconductor
## Uncomment the line below to install/ Follow the instructions provided in bioconductor.org
## to install the bioconductor package
source("https://bioconductor.org/biocLite.R")

library(GEOquery)
BiocManager::install("GEOquery")
#BiocManager::install("biocLite")
#library(biocLite)
#biocLite("biocLite")


## Preprocessing and normalization of Affymetrix expression data

## 1.set directory which has .CEL files (raw data)
## Change to your working directory
#----------------------MODIFIED WORKING DIRECTORY------------------------
setwd("C:/Users/vemma/Documents/Master in Bioinformatics/Spring_2024_Semester/Introduction to Bioinformatics II/Lab 11") 
#--------------------------------------------------------------------------
#set this to YOUR working directory

## 2. Load the "library" that contains the Affymetrix microarray code we will
## Load the "library" that contains the Affymetrix microarray code
## we will want to use with the command

## Uncomment the below three lines to install the affy package (if you have not installed earlier)
 #if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
BiocManager::install("affy")

library(affy)

## read the CEL files (first command below) and then 
## summarize and normalize with "rma" introduced in the lecture class
## The provided input files (brain and liver cel) present in your directory will be loaded
affy.data = ReadAffy() 
eset.rma = rma(affy.data)

## preprocessing with MAS5

eset.mas5 = mas5(affy.data)

## The variable "eset.rma" contains normalized expression
## values for all probesets, along with other information.
## Let's continue by getting the expression matrix(probesets/genes in rows,
## chips in columns)
exprSet.nologs = exprs(eset.mas5)

## List the column (chip) names
colnames(exprSet.nologs)

## Rename the column names if we want
colnames(exprSet.nologs) = c("brain.1","brain.2","fetal.brain.1",
                             "fetal.brain.2",
                             "fetal.liver.1","fetal.liver.2",
                             "liver.1","liver.2")
# check column names
colnames(exprSet.nologs)
# At this time let's log-transform the expression values to get a more normal distribution.
# We have to remember that we have to do this when we calculate ratios.
# Logarithms can use any base, but base 2 is easiest when transforming ratios,
# since transformed 2-fold ratios up or doewn will be +1 or -1.
# As a result, we will do all logs with base 2 to keep thing simplest.

exprSet = log(exprSet.nologs, 2)

############################################## End of getting expresion data

## Statistical Analysis 
## 1. t-test

## t test for (brain.1,brain.2) vs. (fetal.brain.1,fetal.brain.2)
##            (liver.1,liver.2) vs. (fetal.liver.1,fetal.liver.2)
## raw p value
t.p.b =apply(exprSet,1,function(x)t.test(x[1:2],x[3:4])$p.value)
t.p.l = apply(exprSet,1,function(x)t.test(x[5:6],x[7:8])$p.value)

write.csv(t.p.b,"number_1_file.csv", row.names = TRUE)

## 2. ANOVA test: ratio of variance: F statistic
y = c(0,0,1,1)
a.p.b =apply(exprSet,1,function(x)summary(aov(y~x[1:4]))[[1]][1,5])
a.p.l = apply(exprSet,1,function(x)summary(aov(y~x[5:8]))[[1]][1,5])

## 3. multiple correction Benjamini-Hochberg
## Follow documentation to try other methods such as bonferroni
#-------ADDING MORE METHODS AND APPLYING THEM TO THE ANOVA P-VALUES--------
padj.b = p.adjust(t.p.b,method="bonferroni")
padj.l = p.adjust(t.p.l,method="bonferroni")

padj.bBH <- p.adjust(t.p.b, method="BH")
padj.lBH <- p.adjust(t.p.l,method="BH")

ANOVApadj.bBonf <- p.adjust(a.p.b,method="bonferroni")
ANOVApadj.lBonf <- p.adjust(a.p.l,method="bonferroni")

ANOVApadj.bBH <- p.adjust(a.p.b,method="BH")
ANOVApadj.lBH <- p.adjust(a.p.l,method="BH")
#--------------------------------------------------------------------------

#-WE DO A SCATTERPLOT COMPARING THE RAW P-VALUES WITH THE ADJUSTED P-VALUES 
# OF EACH METHOD-----------------------------------------------------------
plot(t.p.b, padj.b, main="Brain Tissue T-Test Original P-Value vs Bonferroni Adjusted P-Value")
plot(t.p.b, padj.bBH, main="Brain Tissue T-Test Original P-Value vs Benjamin & Hochberg Adjusted P-Value")
plot(t.p.l, padj.l, main="Liver Tissue T-Test Original P-Value vs Bonferroni Adjusted P-Value")
plot(t.p.l, padj.lBH, main="Liver Tissue T-Test Original P-Value vs Benjamin & Hochberg Adjusted P-Value")
plot(a.p.b, ANOVApadj.bBonf, main="Brain Tissue ANOVA Original P-Value vs Bonferroni Adjusted P-Value")
plot(a.p.l, ANOVApadj.lBonf, main="Liver Tissue ANOVA Original P-Value vs Bonferroni Adjusted P-Value")
plot(a.p.b, ANOVApadj.bBH, main="Brain Tissue ANOVA Original P-Value vs Benjamin & Hochberg Adjusted P-Value")
plot(a.p.l, ANOVApadj.lBH, main="Liver Tissue ANOVA Original P-Value vs Benjamin & Hochberg Adjusted P-Value")

plot(a.p.l, ANOVApadj.lBH, type = "l", col = "blue", lwd = 2, xlab = "X-axis", ylab = "Y-axis", main = "Line Plot Example")

#---------------------------------------------------------------------------
