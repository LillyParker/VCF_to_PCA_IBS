#Install SNPRelate
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("gdsfmt")
BiocManager::install("SNPRelate")

#Load required libraries 
library(gdsfmt)
library('SNPRelate')
library(vcfR)
library(reshape2)

#read VCF file in and convert to GDS format for SNPRelate 
vcf.fn <- ("swift_bed_DP4_Rm50_RmCoys_RmFoxes_Rm349_miss75_mac2_thin100K.recode.vcf")
snpgdsVCF2GDS(vcf.fn, "Swift2023v7.gds", method="biallelic.only")

#Confirm correct specs in GDS file 
snpgdsSummary("Swift2023v7.gds")

#Load/open GDS file
Swift2023v7.gds <- snpgdsOpen("Swift2023v7.gds")
Swift2023v7.gds

#PCA analysis to confirm species ID - investigate outliers

Swift2023v7_pca <- snpgdsPCA(Swift2023v7.gds, num.thread=2, autosome.only = F)
Swift2023v7_pca

Swiftv7_pc.percent <- Swift2023v7_pca$varprop*100
head(round(Swiftv7_pc.percent, 2))

#Convert to dataframe and plot

Swift2023v7_pca_table <- data.frame(sample.id=Swift2023v7_pca$sample.id,
                                    EV1 = Swift2023v7_pca$eigenvect[,1],
                                    EV2=Swift2023v7_pca$eigenvect[,2],
                                    stringsAsFactors = FALSE)


plot(Swift2023v7_pca_table$EV2, Swift2023v7_pca_table$EV1, xlab = "Eigenvector 2", ylab="Eigenvector 1",main="PCA")

#Pairwise relatedness calculation (identity-by-descent, IBD), using MLE
sample.id <- read.gdsn(index.gdsn(Swift2023v7.gds, "sample.id"))
snp.id <- read.gdsn(index.gdsn(Swift2023v7.gds, "snp.id"))

set.seed(100)
ibdSwift7 <- snpgdsIBDMLE(Swift2023v7.gds, sample.id=sample.id, snp.id=snp.id, num.thread=2, autosome.only=FALSE)
ibdSwift7.coeff <- snpgdsIBDSelection(ibdSwift7)
ibdSwift7.coeff

#write to file 
write.table(ibdSwift7.coeff, file="ibd_Swift7.txt")

# Perform bootstrapping analysis
# Install kinship utils from https://github.com/campanam/kinshipUtils

source("kinshipUtils/kinshipUtils.R")
SwiftBoot <- bootstrap.kinship(Swift2023v7.gds,ibdmethod = "MLE", mlemethod = "EM", resample = 100, autosome.only=FALSE)

Matrix <- write.kinship.matrix(SwiftBoot, cifile="SwiftBootCIv7.csv")
Mean <- write.kinship.matrix(SwiftBoot, meanfile = "SwiftBootMeanv7.csv")

#Estimate pairwise identity-by-state (IBS) to identify recaptures 

ibsSwift7 <- snpgdsIBS(Swift2023v7.gds, sample.id=sample.id, snp.id=snp.id, num.thread=2, autosome.only=FALSE)

#Reshape output 
ibsSwift7_TABLE <- subset(melt[(ibsSwift7$ibs), value!=0])

#Write to table
write.table(ibsSwift7_TABLE, file="IBS_Swift7.csv")

