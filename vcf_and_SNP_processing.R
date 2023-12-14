### AE Melton, 2022
# 


setwd("PATH/SNPs/")

#### Libraries ####
# 
# Load libraries
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("LEA")

library(adegenet)
library(LEA)
library(poppr)
library(SNPRelate)
library(vcfR)

#### Load data ####
# read in VCF file
vcf_file <- "IDT3.IndelGap.BiAllel.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )

# Query vcf file for various metdata
# queryMETA(vcf, element = 'DP')
# queryMETA(vcf, element = 'GT')
vcf

# ***** Object of Class vcfR *****
# 3 samples
# 9 CHROMs
# 40,508,404 variants
# Object size: 16599.3 Mb
# 0 percent missing data
# *****        *****         *****

### IDT3
# ***** Object of Class vcfR *****
# 1 samples
# 9 CHROMs
# 7,299,304 variants
# Object size: 3063.2 Mb
# 0 percent missing data
# *****        *****         *****

### UTT2
# ***** Object of Class vcfR *****
#         1 samples
# 9 CHROMs
# 20,228,009 variants
# Object size: 6427.3 Mb
# 0 percent missing data
# *****        *****         *****

### IDT2
# ***** Object of Class vcfR *****
# 1 samples
# 9 CHROMs
# 24,781,876 variants
# Object size: 8052.2 Mb
# 0 percent missing data
#  *****        *****         *****

#### Process vcf file and filter by coverage ####
# mask by read depth and write a cleaned vcf file
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)
rownames(dp) <- 1:nrow(dp)
tail(dp)
heatmap.bp(dp[1:1000,], rlabels = FALSE)

# 
quants <- apply(dp, MARGIN=2, quantile, probs=c(0.05, 0.95), na.rm = TRUE)
dp2 <- sweep(dp, MARGIN=2, FUN = "-", quants[1,])
dp[dp2 < 0] <- NA

dp2 <- sweep(dp, MARGIN=2, FUN = "-", quants[2,])
dp[dp2 > 0] <- NA

dp[dp < 2] <- NA

vcf@gt[,-1][ is.na(dp) == TRUE ] <- NA
vcf

# ***** Object of Class vcfR *****
# 3 samples
# 9 CHROMs
# 40,508,404 variants
# Object size: 16367.2 Mb
# 18 percent missing data
# *****        *****         **********

### IDT3
# ***** Object of Class vcfR *****
# 1 samples
# 9 CHROMs
# 7,299,304 variants
# Object size: 3020.4 Mb
# 9.774 percent missing data
# *****        *****         *****

### UTT2
# ***** Object of Class vcfR *****
#         1 samples
# 9 CHROMs
# 20,228,009 variants
# Object size: 6358.5 Mb
# 8.103 percent missing data
# *****        *****         *****

### IDT2
# ***** Object of Class vcfR *****
# 1 samples
# 9 CHROMs
# 24,781,876 variants
# Object size: 8052.2 Mb
# 8.423 percent missing data
#  *****        *****         *****
        
heatmap.bp(dp[1:1000,], rlabels = FALSE)

# Omit samples with too much missing data
myMiss <- apply(dp, MARGIN = 2, function(x){ sum( is.na(x) ) } )
myMiss <- myMiss / nrow(dp)
# length(vcf@gt[, c(TRUE, myMiss < 0.9)])

vcf@gt <- as.matrix(vcf@gt[, c(TRUE, myMiss < 0.9)])
vcf

# ***** Object of Class vcfR *****
# 3 samples
# 9 CHROMs
# 40,508,404 variants
# Object size: 16367.2 Mb
# 18 percent missing data
# *****        *****         *****

### IDT3
# ***** Object of Class vcfR *****
# 1 samples
# 9 CHROMs
# 7,299,304 variants
# Object size: 3020.4 Mb
# 9.774 percent missing data
# *****        *****         *****

### UTT2
# ***** Object of Class vcfR *****
#         1 samples
# 9 CHROMs
# 20,228,009 variants
# Object size: 6358.5 Mb
# 8.103 percent missing data
# *****        *****         *****

### IDT2
# ***** Object of Class vcfR *****
# 1 samples
# 9 CHROMs
# 24,781,876 variants
# Object size: 8052.2 Mb
# 8.423 percent missing data
#  *****        *****         *****

dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
heatmap.bp(dp[1:1000,], rlabels = FALSE)

# Omit variants missing from too many samples
myMiss <- apply(dp, MARGIN = 1, function(x){ sum( is.na(x) ) } )
myMiss <- myMiss / ncol(dp)
vcf <- vcf[myMiss < 0.5, ]
vcf

# ***** Object of Class vcfR *****
# 3 samples
# 9 CHROMs
# 36,518,066 variants
# Object size: 15004 Mb
# 11.9 percent missing data
# *****        *****         *****

### IDT3
# ***** Object of Class vcfR *****
# 1 samples
# 9 CHROMs
# 6,585,873 variants
# Object size: 2803.2 Mb
# 0 percent missing data
# *****        *****         *****

### UTT2
# ***** Object of Class vcfR *****
#         1 samples
# 9 CHROMs
# 18,588,960 variants
# Object size: 5932.8 Mb
# 0 percent missing data
# *****        *****         *****

### IDT2
# ***** Object of Class vcfR *****
# 1 samples
# 9 CHROMs
# 22,694,381 variants
# Object size: 7431 Mb
# 0 percent missing data
# *****        *****         *****

dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
heatmap.bp(dp[1:1000,], rlabels = FALSE)

write.vcf(x = vcf, file = "IDT3.cleaned.vcf.gz")

#### Assess SNPs for linkage disequilibrium ####
snpgdsVCF2GDS(vcf.fn = "IDT3.cleaned.vcf.gz", out.fn = "IDT3.ccm.gds") # vcf.fn = vcf_file
genofile <- openfn.gds("IDT3.ccm.gds")
snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2, autosome.only = F)
snpset.id <- unlist(unname(snpset))
head(snpset.id)
write.csv(x = snpset.id, file = "IDT3.snpset_LD_checked.csv", row.names = F) # Save the snps you want to put in PCA. Can use these to filter later if needed so you dont have to do the whoel LD thing again.

vcf.trim <- vcf[snpset.id,]
vcf.trim

# ***** Object of Class vcfR *****
# 3 samples
# 9 CHROMs
# 14,152 variants
# Object size: 8 Mb
# 2.37 percent missing data
# *****        *****         *****

### IDT3
# ***** Object of Class vcfR *****
# 1 samples
# 9 CHROMs
# 7,695 variants
# Object size: 3.7 Mb
# 0 percent missing data
# *****        *****         *****

### UTT2
# ***** Object of Class vcfR *****
#         1 samples
# 9 CHROMs
# 7,764 variants
# Object size: 3.4 Mb
# 0 percent missing data
# *****        *****         *****

### IDT2
# ***** Object of Class vcfR *****
# 1 samples
# 9 CHROMs
# 7,813 variants
# Object size: 3.4 Mb
# 0 percent missing data
# *****        *****         *****

write.vcf(x = vcf.trim, file = "IDT3.cleaned_subset.vcf.gz", APPEND = F) # Save the final cleaned and trimmed VCF file

#### Convert VCF to geno ####
system("gzip -dk IDT3.cleaned_subset.vcf.gz")

snpgdsVCF2GDS("IDT3.cleaned_subset.vcf", "IDT3.cleaned_subset.vcf.gds")
genofile <- snpgdsOpen(filename = "IDT3.cleaned_subset.vcf.gds")

#### Compare SNP calling success across read subsets of UTT2 data ####

# Compare SNP calling rates for subset data. How much data do we need to get ~ all the SNPs?

# UTT2 all variants
# 22,186,504

# read in VCF file
vcf_file <- "UTT2_0.05.bam.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 7,914,361

vcf_file <- "UTT2_0.1.bam.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 11,419,138

vcf_file <- "UTT2_0.25.bam.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 15,946,896

vcf_file <- "UTT2_0.5.bam.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 19,119,821

vcf_file <- "UTT2_0.75.bam.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 20,921,135

vcf_file <- "UTT2_0.9.bam.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 21,724,430

AllVariantCounts <- c(7914361, 11419138, 15946896, 19119821, 20921135, 21724430, 22186504)
AllProportionOriginal <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
AllReadDepth <- c(7.47, 13.94, 33.21, 65.14, 96.96, 116.05, 128.77)

# UTT2 Clean Call Biallic SNPs
# 20,228,009

# read in VCF file
vcf_file <- "UTT2_0.05.bam.IndelGap.BiAllel.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 7,312,958

vcf_file <- "UTT2_0.1.bam.IndelGap.BiAllel.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 10,488,355

vcf_file <- "UTT2_0.25.bam.IndelGap.BiAllel.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 14,587,544

vcf_file <- "UTT2_0.5.bam.IndelGap.BiAllel.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 17,454,263

vcf_file <- "UTT2_0.75.bam.IndelGap.BiAllel.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 19,082,028

vcf_file <- "UTT2_0.9.bam.IndelGap.BiAllel.calls.vcf.gz" # An unprocessed VCF
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf # 19,808,673

CleanVariantCounts <- c(7312958, 10488355, 14587544,17454263, 19082028, 19808673, 20228009)
CleanProportionOriginal <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
CleanReadDepth <- c(7.47, 13.94, 33.21, 65.14, 96.96, 116.05, 128.77)

# Plot by depth
plot(AllVariantCounts ~ AllReadDepth)

pdf(file = "~/Dropbox/Manuscripts/BSU/2_Genome_Comparisons/Figures/SNPsPerReadDepth.pdf")
plot(CleanVariantCounts ~ CleanReadDepth,
     ylab = "SNPs",
     xlab = "Average read depth",
     pch = 19,
     col = "blue",
     ylim = c(0,20500000),
     xlim = c(0, 130))
abline(v = 61.79, col = "red", cex = 5)
text("Read depth >60,\n75% of\nSNPs called", x = 80, y = 10000000, cex = 1, col = "red")
abline(v = 8.6, col = "black", cex = 5)
text("Read depth >8.5,\n50% of\nSNPs called", x = 27, y = 5000000, cex = 1, col = "black")
dev.off()

# nlme::(model = CleanVariantCounts ~ CleanReadDepth, data = df)

glm.out <- glm(CleanVariantCounts ~ CleanReadDepth)
summary(glm.out)

# 0.75 of all SNPs == 15171007, read depth == 65; 10114004 == 0.5 of SNPs
predict(object = glm.out, newdata = data.frame(CleanReadDepth = 8.6), type = "response")
##### ##### #####
