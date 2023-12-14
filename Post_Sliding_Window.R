### AE Melton, 2023
# Post-sliding window analysis
setwd("~/Dropbox/BSU_Research/ReferenceGenomics_2022/EDTA/")

# Libraries
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

# read in some files
IDT3.window <- read.csv("IDT3/IDT3_AllTE_slidingwindow_WITHPERCENT_5mil_AddWindow.csv")# Had to add one sliding window to end of Chrom3, # add window at end of 5 for 5mil base windows
UTT2.window <- read.csv("UTT2/UTT2_AllTE_slidingwindow_WITHPERCENT_5mil_AddWindow.csv")# add window at end of 5 for 5mil base windows
IDT2.window <- read.csv("IDT2/IDT2_AllTE_slidingwindow_WITHPERCENT_5mil.csv")

# Check the number of windows
table(IDT3.window$Chromosome)
table(UTT2.window$Chromosome)
table(IDT2.window$Chromosome)

# Add a pop tag for later analyses
IDT3.window$Pop <- rep("IDT3", nrow(IDT3.window))
UTT2.window$Pop <- rep("UTT2", nrow(UTT2.window))
IDT2.window$Pop <- rep("IDT2", nrow(IDT2.window))

# Let's do a density plot and compare the densities of TEs across the genomes
GlassHouses <- rbind(IDT3.window, UTT2.window, IDT2.window)
head(GlassHouses)


# Density plots
pdf("DensityPlot.pdf")
ggplot(data = GlassHouses,
       aes(x = TE,
           group = Pop,
           fill = Pop)) +
        geom_density(adjust = 1.5, alpha = 0.4) +
        scale_fill_manual(values = c("green", "red", "blue")) # theme_ipsum()
dev.off()

# T.tests with no added window for IDT3
t.test(IDT3.window$TE, IDT2.window$TE)
### Helitron
# data:  IDT3.window$TE and IDT2.window$TE
# t = -7.2377, df = 3914.4, p-value = 5.465e-13
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2464.543 -1413.932
# sample estimates:
#   mean of x mean of y 
# 19095.57  21034.81 

# # Helitron, 5mil window
# Welch Two Sample t-test
# 
# data:  IDT3.window$TE and IDT2.window$TE
# t = -6.5088, df = 1552.5, p-value = 1.02e-10
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#         -6370.414 -3419.962
# sample estimates:
#         mean of x mean of y 
# 47698.95  52594.14

### All
# data:  IDT3.window$TE and IDT2.window$TE
# t = 3.3778, df = 3876.9, p-value = 0.000738
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#         6069.021 22861.114
# sample estimates:
#         mean of x mean of y 
# 2504422   2489956 

# data:  IDT3.window$TE and IDT2.window$TE
# t = 13.502, df = 1412.9, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#         171031.9 229177.0
# sample estimates:
#         mean of x mean of y 
# 6427519   6227415 

t.test(IDT3.window$TE, UTT2.window$TE)
### Helitron
# data:  IDT3.window$TE and UTT2.window$TE
# t = -16.326, df = 3907.8, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -4949.025 -3887.793
# sample estimates:
#   mean of x mean of y 
# 19095.57  23513.98 

# # Helitron, 5mil window
# Welch Two Sample t-test
# 
# data:  IDT3.window$TE and UTT2.window$TE
# t = -14.469, df = 1544.4, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#         -12532.991  -9540.523
# sample estimates:
#         mean of x mean of y 
# 47698.95  58735.71 

# All
# data:  IDT3.window$TE and UTT2.window$TE
# t = -0.32616, df = 3890.8, p-value = 0.7443
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#         -9877.222  7059.603
# sample estimates:
#         mean of x mean of y 
# 2504422   2505830 

# data:  IDT3.window$TE and UTT2.window$TE
# t = 9.9376, df = 1561.8, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#         135338.6 201903.9
# sample estimates:
#         mean of x mean of y 
# 6427519   6258898 

t.test(IDT2.window$TE, UTT2.window$TE)
### Helitron
# data:  IDT2.window$TE and UTT2.window$TE
# t = -9.0094, df = 3918.5, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -3018.674 -1939.670
# sample estimates:
#   mean of x mean of y 
# 21034.81  23513.98 

# Helitron, 5mil window
# Welch Two Sample t-test
# 
# data:  IDT2.window$TE and UTT2.window$TE
# t = -7.774, df = 1559.9, p-value = 1.371e-14
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#         -7691.168 -4591.969
# sample estimates:
#         mean of x mean of y 
# 52594.14  58735.71 

# All
# data:  IDT2.window$TE and UTT2.window$TE
# t = -3.8816, df = 3918.6, p-value = 0.0001055
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#         -23891.568  -7856.186
# sample estimates:
#         mean of x mean of y 
# 2489956   2505830 

# data:  IDT2.window$TE and UTT2.window$TE
# t = -2.1393, df = 1420.8, p-value = 0.03258
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#         -60351.576  -2614.777
# sample estimates:
#         mean of x mean of y 
# 6227415   6258898 

# Calculate TE richness differences across all sliding windows (including the fake one added to IDT3)
idt2.TE.rich.diff <- IDT3.window$TE - IDT2.window$TE
utt2.TE.rich.diff <- IDT3.window$TE - UTT2.window$TE
idt2.utt2.TE.rich.diff <- IDT2.window$TE - UTT2.window$TE

# Generate a plot of the differences

pdf("~/Dropbox/Manuscripts/BSU/1_Genome_Comparisons/Figures/AllDiff.pdf")
par(mfrow = c(3, 1))

plot(utt2.TE.rich.diff,
     ylab = "TE richness difference (per sliding window)",
     xlab = "Sliding window position",
     pch = 19,
     cex = 0.5)
mtext(text = "IDT3 vs. UTT2")
abline(h = 0, col = "red")

plot(idt2.TE.rich.diff,
     ylab = "TE richness difference (per sliding window)",
     xlab = "Sliding window position",
     ylim = c(-200000, 530000),
     pch = 19,
     cex = 0.5)
mtext(text = "IDT3 vs. IDT2")
abline(h = 0, col = "red")

plot(idt2.utt2.TE.rich.diff,
     ylab = "TE richness difference (per sliding window)",
     xlab = "Sliding window position",
     ylim = c(-280000, 300000),
     pch = 19,
     cex = 0.5)
mtext(text = "IDT2 vs. UTT2")
abline(h = 0, col = "red")
dev.off()

# Get a summary of the differences
summary(idt2.TE.rich.diff)
summary(utt2.TE.rich.diff)
summary(idt2.utt2.TE.rich.diff)

### Helitron
# IDT3 - UTT2
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -57306   -7124   -2218   -2199    2843   40830 

# UTT2 - IDT2
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -26194   -1079    2758    2479    6225   48908 

# IDT3 - IDT2
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -46000.0  -5162.0    329.0    280.6   5328.0  44659.0 

### All
# > summary(idt2.TE.rich.diff)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -2686589   -22781    15267    14465    54972   227998 
# 
# > summary(utt2.TE.rich.diff)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -2762675        0        0    -1409        0        0 
# 
# > summary(idt2.utt2.TE.rich.diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -227998  -55314  -15401  -15874   22690  173473

# Make a histogram of the differences
pdf("AllHist.pdf")
par(mfrow = c(3, 1))

hist(utt2.TE.rich.diff,
     main = NULL,
     xlim = c(-40000, 40000),
     xlab = "Difference in TE richness (per sliding window)")
mtext("A")

hist(idt2.TE.rich.diff,
     main = NULL,
     xlim = c(-40000, 40000),
     xlab = "Difference in TE richness (per sliding window)")
mtext("B")

hist(idt2.utt2.TE.rich.diff,
     main = NULL,
     xlim = c(-40000, 40000),
     xlab = "Difference in TE richness (per sliding window)")
mtext("C")

dev.off()

###
# Add in quantiles to identify TE enriched regions; extract the TEs, check types
quantile(idt2.TE.rich.diff, probs = c(0.05, 0.95))
idt2.big.diff.windows <- subset(cbind(IDT2.window, idt2.TE.rich.diff), idt2.TE.rich.diff <= -11868 & idt2.TE.rich.diff <= 6702)


quantile(utt2.TE.rich.diff, probs = c(0.05, 0.95))
utt2.big.diff.windows <- subset(cbind(UTT2.window, utt2.TE.rich.diff), utt2.TE.rich.diff <= -13826 & idt2.TE.rich.diff <= 4537)
sort(x = table(utt2.big.diff.windows$Chromosome), decreasing = T)

quantile(idt2.utt2.TE.rich.diff, probs = c(0.05, 0.95))

# and what kind of SNPs are around them?

# Compare composition of each sliding window - do all the windows have different TEs?
# If so, different genome shaping processes (hybridization? Allopoly? Homeologous chromosome mixing?)

# switching over to SvenLab computer for this bit (ssh and run from terminal below)
# IDT3.gff <- read.csv(file = "EDTA/IDT3/Use_This_One.fasta.mod.EDTA.TEanno.gff3",
#                      header = T,
#                      sep = "\t")
# head(IDT3.gff)
# 
# IDT3.gff.subset <- subset(IDT3.gff,
#                           seqid == "Scaffold_5")
# head(IDT3.gff.subset)
# little.IDT3.window <- subset(IDT3.window,
#                              Chromosome == "Scaffold_5")
# IDT3.Scaff5 <- NULL
# for (i in 1:nrow(IDT3.gff.subset)) {
#         for (j in 1:nrow(little.IDT3.window)) {
#                 if (IDT3.gff.subset$start[i] >= little.IDT3.window$chromStart[j] & IDT3.gff.subset$end[i] <= little.IDT3.window$chromEnd[j]) {
#                         IDT3.Scaff5 <- rbind(subset.test, IDT3.gff.subset[i,])
#                 }
#         }
# }
# head(IDT3.Scaff5)
# 
# # UTT2
# UTT2.gff <- read.csv(file = "UTT2/UTT2_Consensus_ONELINER.fasta.mod.EDTA.TEanno.gff3",
#                      header = T,
#                      sep = "\t")
# head(UTT2.gff)
# 
# UTT2.gff.subset <- subset(UTT2.gff,
#                           seqid == "Scaffold_5")
# head(UTT2.gff.subset)
# little.UTT2.window <- subset(UTT2.window,
#                              Chromosome == "Scaffold_5")
# UTT2.Scaff5 <- NULL
# for (i in 1:nrow(UTT2.gff.subset)) {
#         for (j in 1:nrow(little.UTT2.window)) {
#                 if (UTT2.gff.subset$start[i] >= little.UTT2.window$chromStart[j] & UTT2.gff.subset$end[i] <= little.UTT2.window$chromEnd[j]) {
#                         UTT2.Scaff5 <- rbind(UTT2.Scaff5, UTT2.gff.subset[i,])
#                 }
#         }
# }
# head(UTT2.Scaff5)
# 
