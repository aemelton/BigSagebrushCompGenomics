### AE Melton, 2023
# Edited from S Buerki sliding window script to characterize TE type content

library(tidyverse)
library(rstatix)
library(ggpubr)

##~~
# Set working directory
setwd("~/Dropbox/BSU_Research/ReferenceGenomics_2022/EDTA/")

##~~
# Set Individual ID
ID <- "UTT2"

##~~
# Vector w/ all scaffolds
scaffID <- paste0("Scaffold_", seq(from=1, to=9, by=1))

##~~
# Load Chr length file and subset to ID 
# Csv w/ Chr info
chrL <- read.csv("Sagebrush_Genome_Lengths.csv")
# Subset chrL to ID
chrLInd <- subset(chrL, chrL$Individual == ID)

##~~
# Load GFF file
gffFile <- "UTT2/UTT2_Consensus_ONELINER.fasta.mod.EDTA.TEanno.gff3" # "IDT3/Use_This_One.fasta.mod.EDTA.TEanno.gff3"
gff <- read.csv(gffFile, sep = '\t', header=T)
head(gff)
# colnames(gff) <- c("scaffID", "source", "sequence_ontology", "start", "end", "score", "strand", "phase", "attributes")


##~~
# Declare variable for moving window
#Set size of sliding window
slidwind <- 5000000
##~~
# Loop to infer length of TE within each chunk
# This is done across all scaffolds
OUT <- NULL
for(j in 1:length(scaffID)){
  #Subset to target scaff
  scaff <- scaffID[j]
  
  #Create start points for moving window associated to target scaff
  starts <- seq(1, chrLInd$Length_bp[which(chrLInd$Scaffold == scaff)]-as.numeric(slidwind), by = as.numeric(slidwind))
  
  # Set number of iterations in loop based on starts
  n <- length(starts)
  
  # Create an empty matrix to store TE values
  chunkData <- matrix(ncol=4, nrow=n)
  colnames(chunkData) <- c("Sample",
                           "Chromosome",
                           "chromStart",
                           "chromEnd")
                           
  chunkData <- as.data.frame(chunkData)
  chunkData$Sample <- rep(ID, nrow(chunkData))
  chunkData$Chromosome <- rep(scaff, nrow(chunkData))
  chunkData$chromStart <- starts
  chunkData$chromEnd <- chunkData$chromStart+(as.numeric(slidwind)-1)
  
  #Subset GFF to scaff
  gffscaff <- subset(gff, gff$seqid == scaff)
  
  #Loop to identify # TE (= bp) per sliding window
  pb <- txtProgressBar(min = 0, max = nrow(chunkData), style = 3)
  
  TEcounts <- matrix(ncol=14, nrow=n)
  for(i in 1:nrow(chunkData)){
    #Subset gff file to match sliding window parameters
    tmp <- gffscaff[which(gffscaff$start >= as.numeric(chunkData$chromStart[i]) & gffscaff$end <= as.numeric(chunkData$chromEnd[i])),]
    #Count occurrence of each type of TE
    tmp.df <- as.data.frame(x = table(tmp$sequence_ontology))
    rownames(tmp.df) <- tmp.df$Var1
    tmp.data <- data.frame(t(tmp.df["Freq"]))
    
    if (!is.null(tmp.data$CACTA_TIR_transposon) == T) {
      TEcounts[i,1] <- tmp.data$CACTA_TIR_transposon
    } else {TEcounts[i,1] <- 0 }
    
    if (!is.null(tmp.data$Copia_LTR_retrotransposon) == T) {
      TEcounts[i,2] <- tmp.data$Copia_LTR_retrotransposon
    } else {TEcounts[i,2] <- 0 }
    
    if (!is.null(tmp.data$Gypsy_LTR_retrotransposon) == T) {
      TEcounts[i,3] <- tmp.data$Gypsy_LTR_retrotransposon
    } else {TEcounts[i,3] <- 0 }
    
    if (!is.null(tmp.data$hAT_TIR_transposon) == T) {
      TEcounts[i,4] <- tmp.data$hAT_TIR_transposon
    } else {TEcounts[i,4] <- 0 }
    
    if (!is.null(tmp.data$helitron) == T) {
      TEcounts[i,5] <- tmp.data$helitron
    } else {TEcounts[i,5] <- 0 }
    
    if (!is.null(tmp.data$LINE_element) == T) {
      TEcounts[i,6] <- tmp.data$LINE_element
    } else {TEcounts[i,6] <- 0 }
    
    if (!is.null(tmp.data$long_terminal_repeat) == T) {
      TEcounts[i,7] <- tmp.data$long_terminal_repeat
    } else {TEcounts[i,7] <- 0 }
    
    if (!is.null(tmp.data$LTR_retrotransposon) == T) {
      TEcounts[i,8] <- tmp.data$LTR_retrotransposon
    } else {TEcounts[i,8] <- 0 }
    
    if (!is.null(tmp.data$Mutator_TIR_transposon) == T) {
      TEcounts[i,9] <- tmp.data$Mutator_TIR_transposon
    } else {TEcounts[i,9] <- 0 }
    
    if (!is.null(tmp.data$PIF_Harbinger_TIR_transposon) == T) {
      TEcounts[i,10] <- tmp.data$PIF_Harbinger_TIR_transposon
    } else {TEcounts[i,10] <- 0 }
    
    if (!is.null(tmp.data$polinton) == T) {
      TEcounts[i,11] <- tmp.data$polinton
    } else {TEcounts[i,11] <- 0 }
    
    if (!is.null(tmp.data$repeat_region) == T) {
      TEcounts[i,12] <- tmp.data$repeat_region
    } else {TEcounts[i,12] <- 0 }
    
    if (!is.null(tmp.data$target_site_duplication) == T) {
      TEcounts[i,13] <- tmp.data$target_site_duplication
    } else {TEcounts[i,13] <- 0 }
    
    if (!is.null(tmp.data$Tc1_Mariner_TIR_transposon) == T) {
      TEcounts[i,14] <- tmp.data$Tc1_Mariner_TIR_transposon
    } else {TEcounts[i,14] <- 0 }
    
    # TEcounts <- tmp.data
    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  OUT <- rbind(OUT, TEcounts)
  #OUT
  # OUT <- rbind(OUT, TEcounts)
}
colnames(OUT) <- c("CACTA_TIR_transposon",
                        "Copia_LTR_retrotransposon",
                        "Gypsy_LTR_retrotransposon",
                        "hAT_TIR_transposon",
                        "helitron",
                        "LINE_element",
                        "long_terminal_repeat",
                        "LTR_retrotransposon",
                        "Mutator_TIR_transposon",
                        "PIF_Harbinger_TIR_transposon",
                        "polinton",
                        "repeat_region",
                        "target_site_duplication",
                        "Tc1_Mariner_TIR_transposon")
#Export file as csv
#OUT <- cbind(chunkData[1:4], TEcounts)
write.csv(OUT, file = paste0(ID, "/", ID, "_TEs_per_slidingwindow_5Mil.csv"), row.names = F)


##################
IDT3 <- read.csv("IDT3/IDT3_TEs_per_slidingwindow_5Mil.csv")
UTT2 <- read.csv("UTT2/UTT2_TEs_per_slidingwindow_5Mil.csv")
IDT2 <- read.csv("IDT2/IDT2_TEs_per_slidingwindow_5Mil.csv")

IDT3$Pop <- rep("IDT3", nrow(IDT3))
UTT2$Pop <- rep("UTT2", nrow(UTT2))
IDT2$Pop <- rep("IDT2", nrow(IDT2))

AllDat <- rbind(IDT3, UTT2, IDT2)
head(AllDat)

mydata.long <- AllDat %>%
  pivot_longer(-Pop, names_to = "variables", values_to = "value")
mydata.long %>% sample_n(6)

aov.out <- aov(formula = mydata.long$value ~ mydata.long$Pop)
summary(aov.out)

# Df    Sum Sq Mean Sq F value Pr(>F)
# mydata.long$Pop     2 1.599e+05   79935   1.202  0.301
# Residuals       82345 5.476e+09   66498               


stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ Pop) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

stat.test.df <- data.frame(stat.test)
print(stat.test.df, n = 50)
stat.test.df[order(stat.test.df$p.adj, decreasing = F),]

stat.test.sig <- subset(x = stat.test, stat.test$p.adj <= 0.05)
write.csv(x = stat.test.sig, file = "SigTTestPerType_5Mil.csv", row.names = F)

summary(stat.test.sig$statistic)
quantile(x = stat.test.sig$statistic, probs = c(0.05, 0.95))

BigDiff <- stat.test.sig[which(stat.test.sig$statistic <= -18.28397 | stat.test.sig$statistic >= 15.01195),]
BigDiff
write.csv(x = BigDiff, file = "BigDiff_SigTTestPerType_5Mil.csv", row.names = F)


# Per chunk comp between two genomes
IDT3 <- read.csv("IDT3/IDT3_TEs_per_slidingwindow_5Mil.csv")
UTT2 <- read.csv("UTT2/UTT2_TEs_per_slidingwindow_5Mil.csv")
IDT2 <- read.csv("IDT2/IDT2_TEs_per_slidingwindow_5Mil.csv")

out1 <- NULL
for (i in 1:nrow(IDT3)) {
  out1[i] <- t.test(IDT3[i,],  UTT2[i,])$p.value
}
out1[which(out1 <= 0.05)]
summary(out1)
length(out1)

out2 <- NULL
for (i in 1:nrow(IDT3)) {
  out2[i] <- t.test(IDT3[i,],  IDT2[i,])$p.value
}
out2[which(out2 <= 0.05)]
summary(out2)
length(out2)

out3 <- NULL
for (i in 1:nrow(IDT2)) {
  out3[i] <- t.test(IDT2[i,],  UTT2[i,])$p.value
}
out3[which(out3 <= 0.05)]
summary(out3)
length(out3)

######

