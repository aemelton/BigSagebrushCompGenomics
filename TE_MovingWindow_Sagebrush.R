## Sven Buerki, Associate Professor, Boise State University, 2023
# Modified by AEM, 2023

# 1. Load packages
library(seqinr)
library(RCircos)

##~~
# Set working directory
setwd("~/Dropbox/BSU_Research/ReferenceGenomics_2022/EDTA/")

##~~
# Set Individual ID
ID <- "UTT2"

##~~
# Vector w/ all scaffolds
scaffID <- paste0("Chr", seq(from=1, to=9, by=1)) # paste0("Scaffold_", seq(from=1, to=9, by=1))

##~~
# Load Chr length file and subset to ID 
# Csv w/ Chr info
chrL <- read.csv("Sagebrush_Genome_Lengths.csv")
# Subset chrL to ID
chrLInd <- subset(chrL, chrL$Individual == ID)
for (i in 1:length(chrLInd$Scaffold)) { # number of scaffold names
  chrLInd$Scaffold <- gsub(pattern = paste0("Scaffold_", i), replacement = paste0("Chr", i), x = chrLInd$Scaffold)
}
chrLInd

##~~
# Load GFF file
gffFile <- paste0(ID, "/", ID, "_Consensus_ONELINER.fasta.mod.EDTA.TEanno.gff3") # "IDT3/Use_This_One.fasta.mod.EDTA.TEanno.gff3" # 
gff <- read.csv(gffFile, sep = '\t', header=T)
head(gff)
colnames(gff) <- c("seqid", "source", "sequence_ontology", "start", "end", "score", "strand", "phase", "attributes")

# Simplify the scaffold names
for (i in 1:length(unique(gff$seqid))) { # number of scaffold names
  gff$seqid <- gsub(pattern = paste0("Scaffold_", i), replacement = paste0("Chr", i), x = gff$seqid)
}
head(gff)
tail(gff)

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
  chunkData <- matrix(ncol=5, nrow=n)
  colnames(chunkData) <- c("Chromosome", "chromStart", "chromEnd", "TE", "TE_perc")
  chunkData <- as.data.frame(chunkData)
  
  chunkData$Chromosome <- rep(scaff, nrow(chunkData))
  chunkData$chromStart <- starts
  chunkData$chromEnd <- chunkData$chromStart+(as.numeric(slidwind)-1)
  
  #Subset GFF to scaff
  gffscaff <- subset(gff, gff$seqid == scaff)
  
  #Loop to identify # TE (= bp) per sliding window
  pb <- txtProgressBar(min = 0, max = nrow(chunkData), style = 3)
  for(i in 1:nrow(chunkData)){
    #Subset gff file to match sliding window parameters
    tmp <- gffscaff[which(gffscaff$start >= as.numeric(chunkData$chromStart[i]) & gffscaff$end <= as.numeric(chunkData$chromEnd[i])),]
    #Infer length of TE regions
    chunkData$TE[i] <- sum(tmp$end-tmp$start)
    #Infer SNP percentage
    chunkData$TE_perc[i] <- round(as.numeric(chunkData$TE[i])/as.numeric(slidwind), 4)*100
    
    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  #OUT
  OUT <- rbind(OUT, chunkData)
}
head(OUT)

#Export file as csv
write.csv(OUT, file = paste0(ID, "/", ID, "_AllTE_slidingwindow_WITHPERCENT_5mil.csv"), row.names = F)

#Chromosome data
Van <- data.frame(Chromosome = unique(OUT$Chromosome), ChromStart = rep(0, nrow(chrLInd)), ChromEnd = as.vector(chrLInd$Length_bp), Band = rep("x", nrow(chrLInd)), Stain = rep("gneg", nrow(chrLInd)))
Van$Chromosome <- as.character(Van$Chromosome)
Van$ChromEnd <- as.numeric(Van$ChromEnd)

#Plot data
out.file <- paste0(ID, "/", ID, "_AllTE_density_map_5mil", ".pdf");
pdf(file=out.file, height=10, width=10, compress=TRUE) # height=9, width=9, compress=TRUE
chr.exclude <- NULL
cyto.info <- Van
tracks.inside <- 4
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$track.background <- NULL
rcircos.params$hist.color <- "brown"
rcircos.cyto <- RCircos.Get.Plot.Ideogram()
rcircos.cyto$BandColor <- rcircos.cyto$ChrColor
rcircos.cyto$ChrColor <- rep("black", chrLInd)
RCircos.Reset.Plot.Ideogram(rcircos.cyto)
rcircos.position <- RCircos.Get.Plot.Positions()
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.List.Plot.Parameters()
RCircos.Set.Plot.Area()
par(mai=c(0.25, 0.25, 0.25, 0.25))
plot.new()
plot.window(c(-2.5,2.5), c(-2.5, 2.5))
RCircos.Draw.Chromosome.Ideogram()
RCircos.Chromosome.Ideogram.Plot(tick.interval = 20)

# Number of SNPs
RCircos.Heatmap.Plot(OUT,
                     data.col = 4,
                     track.num = 2,
                     side = "in",
                     min.value = 0,
                     max.value = max(OUT$TE))
text(x=0, y=1, "# TE", cex=.8)

dev.off()

#Simple plot of TE in moving window at chr level
# pdf("UTT2/UTT2_Chr_slidingwindow.pdf")
# for(i in 1:length(scaffID)){
#   tmp <- subset(OUT, OUT$Chromosome == scaffID[i])
#   plot(tmp$chromStart, tmp$TE, type="l")
#   
# }
# dev.off()
