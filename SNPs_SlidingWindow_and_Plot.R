### Borrow from Vanilla Genomics by SB and PE

###
# R code to infer SNP density along moving window and produce a map

# 1. Load packages
library(seqinr)
library(SNPRelate)
library(RCircos)

# 2. Create a data.frame with SNP information
ID <- "UTT2"

genofile <- snpgdsOpen(filename = paste0("SNPs/", ID, ".ccm.gds"))
SNPsInfo <- data.frame(Chr = read.gdsn(index.gdsn(genofile, "snp.chromosome")), Position = read.gdsn(index.gdsn(genofile, "snp.position")), ID = read.gdsn(index.gdsn(genofile, "snp.id")))
SNPsInfo$Chr <- gsub(pattern = "Scaffold_", replacement = "", x = SNPsInfo$Chr)

#Extract chr lengths
##~~
# Vector w/ all scaffolds
scaffID <- paste0("Scaffold_", seq(from=1, to=9, by=1))

##~~
# Load Chr length file and subset to ID 
# Csv w/ Chr info
chrL <- read.csv("EDTA/Sagebrush_Genome_Lengths.csv")
# Subset chrL to ID
chrLength <- subset(chrL, chrL$Individual == ID)
  
# chrLength <- as.matrix(summary(chrFasta)[,1])
colnames(chrLength) <- c("Individual", "Chromosome", "ChrLength")
chrLength$Chromosome <- gsub(pattern = "Scaffold_", replacement = "", x = chrLength$Chromosome)
  
print("Infer SNP density along sliding window")
SNPslidingdat <- NULL
for(j in 1:nrow(chrLength)){
  chr <- chrLength$Chromosome[j]
  slidwind <- 5000000
    
  #Establish foundation for inferring nbr of SNP per sliding window length
  # Start by creating a vector to set the position on the sliding window
  starts <- seq(1, as.numeric(chrLength[match(chr, chrLength$Chromosome),3])-as.numeric(slidwind), by = as.numeric(slidwind))
    
  #Create table with SNP position and ID for target chr
  chrDat <- subset(SNPsInfo, SNPsInfo$Chr == chr)
    
  # Set number of iterations in loop based on starts
  n <- length(starts)
    
  # Create an empty matrix to store GC values
  chunkData <- matrix(ncol=5, nrow=n)
  colnames(chunkData) <- c("Chromosome", "chromStart", "chromEnd", "SNP_nbr", "SNP_perc")
  chunkData <- as.data.frame(chunkData)
  chunkData$Chromosome <- rep(chr, nrow(chunkData))
  chunkData$chromStart <- starts
  chunkData$chromEnd <- as.numeric(chunkData$chromStart)+(as.numeric(slidwind)-1)
    
  ###
  # Sliding window: SNP nbr and percentages
  ###
  print(paste0("Processing chromsome: ", chr))
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for(i in 1:n){
    #Where are we in the loop
    #print(paste("Iteration:", i, "from", n, sep=' '))
    #Infer nbr of SNP in sliding window
    chunkData$SNP_nbr[i] <- length(which(chrDat$Position %in% seq(from = chunkData$chromStart[i], to = chunkData$chromEnd[i], by = 1)))
    #Infer SNP percentage
    chunkData$SNP_perc[i] <- round(as.numeric(chunkData$SNP_nbr[i])/as.numeric(slidwind), 2)
      
    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
    
  SNPslidingdat <- rbind(chunkData, SNPslidingdat)
}
  
#Load data for map
#Chromosome data
Van <- data.frame(Chromosome = chrLength$Chromosome, ChromStart = rep(0, nrow(chrLength)), ChromEnd = as.vector(chrLength$ChrLength), Band = rep("x", nrow(chrLength)), Stain = rep("gneg", nrow(chrLength)))
Van$Chromosome <- as.character(Van$Chromosome)
Van$ChromEnd <- as.numeric(Van$ChromEnd)
  
#Plot data
out.file <- paste0("SNPs/SNP_density_map_5mil", ID, ".pdf");
pdf(file=out.file, height=9, width=9, compress=TRUE)
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
rcircos.cyto$ChrColor <- rep("black", nrow(chrLength))
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
RCircos.Heatmap.Plot(SNPslidingdat, data.col = 4, track.num = 2, side = "in", min.value = 0, max.value = max(SNPslidingdat$SNP_nbr))
text(x=0, y=1, "#SNPs", cex=.8)
  
dev.off()

