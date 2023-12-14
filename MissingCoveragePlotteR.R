### AE Melton, 2023
# Look at bam coverages
# All of this was written with regards to the Helitron gff

# Read in some data. Made with genomeCoverageBed -ibam UTT2.sam.bam -bga | awk '$4 == 0' > UTT2.intervals.with.0.cov.bedg
bam.cvg.file <- read.csv("BAMs/UTT2.intervals.with.0.cov.bedg", sep = "\t", header = F)
head(bam.cvg.file)

# # Let's subset first....
# bam.cvg.file.Scaff2 <- subset(bam.cvg.file, bam.cvg.file$V1 == "Scaffold_2")
# nrow(bam.cvg.file.Scaff2)
# 
# # Next, plot some stuff
# plot(x = 0, y = 0, xlim = c(0, 513745672)) # Make the length of the scaffold the x lim
# for (i in 1:nrow(bam.cvg.file.Scaff2)) {
#   segments(x0 = bam.cvg.file.Scaff2$V2[i], x1 = bam.cvg.file.Scaff2$V3[i], y0 = 0, col = "red", lwd = 1)
# }

# Let's loop over the table to add the gaps together per scaffold and then use table()
sums.out <- NULL
# for (i in 1:9) {
#   scaffID <- paste0("Scaffold_", i)
#   scaff.subset <- subset(bam.cvg.file, bam.cvg.file$V1 == scaffID)
  for (j in 1:nrow(bam.cvg.file)) {
    tmp.out <- bam.cvg.file$V3[j] - bam.cvg.file$V2[j]
    sums.out <- rbind(sums.out, tmp.out)
  }
bam.cvg.file <- cbind(bam.cvg.file, sums.out)
# }
head(sums.out)
hist(bam.cvg.file$sums.out)
range(bam.cvg.file$sums.out) # 1 to 47752 for UTT2
bam.cvg.file[which(bam.cvg.file$sums.out == max(bam.cvg.file$sums.out)),]

# make hist of no read coverage gaps
hist(bam.cvg.file$sums.out)

# Check location against TE gff
UTT2.gff <- read.csv(file = "EDTA/UTT2/UTT2_Consensus_ONELINER.fasta.mod.EDTA.TEanno.gff3",
                     header = T,
                     sep = "\t")
head(UTT2.gff)

UTT2.gff.subset <- subset(UTT2.gff,
                          seqid == "Scaffold_5")
UTT2.gff.subset

UTT2.gff.subset2 <- subset(UTT2.gff.subset,
                           start > 85350000 & start < 85450000)
UTT2.gff.subset2
range(UTT2.gff.subset2$end - UTT2.gff.subset2$start) # Check the sizes of TEs; do any ~match sizes of gaps in mapping?
quantile((UTT2.gff.subset2$end - UTT2.gff.subset2$start), c(0.05, 0.95)) # 

UTT2.gff.subset2[which((UTT2.gff.subset2$end - UTT2.gff.subset2$start) == max(range(UTT2.gff.subset2$end - UTT2.gff.subset2$start))),]

# Check reference
IDT3.gff <- read.csv(file = "EDTA/IDT3/Use_This_One.fasta.mod.EDTA.TEanno.gff3",
                     header = T,
                     sep = "\t")
head(IDT3.gff)

IDT3.gff.subset <- subset(IDT3.gff,
                          seqid == "Scaffold_5")
IDT3.gff.subset

IDT3.gff.subset2 <- subset(IDT3.gff.subset,
                           start > 85380000 & start < 85430000)

bam.cvg.file[which(bam.cvg.file$sums.out == max(bam.cvg.file$sums.out)),] # UTT2 mapping
head(UTT2.gff.subset2) # UTT2 gff
head(IDT3.gff.subset2) # IDT3 gff

range(IDT3.gff.subset2$end - IDT3.gff.subset2$start) # Check the sizes of TEs; do any ~match sizes of gaps in mapping?
quantile((IDT3.gff.subset2$end - IDT3.gff.subset2$start), c(0.05, 0.95)) # 

IDT3.gff.subset2[which((IDT3.gff.subset2$end - IDT3.gff.subset2$start) == max(range(IDT3.gff.subset2$end - IDT3.gff.subset2$start))),]



