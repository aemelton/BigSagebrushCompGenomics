### AE Melton, 2023
# 

# Load libraries
library(stringr)
setwd("~/Dropbox/BSU_Research/ReferenceGenomics_2022/EDTA/")

# Read in gff3. gff3 files have been modified to remove metadata at start of file and change spaces to tabs in headers
# TE.anno.IDT3 <- read.csv(file = "IDT3/IDT3_Consensus_ONELINER.fasta.mod.EDTA.TEanno.gff3", header = T, sep = "\t")
# head(TE.anno.IDT3)

# 
# helitron.IDT3 <- grep(pattern = "helitron", x = TE.anno.IDT3, value = T)
system("cat IDT3/Use_This_One.fasta.mod.EDTA.TEanno.gff3 | grep \"helitron\" > IDT3/IDT3_Helitron.tsv")
IDT3.Helitron <- read.csv(file = "IDT3/IDT3_Helitron.tsv", header = F, sep = "\t")
head(IDT3.Helitron)
nrow(IDT3.Helitron)

# Break down the atributes column (column 9) to get the sequence ontology
SeqOG <- str_extract(string = IDT3.Helitron[,9], pattern = "Name=TE_\\d*")
SeqOG <- gsub(pattern = "Name=TE_", replacement = "", x = SeqOG)
head(SeqOG)
length(unique(SeqOG))

# Extract sequences from fasta
setwd("IDT3/")
IDT3.Helitron <- read.csv("IDT3_Helitron.tsv", header = F, sep = "\t")
head(IDT3.Helitron)
SeqOG <- str_extract(string = IDT3.Helitron[,9], pattern = "Name=TE_\\d*")
SeqOG <- gsub(pattern = "Name=TE_", replacement = "", x = SeqOG)
SeqOG
length(unique(SeqOG))

SeqOG_Headers <- NULL
for (i in 1:length(unique(SeqOG))) {
  SeqOG_Headers[i] <- gsub(pattern = "\\d*", replacement = paste0(">TE_", unique(SeqOG)[i], "#DNA/Helitron"), x = unique(SeqOG)[i])
}
SeqOG_Headers
SeqOG_Headers <- SeqOG_Headers[!is.na(SeqOG_Headers)]

GrepSeq <- NULL
for (i in 1:length(SeqOG_Headers)) {
GrepSeq <- paste0("grep -A 1 \"", SeqOG_Headers[i], "\" Use_This_One.fasta.mod.EDTA.TElib.fa >> IDT3_Unique_Helitron_Seq.fa")
system(GrepSeq)
}

#### Count occurrences of unique TEs and SO numbers
# 
setwd("IDT3/") # setwd to sample folder
IDT3.Helitron <- read.csv("IDT3_Helitron.tsv", header = F, sep = "\t")
head(IDT3.Helitron)
SeqOG <- str_extract(string = IDT3.Helitron[,9], pattern = "Name=TE_\\d*")
sort(x = table(SeqOG), decreasing = T) # Sorted table of count data for each TE#

table(IDT3.Helitron$V1)

df <- grep(pattern = "TE_00037659", x = IDT3.Helitron$V9) # Pull out specific TE#s
head(df)

subset.Helitron <- IDT3.Helitron[df,]
unique(subset.Helitron$V1) # Which chromosomes are the TEs on?

#### 