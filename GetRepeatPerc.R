### AE Melton, 2022
# Count soft-masked repeat nucleotides in an assembly. Code assumes single line per scaffold.

#### Load libraries ####
library(stringr)

#### Load data ####
# Load the assembly (R won't like really big assemblies, so going scaffold by scaffold may be best)
genome <- readLines(con = "FinalScaffs_pseudochrom.fa")
head(genome, 1) # check just first header
length(genome) # Get number of lines in assembly
scaffolds <- genome[!grepl(">", genome)]
summary(scaffolds)
rm(genome) # Free up some memory
#

#### Count soft-masked and total characters per scaffold ####
soft.masked <- str_count(scaffolds, "[a-z]")
not.masked <- str_count(scaffolds, "[A-Z]")
rm(scaffolds)
total.char <- soft.masked + not.masked

#### Get percent masked for samples ####
IDT3.perc.masked <- soft.masked/total.char

# Load the assembly (R won't like really big assemblies, so going scaffold by scaffold may be best)
genome <- readLines(con = "UTT2_consensus.fa") # readLines(con = "FinalScaffs_pseudochrom.fa")
head(genome, 1) # check just first header
length(genome) # Get number of lines in assembly
scaffolds <- genome[!grepl(">", genome)]
summary(scaffolds)
rm(genome) # Free up some memory
#

# Count soft-masked and total characters per scaffold
soft.masked <- str_count(scaffolds, "[a-z]")
not.masked <- str_count(scaffolds, "[A-Z]")
rm(scaffolds)
total.char <- soft.masked + not.masked
UTT2.perc.masked <- soft.masked/total.char

# Load the assembly (R won't like really big assemblies, so going scaffold by scaffold may be best)
genome <- readLines(con = "IDT2_consensus.fa") # readLines(con = "FinalScaffs_pseudochrom.fa")
head(genome, 1) # check just first header
length(genome) # Get number of lines in assembly
scaffolds <- genome[!grepl(">", genome)]
summary(scaffolds)
rm(genome) # Free up some memory
#

# Count soft-masked and total characters per scaffold
soft.masked <- str_count(scaffolds, "[a-z]")
not.masked <- str_count(scaffolds, "[A-Z]")
rm(scaffolds)
total.char <- soft.masked + not.masked
IDT2.perc.masked <- soft.masked/total.char
#

#### Let's make a table of the results ####
Samples <- c("IDT3", "UTT2", "IDT2")
Scaffold.Names <- c("Scaffold 1", "Scaffold 2", "Scaffold 3",
                    "Scaffold 4", "Scaffold 5", "Scaffold 6",
                    "Scaffold 7", "Scaffold 8", "Scaffold 9")
masked.data <- data.frame(IDT3.perc.masked, UTT2.perc.masked, IDT2.perc.masked)

####