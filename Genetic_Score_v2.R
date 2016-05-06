#Read in EBI GWAS
setwd("C:/Users/Jeremy/Documents/Fourth Year/Data Science")
EBI <- read.delim("gwas_catalog_v1.0.1.tsv", header=TRUE,sep="\t",fileEncoding="windows-1252")

#read in files from Genome A
setwd("C:/Users/Jeremy/Documents/Fourth Year/Data Science/Imputed from Website/Genome_A_Imputed")
file_list <- list.files()

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    GenomeA <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    GenomeA<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  
}
#read in files from Genome B
setwd("C:/Users/Jeremy/Documents/Fourth Year/Data Science/Imputed from Website/Genome_B_Imputed")
file_list <- list.files()

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    GenomeB <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    GenomeB<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  
}
#read in files from Genome C
setwd("C:/Users/Jeremy/Documents/Fourth Year/Data Science/Imputed from Website/Genome_C_Imputed")
file_list <- list.files()

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    GenomeC <- read.table(file, header=TRUE, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    GenomeC<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  
}


#Alzheimers
#Step 1 - Store imputed SNPs for Alzheimers 
aSNP <-match(GenomeA, EBI[912:961,"SNPS"])
bSNP <-match(GenomeB, EBI[912:961,"SNPS"])
cSNP <-match(GenomeC, EBI[912:961,"SNPS"])
#Step 2 - find which SNPs for Alzheimers are common between the genomes
commonality <- combine(aSNP,bSNP,cSNP)
aSNP <- match(aSNP,commonality)
bSNP <- match(bSNP,commonality)
cSNP <- match(cSNP,commonality)
#Step 3 - calculate genetic score
riskAllele <- split(EBI[912:961, "STRONGEST SNP-RISK ALLELE"],"-")
numberOfAllelesA <- countCharacter(aSNP, riskAllele(1))
numberOfAllelesB <- countCharacter(bSNP, riskAllele(1))
numberOfAllelesC <- countCharacter(cSNP, riskAllele(1))
ScoreA <- sum(numberOfAllelesA * EBI[912:961,"OR or BETA"])
ScoreB <- sum(numberOfAllelesB * EBI[912:961,"OR or BETA"])
ScoreC <- sum(numberOfAllelesC * EBI[912:961,"OR or BETA"])

#Type 2 Diabetes
#Step 1 - Store imputed SNPs for Type 2 Diabetes
aSNP <-match(GenomeA, EBI[6733-6794,"SNPS"])
bSNP <-match(GenomeB, EBI[6733-6794,"SNPS"])
cSNP <-match(GenomeC, EBI[6733-6794,"SNPS"])
#Step 2 - find which SNPs for Type 2 Diabetes are common between the genomes
commonality <- combine(aSNP,bSNP,cSNP)
aSNP <- match(aSNP,commonality)
bSNP <- match(bSNP,commonality)
cSNP <- match(cSNP,commonality)
#Step 3 - calculate genetic score
riskAllele <- split(EBI[6733-6794, "STRONGEST SNP-RISK ALLELE"],"-")
numberOfAllelesA <- countCharacter(aSNP, riskAllele(1))
numberOfAllelesB <- countCharacter(bSNP, riskAllele(1))
numberOfAllelesC <- countCharacter(cSNP, riskAllele(1))
ScoreA <- sum(numberOfAllelesA * EBI[6733-6794,"OR or BETA"])
ScoreB <- sum(numberOfAllelesB * EBI[6733-6794,"OR or BETA"])
ScoreC <- sum(numberOfAllelesC * EBI[6733-6794,"OR or BETA"])

#Hypertension
#Step 1 - Store imputed SNPs for Hypertension
aSNP <-match(GenomeA, EBI[16117:16141,"SNPS"])
bSNP <-match(GenomeB, EBI[16117:16141,"SNPS"])
cSNP <-match(GenomeC, EBI[16117:16141,"SNPS"])
#Step 2 - find which SNPs for Hypertension are common between the genomes
commonality <- combine(aSNP,bSNP,cSNP)
aSNP <- match(aSNP,commonality)
bSNP <- match(bSNP,commonality)
cSNP <- match(cSNP,commonality)
#Step 3 - calculate genetic score
riskAllele <- split(EBI[16117:16141, "STRONGEST SNP-RISK ALLELE"],"-")
numberOfAllelesA <- countCharacter(aSNP, riskAllele(1))
numberOfAllelesB <- countCharacter(bSNP, riskAllele(1))
numberOfAllelesC <- countCharacter(cSNP, riskAllele(1))
ScoreA <- sum(numberOfAllelesA * EBI[16117:16141,"OR or BETA"])
ScoreB <- sum(numberOfAllelesB * EBI[16117:16141,"OR or BETA"])
ScoreC <- sum(numberOfAllelesC * EBI[16117:16141,"OR or BETA"])

#Coronary Artery Disease
#Step 1 - Store imputed SNPs for Coronary Artery Disease
aSNP <-match(GenomeA, EBI[19296:19325,"SNPS"])
bSNP <-match(GenomeB, EBI[19296:19325,"SNPS"])
cSNP <-match(GenomeC, EBI[19296:19325,"SNPS"])
#Step 2 - find which SNPs for Coronary Artery Disease are common between the genomes
commonality <- combine(aSNP,bSNP,cSNP)
aSNP <- match(aSNP,commonality)
bSNP <- match(bSNP,commonality)
cSNP <- match(cSNP,commonality)
#Step 3 - calculate genetic score
riskAllele <- split(EBI[19296:19325, "STRONGEST SNP-RISK ALLELE"],"-")
numberOfAllelesA <- countCharacter(aSNP, riskAllele(1))
numberOfAllelesB <- countCharacter(bSNP, riskAllele(1))
numberOfAllelesC <- countCharacter(cSNP, riskAllele(1))
ScoreA <- sum(numberOfAllelesA * EBI[19296:19325,"OR or BETA"])
ScoreB <- sum(numberOfAllelesB * EBI[19296:19325,"OR or BETA"])
ScoreC <- sum(numberOfAllelesC * EBI[19296:19325,"OR or BETA"])
