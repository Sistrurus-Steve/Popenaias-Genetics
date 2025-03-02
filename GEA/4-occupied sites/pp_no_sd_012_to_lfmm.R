# Author: Keaka Farleigh 
# Date: December 15th, 2020



##### !!! Things You Will Need !!! #####
# A .geno/.012 file 
# The original VCF used to create the 012 file

### Load your packages 
library(LEA)
library(lfmm)
library(vcfR)

  ########################
##### Read in Data #####
########################
setwd("/path/to/working/directory")

Geno_file <- 'pp_no_sd.012'
VCF_file <- 'pp_no_sd.recode.vcf'
Outputname <- 'pp_no_sd'

# Read in our .012 file and change the missing data value 
Dat = read.table(Geno_file)
Dat[Dat == '-1'] <- 9
Dat = Dat[,-1]
Dat = as.matrix(Dat)
write.geno(output.file = paste(Outputname,'_uncorrected', '.geno', sep = ''), R = Dat)
Dat = read.geno(paste(Outputname,'_uncorrected', '.geno', sep = ''))

# Read in vcf and extract the genotypes 
VCF<-read.vcfR(VCF_file)
gt<-extract.gt(VCF)

# Get the number of individuals in our dataset 
n.ind <- as.numeric(nrow(Dat))

# Pull the SNP names from the VCF and name our geno file SNPs
SNPs<- rownames(gt)
SNPs<-as.data.frame(SNPs)
colnames(Dat) = SNPs$SNPs



######################################################
##### Analyze Polymorphic Sites and Missing Data #####
######################################################

# We will use a stringent threshold that requires at least 50% of the genotypes to be present 

u9 <-  apply(Dat, 2, FUN = function(x) length( unique(x[x != 9 ]) ) ) 
Dat9 <- Dat[, u9 > 1]
na9 <-  apply(Dat9, 2, FUN = function(x) sum( x == 9 ) )

threshold <- round(n.ind/2)
loc_boo <- na9 < threshold
Dat9 <- Dat9[, na9 < threshold]
u9 <-  apply(Dat9, 2, FUN = function(x) length( unique(x[x != 9 ]) ) ) 
Dat9 <- Dat9[, u9 > 1]
loc_boo[loc_boo == TRUE] <- u9 > 1

# What are we left with (Individuals, SNPS)
dim(Dat9)


write.geno(output.file = paste(Outputname,'9', '.geno', sep = ''), R = Dat9 )
# Write out the matrix without the bad eggs

###############################
##### Impute Missing Data #####
###############################

obj.snmf <- snmf(paste(Outputname,'9', '.geno', sep = ''), K = 1:10, rep = 1, alpha = 100, CPU = 4, entropy = T,  project = "new")
plot(obj.snmf, cex = 1.2, col = "navy", pch = 19)

k = 3 # Change me based on cross-entropy plot inspection

obj <- snmf(paste(Outputname,'9', '.geno', sep = ''), K = k, rep = 50, alpha = 100, CPU = 4, entropy = T,  project = "new")
ce = cross.entropy(obj, K = k)
best = which.min(ce)
# Get the cross-entropy for the chosen K value and take the run with the lowest ce

### !!! Impute It !!! ###

dat.complete = Dat9

impute(obj, paste(Outputname,'9', '.geno', sep = ''), method = 'mode', K=k, run = best)
dat.complete = read.lfmm(paste(Outputname,'9', '.lfmm_imputed.lfmm', sep = ''))
colnames(dat.complete) <- colnames(Dat9)
write.lfmm(output.file = paste(Outputname,'c', '.lfmm', sep = ''), R = dat.complete)
lfmm2geno(paste(Outputname,'c', '.lfmm', sep = ''))

# USE the _c file for GEA analysis 

