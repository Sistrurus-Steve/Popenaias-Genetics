setwd("/media/stevehein/SteveHeinDATA/Popenaias/final_data/anaylses/GEA")

library(psych)
library(vegan)
library(trio)
library(LEA)
library(dplyr)
library(vcfR)
library(caret)
library(robust)
library(SoDA)
library(adespatial)

###########################
##### !!! Read ME !!! #####
###########################

# This script was adapted from Forester et al., 2018 
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html
# For rda analysis the genotype matrix has to been complete (No NA's) and values can be 0, 1, or 2.
# This is the same format as in LEA
# So what we'll do is take the imputed lfmm file and use that 
# See the LEA tutorial by Olivier Francois or Frichot et al., 2015 for more information  regarding formatting/imputation

# For rda, we will use an enviromental csv with populations and coordinates for downstream visualization purposes

Gtypes <- 'ppc.lfmm' # The 012 file with no missing data (it's been imputed)
Env <- 'pp_WorldClimData.csv' # The popmap and environmental data 
Snps <- 'pp_SNPs.txt' # A file with your SNPs
Outputprefix <- 'ppRDA' # Output prefix to add to file names 
###############################################
##### Bring in files and check formatting #####
###############################################

Genotypes<-read.lfmm(Gtypes)
EnvVars<-read.csv(Env)

sum(is.na(Genotypes))
9 %in% Genotypes
# or
any(Genotypes == 9)
# checks for any NA's or 9's in your dataframe

str(EnvVars)
EnvVars$Samples<-as.character(EnvVars$Samples)
# check the structure of your environmental file and change Samples to characters

Samples<-EnvVars$Samples
row.names(Genotypes)<- Samples
row.names(Genotypes)
identical(rownames(Genotypes), EnvVars$Samples)
# makes sure that your samples match each other in your genotypes and environmental variables

##########################
##### Make a SNP map #####
##########################
# We need to read in a SNP map that we created during LFMM analysis

SNPs = read.delim(Snps, header = FALSE)

# There ya go, have fun with analysis 

################################################################
##### Check for Correlation between Environmental Variables#####
################################################################

# Find variable set that retains the most variables
# First we identify which variables we should remove 
Cor <- findCorrelation(cor(EnvVars[,6:23]), cutoff = 0.7)
SelectedVars <- EnvVars[,6:23]
SelectedVars <- SelectedVars[-Cor]
# with RDA we have to check for correlation between environmental predicts as to not bias our results

Correlation<-cor.plot(SelectedVars, xlas =  2, cex.axis = 0.75, MAR = 8)
corPlot(Correlation)

# look at correlation between environmental variables, when using 10 or less you can use the paris plot
# otherwise use the correlation plot as it is easier to visualize with more variables


colnames(SelectedVars) <- c("TWaQ","PWQ","ISO","TWeQ")

write.csv(SelectedVars, file = paste(Outputprefix, "_SelectedVars", ".csv", sep = ""), row.names = FALSE)
# write it out as a file

###########################
##### Model Selection #####
###########################
mod0 <- rda(Genotypes ~ 1, data = SelectedVars, scale = TRUE)
mod1 <- rda(Genotypes ~ ., data = SelectedVars, scale = TRUE)

mod <- ordistep(mod0, scope = formula(mod1))

# Which variables should you include?
mod 

########################
##### Run your RDA #####
########################

HL_rda<- rda(Genotypes ~ ., data = SelectedVars, scale = TRUE)
HL_rda
# runs the rda and looks at it

Rsquared<-RsquareAdj(HL_rda)
Rsquared
# look at the rsquared and adjusted rsquared values

summary(eigenvals(HL_rda, model = "constrained"))
# provides a summary of the eigenvalues for each axis

screeplot(HL_rda, main = "Eigenvalues from RDA", col = "navy")
# plots a scree plot of the eigenvalues

Fullsig <- anova.cca(HL_rda, parallel=getOption("mc.cores"))
Fullsig
# test the significance of the rda
# displays results

Axissig <- anova.cca(HL_rda, by="axis", parallel=getOption("mc.cores"))
Axissig
# find which constrained axes are significant
# this can take awhile 

vif.cca(HL_rda)
# checks the variance inflation factors, should be under 10
# the lower the better


#############################
# lets make some cool plots

EnvVars$wc2.1_30s_bio_10 <- as.factor(EnvVars$wc2.1_30s_bio_10)
Pop<- as.factor(EnvVars$Pops)
bg <- c("#FF0000", "#222568")
# pull out populations and get some colors

###############
# axes 1 & 2
plot(HL_rda, type="n", scaling=3)
points(HL_rda, display="species", pch=20, cex=0.7, col="blue1", scaling=3) #displays SNPs
points(HL_rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[EnvVars$wc2.1_30s_bio_10]) #displays individuals
text(HL_rda, scaling=3, display="bp", col="356689", cex=1.25)  #displays predictors
legend("bottomright", legend=levels(EnvVars$wc2.1_30s_bio_10), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend

###############
# axes 1 & 3

plot(HL_rda, type="n", scaling=3, choices = c(1,3))
points(HL_rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices = c(1,3)) #displays SNPs
points(HL_rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[EnvVars$Habitat], choices = c(1,3)) #displays individuals
text(HL_rda, scaling=3, display="bp", col="#0868ac", cex=1, choices = c(1,3))  #displays predictors
legend("bottomright", legend=levels(EnvVars$Habitat), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend

################################
##### Testing for Outliers #####
################################

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# writes a function identify outliers
# Your choices should be your significant axes 
load.rda <- scores(HL_rda, choices=c(1:4), display="species") 
# extract SNP loadings from the significant axes

hist(load.rda[,1], main="Loadings on Axis 1", xlab = NULL)
hist(load.rda[,2], main="Loadings on Axis 2", xlab = NULL)
hist(load.rda[,3], main="Loadings on Axis 3", xlab = NULL)
hist(load.rda[,4], main="Loadings on Axis 4", xlab = NULL)
# plot histograms of the loadings, make sure that they are relatively normal
# If not, correct them 

cand1 <- outliers(load.rda[,1],2.5)
cand2 <- outliers(load.rda[,2],2.5)
cand3 <- outliers(load.rda[,3],2.5)
cand4 <- outliers(load.rda[,4],2.5)
# apply the function you wrote to each significant axis
# the number by itself indicates how many standard deviations to go
# choosing higher numbers i.e. 3 looks for loci under strong selection, lesser i.e. 2 will be more lenient

ncand <- length(cand1) + length(cand2) +length(cand3) + length(cand4)
ncand
# lets see how many outliers we have

###################################
##### Organize into Dataframe #####
###################################

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(1,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(1,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(1,times=length(cand4)), names(cand4), unname(cand4))
colnames(cand1) <- colnames(cand2) <- colnames(cand3)<- colnames(cand4)<- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3, cand4)
cand$snp <- as.character(cand$snp)

cand$snp

# Now we combine this data with our environmental variables

foo <- matrix(nrow=(ncand), ncol=4)
colnames(foo) <- c("TWaQ","PWQ","ISO","TWeQ")
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- Genotypes[,nam]
  foo[i,] <- apply(SelectedVars,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)
cand

#################################
##### Investigate Canidates #####
#################################

length(cand$snp[duplicated(cand$snp)])
# look for any canidate SNPs

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2])
# check for duplicates on axis 1

table(foo[foo[,1]==2,2])
# check for duplicates on axis 2

cand <- cand[!duplicated(cand$snp),]
# remove any duplicate SNPs

for (i in 1:length(cand$snp)) {
      bar <- cand[i,]
      cand[i,8] <- names(which.max(abs(bar[4:7]))) # makes a column with predictor
      cand[i,9] <- max(abs(bar[4:7])) } #makes a column with correlation
# Will tell us which predictor each SNP is associated with

colnames(cand)[8] <- "predictor"
colnames(cand)[9] <- "correlation"
# Adds names to the columns we just made

table(cand$predictor) 
# Lets see which variables our data are most associated with

write.csv(cand, "RefCanidateSNPs_Full.csv")
##############################
##### Compare to SNP Map #####
##############################

## I just did by hand because fixing code would take longer... NO ISSUE WITH DATA ITSELF

#candSNP<- cand$snp
#candSNP<-gsub(pattern = "V", "", candSNP)
#candSNP<-as.data.frame.integer(candSNP)
#colnames(candSNP)<- "Position"
#candSNP$Position<-as.double(candSNP$Position)
# takes the V off, you are left with numbers only
# turns into a dataframe for data manipulation

#Mapped<- inner_join(candSNP, SNPs, by = "Position")
# maps our SNPs that were identified as candidates

#Mapped

#write.csv(Mapped, file = "MappedRefSNPs_RDA_Full.csv")


###############################
##### Let's Plot the SNPs #####
###############################

sel <- cand$snp
env <- cand$predictor

# take out the SNPs and predictors

env[env=="TWaQ"] <- '#e31a1c'
env[env=="PWQ"] <- '#ffa500'
env[env=="ISO"] <- '#a6cee3'
env[env=="TWeQ"] <-  '#2b94db'

# add colors for plotting 

col.pred <- rownames(HL_rda$CCA$v)
# get SNP Names

for (i in 1:length(sel)) {          
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

# color code the canidate SNPs

col.pred[grep("V",col.pred)] <- '#f1eef6'
# add a color to the non-canidate SNPs

empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#e31a1c', '#ffa500','#a6cee3','#2b94db')
# more visualizaton prep


#####################
##### Plot SNPs #####
#####################

plot(HL_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(HL_rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(HL_rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(HL_rda, scaling=3, display="bp", col="356689", cex=1.25)  #displays predictors

legend("topleft", legend= c("TWaQ","PWQ","ISO","TWeQ"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend
# adds legend

##################
# axis 3 & 4
plot(HL_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices = c(3,4))
points(HL_rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices = c(1,3))
points(HL_rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices = c(1,3))
text(HL_rda, scaling=3, display="bp", col="356689", cex=1.25, choices = c(1,3))  #displays predictors

legend("topleft", legend= c("TWaQ","PWQ","ISO","TWeQ"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend
# adds legend

HL_rda

##################
# axis 1 & 4
plot(HL_rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices = c(1,4))
points(HL_rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices = c(1,4))
points(HL_rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices = c(1,4))
text(HL_rda, scaling=3, display="bp", col="356689", cex=1.25, choices = c(1,4))  #displays predictors

legend("topleft", legend= c("AMT", "TS", "PWQ", "PCQ"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend
# adds legend


# !!!!!!!!!!!!!!!!!!!!!!! Your Done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# for more information go to https://popgen.nescent.org/2018-03-27_RDA_GEA.html
