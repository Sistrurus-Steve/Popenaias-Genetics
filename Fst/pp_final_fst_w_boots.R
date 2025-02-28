library("hierfstat")
library('adegenet')

#setwd("/media/stevehein/SteveHeinDATA/Lampsilis_raf/genetics/stacks/output/M2n2m30_full_set/populations")
setwd("/media/stevehein/SteveHeinDATA/Popenaias/final_data/anaylses/Fst/")

# read genpop file in as genind
data<-read.genepop("neutral_pp_snps.gen")

#check order of populations 
data@pop

# create population vector **IN CORRECT ORDER**
PopNames <- c('Laredo',"SanDeigo",'LowerCanyon','Devils','Black')
# set pop names in genind object 
popNames(data) <- PopNames
popNames(data)
# try again 
data@pop
# summary of the data
data

# convert genind to hierstat data frame
hierdata<-genind2hierfstat(data, pop = data@pop)

# generate pairwise Fst values; this will take a few minutes with a standard RADseq SNP dataset
fst <- pairwise.WCfst(hierdata [,-2], diploid = TRUE)
fst

# bootstrap fst values for static significance; below 95CI. This might also take a few minutes to run. 
bootsfst<-hierfstat::boot.ppfst (hierdata [,-2], nboot = 1000, quant = c(0.025, 0.975), diploid = TRUE)
bootsfst


write.csv(bootsfst, file="pp_neutral_fst_boots.csv")

write.csv(fst, file="pp_nuetral_fst.csv")
  