library(adegenet)
library(pegas)
library(ape)
library(vcfR)

setwd("C:/Users/steve/OneDrive/Desktop/manuscripts/Popenaias/final_data-20250222T151229Z-001/final_data/anaylses/DAPC")
getwd()

## Neutral ONLY dataset ##
vcf <- read.vcfR("./neutral_pp_snps.vcf")

# check your data
vcf 
# convert to genind
pp_genind <- vcfR2genind(vcf)
# double check the data again
pp_genind
# find clusters, saving all PC
set.seed(7177)
grp1<-find.clusters(seed = 123, pp_genind, max.n.clust = 25)
90
4
# Kstat, this will return the BIC value for each of the 25 possible clusters you created using find.cluster. 
grp1$Kstat
# stat, the BIC for your choosen k value
grp1$stat
# grp, which cluster (1 -4) each indivudal is assigned to 
grp1$grp
# size, number of individuals assigned to each cluster 
grp1$size
# perform a discriminate analyses of principal components (DAPC) using the clusters you just defined. i.e. grp1$grp
# determining k-means clusters benefits from retaining all variation a DAPC does not
# you will again be asked to reatin a certain amount of PCs and again retain 120
# then you'll be asked to define how many discriminant functions to retain, retain all 4
dapc.NoPrio.4<-dapc(seed = 123, pp_genind, grp1$grp)
90
3
# the "optim.a.score" will take the information you just saved into the dapc.np.5 object and determine the 
# optimimum number of PCs to retain 
temp1<-optim.a.score(seed = 99, dapc.NoPrio.4, n.sim = 100, smart = T)
# the optimimum output is displayed graphifically or can be display 
# the optimimum output is displayed graphifically or can be display hitting enter below
temp1$best
# see the percent variation explained by each PC. remember you are only keeping the first 9 though for 
# this analysis
dapc.NoPrio.4$pca.eig/sum(dapc.NoPrio.4$pca.eig)*100
# rerun the DAPC using the optimial number of PCs as determined in "temp1" i.e. temp1$best
# again retain all discriminant functions, which should be 3
dapc.NoPrio.4<-dapc(pp_genind, grp1$grp, n.pca = temp1$best)
2
# plot the results of the dapc. 
set.seed(909)
scatter(dapc.NoPrio.4, bg="white", scree.da=F, scree.pca = F,
        posi.da = "bottomright", legend=TRUE, posi.leg = "topleft", solid=1, cstar = 1, 
        clabel = 0, cellipse = 0)
dapc.NoPrio.4$grp


#########################################################################

vcf <- read.vcfR("./neutral_pp_snps.vcf")
# check your data
vcf 
# convert to genind
pp_genind <- vcfR2genind(vcf)
# double check the data again
pp_genind
# find clusters, saving all PC
set.seed(789)
grp2<-find.clusters(pp_genind, max.n.clust = 25)
90
5
# Kstat, this will return the BIC value for each of the 25 possible clusters you created using find.cluster. 
grp2$Kstat
# stat, the BIC for your choosen k value
grp2$stat
# grp, which cluster (1 -4) each indivudal is assigned to 
grp2$grp
        # size, number of individuals assigned to each cluster 
grp2$size
# perform a discriminate analyses of principal components (DAPC) using the clusters you just defined. i.e. grp1$grp
# determining k-means clusters benefits from retaining all variation a DAPC does not
# you will again be asked to reatin a certain amount of PCs and again retain 120
# then you'll be asked to define how many discriminant functions to retain, retain all 4
dapc.NoPrio.5<-dapc(seed = 99, pp_genind, grp2$grp)
90
4
# the "optim.a.score" will take the information you just saved into the dapc.np.5 object and determine the 
# optimimum number of PCs to retain 
temp2<-optim.a.score(seed = 99, dapc.NoPrio.5, n.sim = 100)
# the optimimum output is displayed graphifically or can be display hitting enter below
temp2$best
# see the percent variation explained by each PC. remember you are only keeping the first 9 though for 
# this analysis
dapc.NoPrio.5$pca.eig/sum(dapc.NoPrio.5$pca.eig)*100
# rerun the DAPC using the optimial number of PCs as determined in "temp1" i.e. temp1$best
# again retain all discriminant functions, which should be 3
dapc.NoPrio.5<-dapc(pp_genind, grp2$grp, n.pca = temp2$best)
4
# plot the results of the dapc. 
set.seed(909)
scatter(dapc.NoPrio.5, bg="white", scree.da=FALSE, scree.pca = FALSE,
        posi.da = "bottomright", legend=TRUE, posi.leg = "topleft", solid=1, cstar = 1, 
        clabel = 0, cellipse = 0)
dapc.NoPrio.5$grp

vcf

####################### K3 ################
vcf <- read.vcfR("./populations.snps.vcf")
# check your data
vcf 
# convert to genind
pp_genind <- vcfR2genind(vcf)
# double check the data again
pp_genind
# find clusters, saving all PC
set.seed(456)
grp3<-find.clusters(pp_genind, max.n.clust = 25)
120
3
# Kstat, this will return the BIC value for each of the 25 possible clusters you created using find.cluster. 
grp3$Kstat
# stat, the BIC for your choosen k value
grp3$stat
# grp, which cluster (1 -4) each indivudal is assigned to 
grp3$grp
# size, number of individuals assigned to each cluster 
grp3$size
# perform a discriminate analyses of principal components (DAPC) using the clusters you just defined. i.e. grp1$grp
# determining k-means clusters benefits from retaining all variation a DAPC does not
# you will again be asked to reatin a certain amount of PCs and again retain 120
# then you'll be asked to define how many discriminant functions to retain, retain all 4
dapc.NoPrio.3<-dapc(seed = 99, pp_genind, grp3$grp)
120
2
# the "optim.a.score" will take the information you just saved into the dapc.np.5 object and determine the 
# optimimum number of PCs to retain 
temp3<-optim.a.score(seed = 99, dapc.NoPrio.3, n.sim = 20)
# the optimimum output is displayed graphifically or can be display hitting enter below
temp3$best
# see the percent variation explained by each PC. remember you are only keeping the first 9 though for 
# this analysis
dapc.NoPrio.3$pca.eig/sum(dapc.NoPrio.3$pca.eig)*100
# rerun the DAPC using the optimial number of PCs as determined in "temp1" i.e. temp1$best
# again retain all discriminant functions, which should be 3
dapc.NoPrio.3<-dapc(pp_genind, grp3$grp, n.pca = 2)
2
# plot the results of the dapc. 
set.seed(909)
scatter(dapc.NoPrio.3, bg="white", scree.da=FALSE, scree.pca = FALSE,
        posi.da = "bottomright", legend=TRUE, posi.leg = "topleft", solid=1, cstar = 0, 
        clabel = 0, cellipse = 2.5)
dapc.NoPrio.3$grp


########## K=6 ###############



vcf <- read.vcfR("./populations.snps.vcf")
# check your data
vcf 
# convert to genind
pp_genind <- vcfR2genind(vcf)
# double check the data again
pp_genind
# find clusters, saving all PC
set.seed(456)
grp4<-find.clusters(pp_genind, max.n.clust = 25)
120
6
# Kstat, this will return the BIC value for each of the 25 possible clusters you created using find.cluster. 
grp4$Kstat
# stat, the BIC for your choosen k value
grp4$stat
# grp, which cluster (1 -4) each indivudal is assigned to 
grp4$grp
# size, number of individuals assigned to each cluster 
grp4$size
# perform a discriminate analyses of principal components (DAPC) using the clusters you just defined. i.e. grp1$grp
# determining k-means clusters benefits from retaining all variation a DAPC does not
# you will again be asked to reatin a certain amount of PCs and again retain 120
# then you'll be asked to define how many discriminant functions to retain, retain all 4
dapc.NoPrio.6<-dapc(seed = 99, pp_genind, grp4$grp)
120
6
# the "optim.a.score" will take the information you just saved into the dapc.np.5 object and determine the 
# optimimum number of PCs to retain 
temp6<-optim.a.score(seed = 99, dapc.NoPrio.6, n.sim = 20)
# the optimimum output is displayed graphifically or can be display hitting enter below
temp6$best
# see the percent variation explained by each PC. remember you are only keeping the first 9 though for 
# this analysis
dapc.NoPrio.6$pca.eig/sum(dapc.NoPrio.6$pca.eig)*100
# rerun the DAPC using the optimial number of PCs as determined in "temp1" i.e. temp1$best
# again retain all discriminant functions, which should be 3
dapc.NoPrio.6<-dapc(pp_genind, grp4$grp, n.pca = temp6$best)
5
# plot the results of the dapc. 
set.seed(909)
scatter(dapc.NoPrio.6, bg="white", scree.da=FALSE, scree.pca = FALSE,
        posi.da = "bottomright", legend=TRUE, posi.leg = "topleft", solid=1, cstar = 0, 
        clabel = 0, cellipse = 2.5)
dapc.NoPrio.6$grp




