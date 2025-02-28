library('pcadapt')

###################################################
setwd("C:/Users/steve/OneDrive/Desktop/manuscripts/Popenaias/final_data-20250222T151229Z-001/final_data/anaylses/pcadapt/")
getwd()

## read in .bed file
pp_pcadapt <- read.pcadapt("./pp_FINAL.bed", type = "bed")

##Add pop names
poplist.names <- c(rep("Laredo", 22),rep("SanDiego", 23),rep("LowCan", 17),rep("Devils", 17),rep("Black", 26))
print(poplist.names)

pp_pca <- pcadapt(input = pp_pcadapt, K = 20)

plot(pp_pca, option = "screeplot")
plot(pp_pca, option = "screeplot", K = 10)
write.csv(pp_pca$scores, file = "pp_pca_pcaadapt_PCA_scores.csv")

plot(pp_pca, option = "screeplot", addlabels = TRUE, ylim = c(0, 50))

plot(pp_pca, option = "scores", pop = poplist.names)
plot(pp_pca, option = "scores", i = 3, j = 4, pop = poplist.names)

# Computing the test statistic based on PCA
pp_pca_K4 <- pcadapt(input = pp_pcadapt, K = 4)
summary(pp_pca_K4)

snp_pc<-get.pc(pp_pca_K4, outliers_qval)
snp_pc

plot(pp_pca_K4, option = "scores", pop = poplist.names)
plot(pp_pca_K4, option = "scores", i = 3, j = 4, pop = poplist.names)

############ K = 4 from above ####
#x <- pcadapt(filename, K = 4)
#summary(x)
#plot(x , option = "manhattan")
#plot(x, option = "qqplot")
#hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
#plot(x, option = "stat.distribution")

########## qvalue #################
library(qvalue)
qval <- qvalue(pp_pca_K4$pvalues)$qvalues
alpha <- 0.01
outliers_qval <- which(qval < alpha)
length(outliers_qval) #228
outliers_qval

#Benjamini-Hochberg Procedure
padj <- p.adjust(pp_pca_K4$pvalues,method="BH")
alpha <- 0.01
outliers_BH <- which(padj < alpha)
length(outliers_BH) #48

#Bonferroni correction
padj <- p.adjust(pp_pca_K4$pvalues,method="bonferroni")
alpha <- 0.01
outliers_Bon <- which(padj < alpha)
length(outliers_Bon) #2


#Association between PCs and outliers
snp_pc_Bon <- get.pc(pp_pca_K4, outliers_Bon)
snp_pc_Bon

snp_pc_BH <- get.pc(pp_pca_K4, outliers_BH)
snp_pc_BH

# Manhattan Plot - A Manhattan plot displays -log10 of the p-values.
pdf("pp_pcadapt_K4_manhattan.pdf")
plot(pp_pca_K4, option = "manhattan", threshold = -log10(0.01)) 




