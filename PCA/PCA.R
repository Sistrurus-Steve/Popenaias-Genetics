library(SNPRelate)

setwd("/path/to/wd")
getwd()

############################ combined outliers ##################

Com.vcf.fn <- "combined.vcf"
snpgdsVCF2GDS(Com.vcf.fn, "C.gds", method = "biallelic.only")

Com <- snpgdsOpen("C.gds")
C.pca <- snpgdsPCA(Com, num.thread=4)

pc.percent <- C.pca$varprop*100
head(round(pc.percent, 2))

plot(C.pca)

tab <- data.frame(sample.id = C.pca$sample.id,
                  EV1 = C.pca$eigenvect[,1], # the first eigenvector
                  EV2 = C.pca$eigenvect[,2], # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)


# draw
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

#add pop info
pop_code <- scan("pop_order.txt", what=character())

# get sample id
samp.id <- read.gdsn(index.gdsn(Com, "sample.id"))

# assume the order of sample IDs is as the same as population codes
cbind(samp.id, pop_code)

tab <- data.frame(sample.id = C.pca$sample.id,
                  pop = factor(pop_code)[match(C.pca$sample.id, samp.id)],
                  EV1 = C.pca$eigenvect[,1], # the first eigenvector
                  EV2 = C.pca$eigenvect[,2], # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

plot(tab$EV1, tab$EV2, col=as.integer(tab$pop),
     xlab="eigenvector 2", ylab="eigenvector 1")

legend("bottomleft", legend=levels(tab$pop), pch="o", col=1:4)

############################### PCAdapt outlier ################3

PC.vcf.fn <- "PCadapt.vcf"
snpgdsVCF2GDS(P.vcf.fn, "PC.gds", method = "biallelic.only")

PC <- snpgdsOpen("PC.gds")
PC.pca <- snpgdsPCA(PC, num.thread=4)

pc.percent <- PC.pca$varprop*100
head(round(pc.percent, 2))

plot(PC.pca)

tab <- data.frame(sample.id = PC.pca$sample.id,
                  EV1 = PC.pca$eigenvect[,1], # the first eigenvector
                  EV2 = PC.pca$eigenvect[,2], # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)


# draw
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

#add pop info
pop_code <- scan("pop_order.txt", what=character())

# get sample id
samp.id <- read.gdsn(index.gdsn(PC, "sample.id"))

# assume the order of sample IDs is as the same as population codes
cbind(samp.id, pop_code)

tab <- data.frame(sample.id = PC.pca$sample.id,
                  pop = factor(pop_code)[match(C.pca$sample.id, samp.id)],
                  EV1 = PC.pca$eigenvect[,1], # the first eigenvector
                  EV2 = PC.pca$eigenvect[,2], # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

plot(tab$EV1, tab$EV2, col=as.integer(tab$pop),
     xlab="eigenvector 2", ylab="eigenvector 1")

legend("bottomleft", legend=levels(tab$pop), pch="o", col=1:4)


############################### RDA outliers ###############


RDAvcf.fn <- "ppRDA.recode.vcf"
snpgdsVCF2GDS(RDAvcf.fn, "R.gds", method = "biallelic.only")

RDA <- snpgdsOpen("R.gds")
R.pca <- snpgdsPCA(RDA, num.thread=4)

pc.percent <- R.pca$varprop*100
head(round(pc.percent, 2))

plot(R.pca)

tab <- data.frame(sample.id = R.pca$sample.id,
                  EV1 = R.pca$eigenvect[,1], # the first eigenvector
                  EV2 = R.pca$eigenvect[,2], # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)


# draw
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

#add pop info
pop_code <- scan("pop_order.txt", what=character())

# get sample id
samp.id <- read.gdsn(index.gdsn(RDA, "sample.id"))

# assume the order of sample IDs is as the same as population codes
cbind(samp.id, pop_code)

tab <- data.frame(sample.id = R.pca$sample.id,
                  pop = factor(pop_code)[match(R.pca$sample.id, samp.id)],
                  EV1 = R.pca$eigenvect[,1], # the first eigenvector
                  EV2 = R.pca$eigenvect[,2], # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),
     xlab="eigenvector 2", ylab="eigenvector 1")
legend("topright", legend=levels(tab$pop), pch="o", col=1:4)

################################ MSOD MSR outliers ###################################

M.vcf.fn <- "MSOD.vcf"
snpgdsVCF2GDS(M.vcf.fn, "M.gds", method = "biallelic.only")

M <- snpgdsOpen("M.gds")
M.pca <- snpgdsPCA(M, num.thread=4)

pc.percent <- M.pca$varprop*100
head(round(pc.percent, 2))

plot(M.pca, col=as.integer(tab$pop))




tab <- data.frame(sample.id = M.pca$sample.id,
                  EV1 = M.pca$eigenvect[,1], # the first eigenvector
                  EV2 = M.pca$eigenvect[,2], # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)


# draw
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

#add pop info
pop_code <- scan("pop_order.txt", what=character())

# get sample id
samp.id <- read.gdsn(index.gdsn(M, "sample.id"))

# assume the order of sample IDs is as the same as population codes
cbind(samp.id, pop_code)

tab <- data.frame(sample.id = M.pca$sample.id,
                  pop = factor(pop_code)[match(M.pca$sample.id, samp.id)],
                  EV1 = M.pca$eigenvect[,1], # the first eigenvector
                  EV2 = M.pca$eigenvect[,2], # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

plot(tab$EV1, tab$EV2, col=as.integer(tab$pop),
     xlab="eigenvector 2", ylab="eigenvector 1")

legend("topright", legend=levels(tab$pop), pch="o", col=1:4)






