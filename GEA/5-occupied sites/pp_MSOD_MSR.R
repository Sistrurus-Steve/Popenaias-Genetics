setwd("/media/stevehein/SteveHeinDATA/Popenaias/final_data/anaylses/GEA")

library(spdep)
library(adespatial)
library(LEA)
library(raster)
library(viridis)

Samp_log <- read.csv('pp_WorldClimData.csv')
Coord <- Samp_log[,3:4]
Coord <- as.matrix(Coord)
SNPs <- read.delim('pp_SNPs.txt', header = FALSE)
Env <- Samp_log[,6:23]
Loci <- read.lfmm(input.file = 'ppc.lfmm')


# Calculate MEM vectors and values
# --------------------------------
nb      <- graph2nb(gabrielneigh(Coord, 10), sym=TRUE)  # Gabriel graph: neighbor definition
listW   <- nb2listw(nb, style="W")                  # Spatial weights matrix
disttri <- nbdists(nb, Coord)                       # Add longlat=T for lat/long coordinates
fdist   <- lapply(disttri, function(x) x^(-1))      # Use inverse distance weights
listW   <- nb2listw(nb, glist=fdist, style="W")     # Revised spatial weights matrix
tmp     <- scores.listw(listW, MEM.autocor = "all") # Eigenanalysis

mem     <- list(vectors = as.matrix(tmp), values = attr(tmp, "values")) # MEM eigenvectors and eigenvalues
mem$values <- mem$values / abs(sum(mem$values))     # Rescale eigenvalues to Moran's I

R.YV <- cor(Loci, mem$vectors, use="pairwise.complete.obs") # R.YV = Correlation with MEM axes     
S    <- apply(R.YV^2, 2, mean)  

cutoffs <- abs(qnorm(c(0.05, 0.01, 0.001)/2))  # Cutoffs (can be modified!)
cutoffs

# Calculate z-scores for power spectra
# ------------------------------------
Dev <- sweep(R.YV^2, 2, S, "/") - 1            # Subtract average power spectrum from each locus.
Dev[Dev > 0] <- 0                              # Set positive deviations to zero.
Dev <- apply(Dev, 1, sum)                      # Sum of negative deviations
z <- scale(Dev)

# Calculate p-values from z-scores
pval <- pnorm(abs(z), lower.tail = FALSE)

cutoff.msod <- cutoffs[1]                              # Just the middle cutoff of 0.01

Candidates.msod <- c(1:ncol(Loci))[abs(z)>cutoff.msod]

# Get the candidates and their z-scores 
Cands <- data.frame(SNPs[Candidates.msod,])
Cands$Z.score <- z[Candidates.msod,]
Cands$Pval <- pval[Candidates.msod,]
Cands$locus <- gsub('_.*', '', Cands$SNPs.Candidates.msod...)
Cand_loc <- Cands[!duplicated(Cands$locus),]

# Set a cutoff & number of permutations for MSR
# ---------------------------------------------
cutoff.msr <- 0.05    # Set a less stringent cutoff
nPerm <- 999          # Set number of permutations for MSR test (may choose e.g. 499 or 999)

# MEM correlations for Env and coordinates (as spu  rious predictors)
# -----------------------------------------------------------------
get.pvalue.msr <- function(r.XV=R.XV, r.YV=R.YV, nPerm=999)
{
  R.XV.rand <- matrix(r.XV, nPerm, ncol(r.XV), byrow=TRUE) 
  R.XV.rand <- R.XV.rand * sample(c(-1,1), length(R.XV.rand), replace=TRUE)
  Cor.obs <- abs(as.vector(r.YV %*% t(r.XV)))
  Cor.rand <- abs(r.YV %*% t(R.XV.rand))
  P.values.MSR <- apply((cbind(Cor.obs,Cor.rand) >= Cor.obs), 1, mean)
  P.values.MSR
}


MSR.outliers <- list()
N.Env <- length(Env)
for (i in seq(1:length(Env))) {
R.XV.Env <- cor(Env[i], mem$vectors)
R.XV.xcoord <- cor(Coord[,1], mem$vectors)
R.XV.ycoord <- cor(Coord[,2], mem$vectors)


b.Env <- get.pvalue.msr(r.XV=R.XV.Env, r.YV=R.YV[Candidates.msod,], nPerm=nPerm)
b.X <- get.pvalue.msr(r.XV=R.XV.xcoord, r.YV=R.YV[Candidates.msod,], nPerm=nPerm)
b.Y <- get.pvalue.msr(r.XV=R.XV.ycoord, r.YV=R.YV[Candidates.msod,], nPerm=nPerm)

# Get the outliers 
MSR.outliers[[i]] <- data.frame(names(b.Env)[b.Env < cutoff.msr], b.Env[b.Env < cutoff.msr],
                                colnames(Env[i]), TRUE)
colnames(MSR.outliers[[i]]) <- c('SNP', 'P-value','Predictor', 'Candidate')
}

# Combine candidates into a single dataframe
MSR <- data.table::rbindlist(MSR.outliers, fill = TRUE)
MSR <- MSR[,-4]

# Reorder before removing any duplicates
MSR <- MSR[order(MSR$`P-value`, decreasing = FALSE)]

# Remove duplicates
MSR <- MSR[!duplicated(MSR$SNP)]

# Change SNP names for comparison with SNP map 
MSR_cands <- as.numeric(gsub('V', '', MSR$SNP))

# Replace MSR candidates with real locus names
MSR$SNP <- SNPs[MSR_cands,]

MSR$locus <- gsub('_.*', '', MSR$SNP)

MSR <-  MSR[!duplicated(MSR$locus)]
MSR <- MSR[,-4]

table(MSR$Predictor)

write.csv(MSR, 'pp_MSPcandidates.csv', row.names = FALSE)

length(MSR)
MSR




