## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  tidy.opts = list(width.cutoff = 100), tidy = TRUE
)

## ----setup--------------------------------------------------------------------
library(hsrecombi) 
library(AlphaSimR)
library(doParallel)
library(ggplot2)
library(rlist)
Sys.time()

## ----global-------------------------------------------------------------------
# Number of SNPs (e.g., p = 1500 resembles an average bovine chromosome)
p <- 150
# Number of progeny in each half-sib family 
n <- 1000
# Directory for (simulated) data
path.dat <- 'data'
dir.create(path.dat, showWarnings = FALSE)
# Directory for output
path.res <- 'results'
dir.create(path.res, showWarnings = FALSE)
# Number of computing clusters allocated
nclust <- 2

## ----alphasim-----------------------------------------------------------------
founderPop <- runMacs2(nInd = 1000, nChr = 2, segSites = p)
SP <- SimParam$new(founderPop)
SP$setSexes("yes_sys")
# Enable tracing location of recombination events
SP$setTrackRec(TRUE)  

pop <- newPop(founderPop)
N <- 10; ntotal <- N * n
my_pop <- selectCross(pop = pop, nFemale = 500, nMale = N, use = "rand", nCrosses = ntotal)

save(list = c('SP', 'founderPop', 'pop', 'my_pop', 'ntotal'), file = file.path(path.dat, 'pop.RData'))

## ----genetic-data-------------------------------------------------------------
PAT <- my_pop@father
rown <- paste(rep(unique(PAT), each = 2), c(1, 2), sep = '_')
H.pat <- pullSegSiteHaplo(pop)[rown, ]
X <- pullSegSiteGeno(my_pop)
# Physical position of markers in Mbp
map.snp <- lapply(founderPop@genMap, function(z) z * 100)

## ----plink-format-------------------------------------------------------------
map <- data.frame(Chr = rep(1:length(map.snp), unlist(lapply(map.snp, length))), 
                  Name = paste0('SNP', 1:length(unlist(map.snp))),
                  locus_Mb = unlist(map.snp), 
                  locus_bp = unlist(map.snp) * 1e+6)

colnames(X) <- map$Name
FID <- 'FAM001'
IID <- my_pop@id
MAT <- my_pop@mother
SEX <- 2
PHENOTYPE <- -9

for(chr in 1:2){
  write.table(map[map$Chr == chr, ], file.path(path.dat, paste0('map', chr, '.map')), 
              col.names = F, row.names = F, quote = F)
  write.table(cbind(FID, IID, PAT, MAT, SEX, PHENOTYPE, X[, map$Chr == chr]), 
              file.path(path.dat, paste0('hsphase_input', chr, '.raw')), col.names = T, row.names = F, quote = F) 
}

## ----parallel-computing-------------------------------------------------------
cl <- makeCluster(nclust)
registerDoParallel(cl)

## ----recombination-rate-------------------------------------------------------
out <- foreach(chr = 1:2, .packages = 'hsrecombi') %dopar% {
  
  # 1: Physical  map
  map <- read.table(file.path(path.dat, paste0('map', chr, '.map')), col.names = c('Chr', 'Name', 'locus_Mb', 'locus_bp'))
  map$SNP <- 1:nrow(map)
  locus_Mb <- map$locus_Mb

  # 2: Genotype matrix
  genomatrix <- data.table::fread(file.path(path.dat, paste0('hsphase_input', chr, '.raw')))
  X <- as.matrix(genomatrix[, -c(1:6)])
  
  # 3: Assign daughters to sire IDs
  daughterSire <- genomatrix$PAT
  
  # 4: Estimate sire haplotypes and format data 
  hap <- makehappm(unique(daughterSire), daughterSire, X)
  save('hap', file = file.path(path.res, paste0('hsphase_output_chr', chr, '.Rdata')))
  
  # Check order and dimension
  io <- sapply(1:nrow(map), function(z){grepl(x = colnames(X)[z], pattern = map$Name[z])})
  if(sum(io) != nrow(map)) stop("ERROR in dimension")
  
  # 5: Estimate recombination rates
  res <- hsrecombi(hap, X, map$SNP)
  final <- editraw(res, map)
  save(list = c('final', 'locus_Mb'), file = file.path(path.res, paste0("Results_chr", chr, ".RData")))
  
  ifelse(nrow(final) > 0, 'OK', 'no result')
}
print(which(unlist(out) == 'OK'))

## ----misplaced----------------------------------------------------------------
# 6a: Filter SNPs with unusually large recombination rate to neighbouring (30) SNPs
excl <- foreach(chr = 1:2, .packages = 'hsrecombi') %dopar% {
  load(file.path(path.res, paste0("Results_chr", chr, ".RData")))
  checkCandidates(final)
}

# 6b: Heatmap plot of recombination rates for visual verification, e.g.:
chr <- 2
load(file.path(path.res, paste0("Results_chr", chr, ".RData")))
cand <- excl[[chr]][1]
win <- cand + (-100:100)
win <- win[(win >= 1) & (win <= max(final$SNP2))]

target <- final[(final$SNP1 %in% win) & (final$SNP2 %in% win), ]

ggplot(data = target, aes(SNP2, SNP1, fill = theta)) + 
  geom_tile() +
  xlab("Locus 2") +
  ylab("Locus 1") +
  coord_equal() + 
  scale_y_continuous(trans = "reverse") +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_line(colour = "grey", size = 0.1),
        panel.grid.minor = element_line(colour = "grey")) +
  theme(text = element_text(size = 18)) +
  scale_fill_gradientn(colours = c('yellow', 'red'), limits = c(0, 1+1e-10), na.value = 'white')
# -> nothing conspicious

excl[[1]] <- excl[[2]] <- NA

## ----genetic-position---------------------------------------------------------
pos <- foreach(chr = 1:2, .packages = 'hsrecombi') %dopar% {
  load(file.path(path.res, paste0("Results_chr", chr, ".RData")))
  out <- geneticPosition(final, exclude = excl[[chr]])
  list(pos.cM = c(out, rep(NA, length(locus_Mb) - length(out))), pos.Mb = locus_Mb)
}

## ----stop-parallel------------------------------------------------------------
stopCluster(cl)

## ----plot-1-------------------------------------------------------------------
data <- data.frame(Chr = rep(1:length(pos), times = unlist(list.select(pos, length(pos.Mb)))), 
                   Mbp = unlist(list.select(pos, pos.Mb)), cM = unlist(list.select(pos, pos.cM)))
ggplot(data, aes(x = Mbp, y = cM)) + geom_point(na.rm = T) + facet_grid(Chr ~ .)

## ----selection----------------------------------------------------------------
chr <- 2

## ----check-simulated----------------------------------------------------------
load(file.path(path.dat, 'pop.RData'))
co.pat <- matrix(0, ncol = p, nrow = ntotal)
for(i in 1:ntotal){
  if(nrow(SP$recHist[[1000 + i]][[chr]][[2]]) > 1){
    # 1. line contains 1 1 by default
    loci <- SP$recHist[[1000 + i]][[chr]][[2]][-1, 2] 
    co.pat[i, loci] <- 1 
  }
}
probRec <- colMeans(co.pat)
sim.cM <- cumsum(probRec) * 100

## ----check-deterministic------------------------------------------------------
load(file.path(path.res, paste0('hsphase_output_chr', chr, '.Rdata')))
hsphase.cM <- c(0, cumsum(hap$probRec)) * 100

## ----plot-2-------------------------------------------------------------------
par(mar=c(5.1, 4.1, 4.1, 9.1), xpd = TRUE)
plot(pos[[chr]]$pos.Mb, pos[[chr]]$pos.cM, xlab = 'physical position (Mbp)', ylab = 'genetic position (cM)', 
     ylim = range(c(sim.cM, hsphase.cM, pos[[chr]]$pos.cM), na.rm = T), pch = 20)
points(pos[[chr]]$pos.Mb, sim.cM, pch = 20, col = 4)
points(pos[[chr]]$pos.Mb, hsphase.cM, pch = 20, col = 8)
legend('topleft', inset=c(1.01,0), legend = c('simulated', 'deterministic', 'likelihood-based'), pch = 20, col = c(4, 8, 1), bty = 'n')

## ----cleanup------------------------------------------------------------------
unlink(path.dat, recursive = TRUE)
unlink(path.res, recursive = TRUE)
Sys.time()

