## ----message=FALSE, warning=FALSE, echo = FALSE-------------------------------
library("knitr")
require("Ravages")

## -----------------------------------------------------------------------------
#GRR calculated using the same formula as in SKAT,
#with values in the second group of cases being twice the values from the first one

GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov, n.case.groups = 2,
                      GRR.multiplicative.factor=2, select.gene = "R1")
GRR.del[,1:5]

#Calculation of genotype frequencies in the two groups of cases and the controls group
#The previous GRR matrix is used with a multiplicative model of the disease
#The prevalence in each group of cases is 0.001

geno.freq.groups <- genotypic.freq(genes.maf = Kryukov, select.gene="R1", 
                                   GRR.het = GRR.del, prev = c(0.001, 0.001),
                                   genetic.model = "multiplicative")
str(geno.freq.groups)

#frequencies of the homozygous alternative genotype in the different groups
geno.freq.groups$freq.homo.alt[,1:5]

## -----------------------------------------------------------------------------
#MAF calculation for the first five variants
geno.freq.groups$freq.homo.alt[,1:5] + 0.5*geno.freq.groups$freq.het[,1:5]

## -----------------------------------------------------------------------------
x <- random.bed.matrix(genes.maf = Kryukov, size = c(1000, 500, 500), replicates = 5,
                       prev = c(0.001, 0.001), GRR.matrix.del = GRR.del, 
                       p.causal = 0.5, p.protect = 0, same.variant = FALSE,
                       genetic.model = "multiplicative", select.gene = "R1")
x
table(x@ped$pheno)
table(x@snps$genomic.region)

## -----------------------------------------------------------------------------
#Load LCT dataset for haplotype matrix
data(LCT.haplotypes)

#Simulation based on haplotypes frequencies
#for the variants in the LCT gene in the EUR population
LCT.gene.hap <- LCT.hap[which(LCT.sample$super.population=="EUR"),
                        which(LCT.snps$pos>=136545410 & LCT.snps$pos<=136594750)]
#Individuals from EUR
LCT.sample.EUR <- subset(LCT.sample, super.population=="EUR")
#Matrix of haplotype frequencies
LCT.freqs <- sapply(unique(LCT.sample.EUR$population), function(z) 
                    ifelse(LCT.sample.EUR$population==z, 
                           1/table(LCT.sample.EUR$population)[z], 0))
#Simulation of genetic data for five groups of 50 individuals
x <- rbm.haplos.freqs(haplos=LCT.gene.hap, freqs=LCT.freqs, size=rep(50,5), replicates=5)

#Simulation of 100 controls, and two groups of 50 cases with 30 causal variants
#and with the second group having half h2 and twice the prevalence
#compared to the first one
#5 replicates are performed and causal variants are sampled once
x <- rbm.haplos.thresholds(haplos=LCT.gene.hap, nb.causal = 30, h2=c(0.01,0.01,0.02),
                           prev=c(1,0.01,0.005), size=c(100, 50, 50), replicates=5,
                           rep.by.causal = 5)

## ---- eval = F----------------------------------------------------------------
#  #Simulations using GRR values
#  #The two groups of cases are genetically heterogeneous
#  GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov, n.case.groups = 2,
#                        GRR.multiplicative.factor=2, select.gene = "R1")
#  
#  x.GRR <- random.bed.matrix(genes.maf = Kryukov, size = c(200, 100, 100),
#                            prev = c(0.001, 0.001), GRR.matrix.del = GRR.del,
#                            p.causal = 0.3, p.protect = 0, same.variant = TRUE,
#                            replicates = 100, genetic.model = "multiplicative",
#                            select.gene = "R1")
#  table(x.GRR@ped$pheno)
#  #Selection of rare variants: MAF<1% in the entire sample
#  x.GRR.filter <- filter.rare.variants(x.GRR, maf.threshold=0.01, filter="whole")
#  
#  #Comparing each group of cases to the controls with WSS
#  x.GRR.1 <- select.inds(x.GRR.filter, x.GRR.filter@ped$pheno %in% c(0,1))
#  x.GRR.2 <- select.inds(x.GRR.filter, x.GRR.filter@ped$pheno %in% c(0,2))
#  #Construction of null models
#  H0.cas1 <- NullObject.parameters(x.GRR.1@ped$pheno, ref.level = 0,
#                                   RVAT = "burden", pheno.type = "categorial")
#  H0.cas2 <- NullObject.parameters(x.GRR.2@ped$pheno, ref.level = 0,
#                                   RVAT = "burden", pheno.type = "categorial")
#  #Run WSS
#  wss.cas1 <- burden(x.GRR.1, H0.cas1, burden="WSS", cores = 1)
#  mean(wss.cas1$p.value<0.05)
#  wss.cas2 <- burden(x.GRR.2, H0.cas2, burden="WSS", cores = 1)
#  mean(wss.cas2$p.value<0.05)
#  
#  
#  #Power of SKAT with each group considered separately
#  #Null model
#  skat.null.3gps <- NullObject.parameters(x.GRR.filter@ped$pheno,
#                                          RVAT = "SKAT", pheno.type = "categorial")
#  skat.3gps <- SKAT(x.GRR.filter, skat.null.3gps,
#                    params.sampling = list(perm.target = 100, perm.max = 5e4))
#  mean(skat.3gps$p.value<0.05)

