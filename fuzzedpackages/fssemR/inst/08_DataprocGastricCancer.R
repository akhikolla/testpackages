library(GEOquery)
library(illuminaHumanv3.db)
library(stringr)
library(biomaRt)
library(MatrixEQTL)
library(org.Hs.eg.db)
library(clusterProfiler)
library("BSgenome.Hsapiens.UCSC.hg19")
library(fssemR)
library(igraph)
library(limma)
library(synbreed)      #### Beagle are used for missing genotype imputation
library(gage)
library(httr)
library(XML)

### The distinction between primary and secondary ovarian tumors
gse1 =
  GEOquery:::parseGSEMatrix(
    "/media/xinchou/Storage/SMLfl/exp/GSE29999-GPL6947_series_matrix.txt.gz",
    destdir = "exp",
    AnnotGPL = FALSE,
    getGPL = F
  )
gse2 =
  GEOquery:::parseGSEMatrix(
    "/media/xinchou/Storage/SMLfl/exp/GSE29999-GPL6801_series_matrix.txt.gz",
    destdir = "exp",
    AnnotGPL = FALSE,
    getGPL = F
  )


# step-2 Retrieve prob expression levels and log-transform them
Exprdat = read.ilmn(files = "/media/xinchou/Storage/SMLfl/exp/GSE29998_non-normalized.txt", probeid = "ID_REF",
                    expr = "AVG_Signal")
title2GSM = as.character(gse1$eset@phenoData@data$geo_accession)
names(title2GSM) = as.character(gse1$eset@phenoData@data$title)
colnames(Exprdat) = title2GSM[colnames(Exprdat)]
ExprNorm = neqc(Exprdat)
## mapping Illumina ID to Gene ID
addrIllumina = toTable(illuminaHumanv3ARRAYADDRESS)[, c("ArrayAddress", "IlluminaID")]
colnames(addrIllumina) = c("ArrayAddrID", "IlluminaID_1")
illuminaToSymbol = toTable(illuminaHumanv3ENTREZREANNOTATED)
addrToSymbol = merge(addrIllumina, illuminaToSymbol, by.x="IlluminaID_1", by.y="IlluminaID")
addrToLocation = toTable(illuminaHumanv3ENSEMBLREANNOTATED)
addrToLocation = merge(addrIllumina, addrToLocation, by.x="IlluminaID_1", by.y="IlluminaID")
rownames(addrToLocation) = addrToLocation$IlluminaID_1
### 1 to 1 Illumina ID to ENSEMBL ID
expressedProbe = rowSums(ExprNorm$other$Detection < 0.05) > 2
Exprvarmat = ExprNorm$E[expressedProbe,]
exprIllumina = intersect(rownames(Exprvarmat), rownames(addrToLocation))
Exprvarmat = Exprvarmat[exprIllumina,]
rownames(Exprvarmat) = addrToLocation[rownames(Exprvarmat), 3]
## [1] 26486    99  ENSEMBL --> TRANSCRIPT ID

## step-3 Retrieve SNP of corresponding patients
SNPlib = read.csv(
  "/media/xinchou/Storage/SMLfl/exp/GenomeWideSNP_6.na29.annot.csv",
  comment.char = "#",
  sep = ",",
  stringsAsFactors = F
)
### only keep SNP with location information
SNPlib = SNPlib[SNPlib[,3] != "---" & SNPlib[,4] != "---",]
### build snp's has table
SNPmap = SNPlib[, c(2, 3, 4)]
SNPmap = unique(SNPmap)
rownames(SNPmap) = SNPmap[,1]
SNPhash = SNPlib[, 2]
names(SNPhash) = SNPlib[, 1]
### filter SNP's hash table
SNPhash = SNPhash[!duplicated(SNPhash)]
SNPvarmat = gse2$eset@assayData$exprs
SNPID = rownames(SNPvarmat)
FilterSNP = intersect(SNPID, names(SNPhash))
SNPvarmat = SNPvarmat[FilterSNP, , drop = F]
snpnames = SNPhash[rownames(SNPvarmat)]
names(snpnames) = NULL
rownames(SNPvarmat) = snpnames  ## SNPID --> RS____ ID
### transform character of SNP to numeric
## SNPvarmat[SNPvarmat == "AA"] = 01
## SNPvarmat[SNPvarmat == "AB"] = 02
## SNPvarmat[SNPvarmat == "BB"] = 03
SNPvarmat[SNPvarmat == "NoCall" | SNPvarmat == "NC"] = NA
### remove unchanged SNP and all Missing NA
### impute missing NA in SNP matrix
SNPvarmat = t(SNPvarmat)
SNPmap = SNPmap[colnames(SNPvarmat),c(2,3)]
colnames(SNPmap) = c("chr", "pos")
SNPmap[,2] = as.numeric(SNPmap[,2])
## dim(SNPvarmat) ## [1]    122 930002
PData2 = phenoData(gse2$eset)                     # SNP
SNPPheno = PData2@data[rownames(SNPvarmat), 10, drop = F]
SNPPheno[,1] = as.numeric(SNPPheno[,1]) - 1
colnames(SNPPheno) = c("Status")
SNPData = create.gpData(pheno = SNPPheno, geno = SNPvarmat, map = SNPmap, map.unit = "bp")
## SNPImputed = codeGeno(SNPData, impute=TRUE, impute.type="beagle", cores = 4)
SNPImputed = readRDS("./data2/ImputedSNP.rds")
## SNPvarmat = t(SNPImputed$geno)

## map back to "AA", "AB", "BB"
ImputeData = SNPImputed$geno
RawData = SNPData$geno[,colnames(ImputeData)]
RawData[RawData == "AA"] = "0"
RawData[RawData == "AB"] = "1"
RawData[RawData == "BB"] = "2"
mode(ImputeData) = "character"
for (i in 1:ncol(RawData)) {
  genomap = na.omit(unique(cbind(ImputeData[,i], RawData[,i])))
  t = genomap[,2]
  names(t) = genomap[,1]
  ImputeData[,i] = t[ImputeData[,i]]
}
mode(ImputeData) = "numeric"
SNPvarmat = t(ImputeData)

## step-4 pair-data of lung cancer and normal data
## remove unpaired "08259T2"
PData1 = phenoData(gse1$eset)                     # GE
PData2 = phenoData(gse2$eset)                     # SNP
## samples have both Gene exprs and SNP vars
sampidGE = as.character(PData1@data$source_name_ch1)
sampidSNP = as.character(PData2@data$source_name_ch1)
## 49 paired sample
sampleID = setdiff(intersect(sampidGE, sampidSNP), "08259T2")
GEix = sapply(sampleID, function(id) {
  which(sampidGE == id)
})
SNPix = sapply(sampleID, function(id) {
  which(sampidSNP == id)
})
Exprvarmat = Exprvarmat[, GEix, drop = F]
SNPvarmat = SNPvarmat[, SNPix, drop = F]

## step-5 search eQTL w.r.t gene expression by
### MatrixEQTL detect eQTL
Libensembl = useDataset("hsapiens_gene_ensembl", mart = useMart("ensembl"))
GeneLocation = getBM(
  attributes = c(
    "ensembl_gene_id",
    "chromosome_name",
    "start_position",
    "end_position"
  ),
  mart = Libensembl
)
### Collect entrezID with ensembl
Libid = bitr(
  rownames(Exprvarmat),
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = "org.Hs.eg.db"
)
### Remove 7.6% of input gene IDs are fail to map
ensembl2id = Libid$SYMBOL
names(ensembl2id) = Libid$ENSEMBL
GeneLocation = GeneLocation[GeneLocation$ensembl_gene_id %in% names(ensembl2id), , drop = F]
GeneLocation = GeneLocation[complete.cases(GeneLocation), , drop = F]
GeneLocation = unique(GeneLocation)

### Build object for MatrixEQTL
geneExprData1 = Exprvarmat[as.character(unique(GeneLocation$ensembl_gene_id)), , drop = F]
Exprmat = SlicedData$new()
colnames(geneExprData1) = NULL
Exprmat$initialize(geneExprData1)
#### SNP build
SNPmat = SlicedData$new()
SNPvarmat1 = SNPvarmat
colnames(SNPvarmat1) = NULL
SNPmat$initialize(SNPvarmat1)
GeneLoc = GeneLocation[as.character(GeneLocation$ensembl_gene_id) %in% rownames(geneExprData1), ]
GeneLoc = GeneLoc[, c(1, 2, 3, 4)]
colnames(GeneLoc) = c("geneid", "chr", "left", "right")
rownames(GeneLoc) = NULL
SNPLoc = data.frame(SNPlib[,c(2, 3, 4)])
SNPLoc[,3] = as.numeric(SNPLoc[,3])
SNPLoc = SNPLoc[complete.cases(SNPLoc), , drop = F]
SNPLoc = unique(SNPLoc)
colnames(SNPLoc) = c("snpid", "chr", "pos")
SNPLoc = SNPLoc[SNPLoc$snpid %in% rownames(SNPmat),]
Covariates = NULL
PData = PData1@data[colnames(Exprvarmat), ]
Status = as.character(PData$characteristics_ch1)
Status[Status == "tissue: Tumor"] = "tumor"
Status[Status != "tumor"] = "normal"
Covariates = rbind(Covariates,
                   status = ifelse(Status == "normal", 0, 1))
Covmat = SlicedData$new()
Covmat$initialize(Covariates)

### Prepare for MatrixEQTL
SNPInfo = SNPvarmat1[SNPLoc$snpid[SNPLoc$chr %in% c(as.character(seq(1,22)), "X", "Y", "MT")], ]
SNPS    = data.frame(id = rownames(SNPInfo), SNPInfo)
colnames(SNPS) = c("id", paste("Sample_", seq(1, ncol(SNPInfo)), sep=""))
write.table(SNPS, "./data2/SNP.txt", quote = F, row.names = F, col.names = TRUE, sep = "\t")
GE  = data.frame(id = rownames(geneExprData1), geneExprData1)
colnames(GE) = c("id", paste("Sample_", seq(1, ncol(Exprvarmat)), sep=""))
write.table(GE, "./data2/GE.txt", quote = F, row.names = F, col.names = TRUE, sep = "\t")
COV = data.frame(id = rownames(Covariates), Covariates)
colnames(COV) = c("id", paste("Sample_", seq(1, ncol(Covariates)), sep=""))
write.table(COV, "./data2/Covariates.txt", quote = F, row.names = F, col.names = TRUE, sep = "\t")
##### Location information of SNP and Gene
saveRDS(SNPLoc, "./data2/SNP.rds")
saveRDS(GeneLoc, "./data2/GE.rds")

##########################################################
## Extract Significant eQTL from MatrixEQTL result
cis_eQTL = read.csv("./data2/cis_eQTL_results_R.txt", sep = "\t", stringsAsFactors = F)
Normal = which(Status == "normal")
Tumor  = which(Status == "tumor")
### SNP with MAF > 0.05 are considered in HapMap, so filter eQTL with MAP (minor allele frequency) > 0.05
FilterByMAF = apply(SNPvarmat, 1, function(x) {
  MAF_N = min(c(sum(x[Normal] == 0), sum(x[Normal] == 1), sum(x[Normal] == 2))) / length(Normal)
  MAF_T = min(c(sum(x[Tumor] == 0), sum(x[Tumor] == 1), sum(x[Tumor] == 2))) / length(Tumor)
  (MAF_N > 0.05 & MAF_T > 0.05)
})
SNPFilterByMAF = names(which(FilterByMAF))
##CandidateGenes = union(HumanBaseRefGRN$G1, HumanBaseRefGRN$G2)
Significant_eQTLs = cis_eQTL[(cis_eQTL$FDR < 0.01 & cis_eQTL$SNP %in% SNPFilterByMAF), , drop = F]
##rownames(SNPLoc) = SNPLoc[,1]
##rownames(GeneLoc) = GeneLoc[,1]
##Distance = apply(Significant_eQTLs, 1, function(x) {
##  min(abs(SNPLoc[x[1], 3] - GeneLoc[x[2], 3]), abs(SNPLoc[x[1], 3] - GeneLoc[x[2], 4]))
##})
## cis-eQTL mapping 585 genes with cis-eQTL
Gene2eQTL = split(Significant_eQTLs$SNP, Significant_eQTLs$gene)
## center function
center = function(X) {
  apply(X, 1, function(x) {
    x - mean(x)
  })
}
## filtered out eQTL, if two eQTL have the same pattern, remove one of them
eQTLRank = lapply(Gene2eQTL, function(g) {
  t = unique(SNPvarmat[g, Tumor, drop = F])
  n = unique(SNPvarmat[g, Normal, drop = F])
  min(qr(crossprod(center(t)))$rank, qr(crossprod(center(n)))$rank)
})
## only keep non-repeated eQTL
Gene2eQTL2 = Gene2eQTL
for (i in 1:length(Gene2eQTL)) {
  FDR = NULL
  g = as.numeric(names(Gene2eQTL[i]))
  for (s in Gene2eQTL[[i]]) {
    FDR = c(FDR, Significant_eQTLs[Significant_eQTLs$SNP == s & Significant_eQTLs$gene == g, 6])
  }
  if (eQTLRank[[i]] != 0) {
    Gene2eQTL[[i]] = Gene2eQTL[[i]][which.min(FDR)]
  } else {
    Gene2eQTL[[i]] = NA
  }
}

for(i in 1:length(Gene2eQTL)) {
  if (eQTLRank[[i]] == length(Gene2eQTL[[i]])) {
    Gene2eQTL[[i]] = Gene2eQTL[[i]]
  } else if (eQTLRank[[i]] == 1 & length(Gene2eQTL[[i]]) > 1) {
    FDR = NULL
    g = as.numeric(names(Gene2eQTL[i]))
    for (s in Gene2eQTL[[i]]) {
      FDR = c(FDR, Significant_eQTLs[Significant_eQTLs$SNP == s & Significant_eQTLs$gene == g, 6])
    }
    Gene2eQTL[[i]] = Gene2eQTL[[i]][which.min(FDR)]
  } else if (eQTLRank[[i]] > 1 & length(Gene2eQTL[[i]]) > 1) {
    FDR = NULL
    g = as.numeric(names(Gene2eQTL[i]))
    for (s in Gene2eQTL[[i]]) {
      FDR = c(FDR, Significant_eQTLs[Significant_eQTLs$SNP == s & Significant_eQTLs$gene == g, 6])
    }
    index = sort.int(FDR, decreasing = F, index.return = T)$ix
    n = 0
    j = 1
    s = c(Gene2eQTL[[i]][index[j]])
    j = j + 1
    n = 1
    while (n < eQTLRank[[i]] & j < length(index)) {
      tmr = SNPvarmat[c(s, Gene2eQTL[[i]][index[j]]), Tumor, drop = F]
      nml = SNPvarmat[c(s, Gene2eQTL[[i]][index[j]]), Normal, drop = F]
      if (qr(crossprod(center(tmr)))$rank == n + 1 && qr(crossprod(center(nml)))$rank == n + 1) {
        s = c(s, Gene2eQTL[[i]][index[j]])
        n = n + 1
      }
      j = j + 1
    }
    Gene2eQTL[[i]] = s
  } else {
    Gene2eQTL[[i]] = NA
  }
}
Gene2eQTL = Gene2eQTL[!is.na(Gene2eQTL)]

## filter only protein-coding gene
annotation2Entrez = read.table(
  "./data/geneAnnotation.txt",
  sep = "\t",
  quote = "\"",
  na.strings = "-",
  fill = TRUE,
  col.names = c("GeneID", "Symbol", "TypeOfGene"), stringsAsFactors = FALSE
)
ensembl2entrez = bitr(
  names(Gene2eQTL),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)
gene4mRNA = annotation2Entrez$GeneID[annotation2Entrez$TypeOfGene == "protein-coding"]
gene4mRNA = intersect(ensembl2entrez$ENTREZID, gene4mRNA)
ensembl2entrez = ensembl2entrez[ensembl2entrez$ENTREZID %in% gene4mRNA, ]
Gene2eQTL = Gene2eQTL[unique(ensembl2entrez$ENSEMBL)]
Gene2eQTL = Gene2eQTL[sapply(Gene2eQTL, length) != 0]

## Gene feature selection, only keep those genes who are significant differentially expressed under tumor vs normal
## design = model.matrix(~ age + status, as.data.frame(t(Covariates)))
## fit0 = lmFit(GE[,-1], design)
## fite = eBayes(fit0)
## limres = topTable(fite, coef = 3, number = Inf)
## sigGene = limres[limres$adj.P.Val < 0.01, ]

## We only consider the differential GRN of genes who was detected differential expressed under tumor status
## Gene2eQTL0 = Gene2eQTL[intersect(rownames(sigGene), names(Gene2eQTL))]

### build input of FSSEM
seed = as.numeric(Sys.time())
N = sum(Status == "normal")
Ng = length(Gene2eQTL)                                                          ## 1459 genes
Nk = sum(sapply(Gene2eQTL, length))                                             ## 2146 eQTLs
set.seed(seed)
Sk = list()
## build Sk for FSSEM
index = 0
CandidateEQTLs = NULL
for (i in 1:length(Gene2eQTL)) {
  Sk[[i]] = index + seq(1:length(Gene2eQTL[[i]]))
  index = max(Sk[[i]])
  CandidateEQTLs = c(CandidateEQTLs, Gene2eQTL[[i]])
}
CandidateGenes = names(Gene2eQTL)
Y = vector("list", 2)  ## Y[[1]] normal; Y[[2]] tumor
Y[[1]] = Exprvarmat[CandidateGenes, Status == "normal"]
Y[[2]] = Exprvarmat[CandidateGenes, Status == "tumor"]
rownames(Y[[1]]) = rownames(Y[[2]]) = NULL
X = vector("list", 2)
X[[1]] = SNPvarmat[CandidateEQTLs, Status == "normal"]
X[[2]] = SNPvarmat[CandidateEQTLs, Status == "tumor"]
rownames(X[[1]]) = rownames(X[[2]]) = NULL

data = list(
  Data = list(
    X = X, Y = Y, Sk = Sk
  ),
  Vars = list(
    Genes = CandidateGenes, eQTLs = CandidateEQTLs,
    n = N, p = Ng, k = Nk
  )
)

saveRDS(data, "./data2/gastric0.01.rds")


seed = as.numeric(Sys.time())
set.seed(seed)
data = readRDS("./script/data2/gastric0.01.rds")
## data$Data$Y = lapply(data$Data$Y, function(D){ 2 * log2(D) })
L2lamax(data$Data$X, data$Data$Y, data$Data$Sk, data$Vars$n, data$Vars$p, data$Vars$k)
gamma = cv.multiRegression(data$Data$X, data$Data$Y, data$Data$Sk, ngamma = 20, nfold = 5, data$Vars$n, data$Vars$p, data$Vars$k)
ifit  = multiRegression(data$Data$X, data$Data$Y, data$Data$Sk, gamma, data$Vars$n, data$Vars$p, data$Vars$k, trans = FALSE)
Xs    = data$Data$X
colnames(Xs[[1]]) = colnames(Xs[[2]]) = NULL
Ys    = data$Data$Y
colnames(Ys[[1]]) = colnames(Ys[[2]]) = NULL
Sk    = data$Data$Sk

cvfit = opt.multiFSSEMiPALM2(Xs = Xs, Ys = Ys, Bs = ifit$Bs, Fs = ifit$Fs, Sk = Sk,
                            sigma2 = ifit$sigma2, nlambda = 20, nrho = 20,
                            p = data$Vars$p, q = data$Vars$k, wt = T)


fit = multiFSSEMiPALM2(Xs = Xs, Ys = Ys, Bs = ifit$Bs, Fs = ifit$Fs, Sk = Sk,
                       sigma2 = ifit$sigma2, lambda = cvfit$lambda, rho = cvfit$rho,
                       Wl = inverseB(ifit$Bs), Wf = flinvB(ifit$Bs),
                       p = data$Vars$p, maxit = 1000, threshold = 1e-5, sparse = T,
                       verbose = T, trans = T, strict = T)

saveRDS(fit, "./data2/gastricfitFSSEM0.01.rds")

####
fit = readRDS("./script/data2/gastricfitFSSEM0.01.rds")
data = readRDS("./script/data2/gastric0.01.rds")

filterDiffNet =  function(fit, data, cutoff = 0.01) {
  Ci1 = rowMeans(qnorm(1 - cutoff, mean = 0, sd = sqrt(fit$sigma2)) / data$Data$Y[[1]])
  Ci2 = rowMeans(qnorm(1 - cutoff, mean = 0, sd = sqrt(fit$sigma2)) / data$Data$Y[[2]])
  BThreshold = quantile(c(abs(fit$Bs[[1]])[fit$Bs[[1]] != 0], abs(fit$Bs[[2]])[fit$Bs[[2]] != 0]), seq(0, 1, by = 0.1))[2 + 1]
  B1 = fit$Bs[[1]]
  B2 = fit$Bs[[2]]
  ## Threshold filter
  for(i in 1:data$Vars$p) {
    B1[,i] = B1[,i] * (abs(B1[,i]) >= max(BThreshold, Ci1[i]))
    B2[,i] = B2[,i] * (abs(B2[,i]) >= max(BThreshold, Ci2[i]))
  }
  DiffB = B2 - B1
  for(i in 1:data$Vars$p) {
    DiffB[,i] = DiffB[,i] * ((abs(DiffB[,i]) >= max(BThreshold, Ci1[i], Ci2[i])) & abs(DiffB[,i]) >= pmin(abs(B1[,i]), abs(B2[,i])))
  }
  list(B1 = B1, B2 = B2, DiffB = DiffB)
}

FilteredB = filterDiffNet(fit, data, cutoff = 0.01)
##Sigma2Threshold = max(max(qnorm(1 - 0.05, mean = 0, sd = sqrt(fit$sigma2)) / abs(data$Data$Y[[1]])), max(qnorm(1 - 0.05, mean = 0, sd = sqrt(fit$sigma2)) / abs(data$Data$Y[[2]])))
##BThreshold = quantile(c(abs(fit$Bs[[1]])[fit$Bs[[1]] != 0], abs(fit$Bs[[2]])[fit$Bs[[2]] != 0]), seq(0, 1, by = 0.1))[2 + 1]
##BThreshold = max(Sigma2Threshold, BThreshold, 1e-4)
##B1 = fit$Bs[[1]] * (abs(fit$Bs[[1]]) > BThreshold)
##B2 = fit$Bs[[2]] * (abs(fit$Bs[[2]]) > BThreshold)
##DiffB = B2 - B1
B1 = FilteredB$B1
B2 = FilteredB$B2
DiffB = FilteredB$DiffB
adjNormalGRN = as.matrix(B1 != 0)
adjTumorGRN  = as.matrix(B2 != 0)
## adjDifferentialGRN  = as.matrix((abs(DiffB) > BThreshold) & (abs(DiffB) >= pmin(abs(B1), abs(B2))))
adjDifferentialGRN = as.matrix(DiffB != 0)
NormalGRN = graph_from_adjacency_matrix(t(adjNormalGRN)) %>% set_vertex_attr("name", value = data$Vars$Genes)
TumorGRN = graph_from_adjacency_matrix(t(adjTumorGRN)) %>% set_vertex_attr("name", value = data$Vars$Genes)
DifferentialGRN  = graph_from_adjacency_matrix(t(adjDifferentialGRN)) %>% set_vertex_attr("name", value = data$Vars$Genes)
NormalGRN = delete.vertices(igraph::simplify(NormalGRN), degree(NormalGRN) == 0)
TumorGRN = delete.vertices(igraph::simplify(TumorGRN), degree(TumorGRN) == 0)
DifferentialGRN  = delete.vertices(igraph::simplify(DifferentialGRN), degree(DifferentialGRN) == 0)

### Entrez ID to Gene Symbol
DiffGene = bitr(
  names(V(DifferentialGRN)),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

DiffGname = bitr(
  names(V(DifferentialGRN)),
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = "org.Hs.eg.db"
)
## remove non protein coding RNA related mapping
mRNAgenes = annotation2Entrez[annotation2Entrez$TypeOfGene == "protein-coding", 2]
DiffGname = DiffGname[DiffGname$SYMBOL %in% mRNAgenes,]
rownames(DiffGname) = DiffGname[,1]


UniGene = bitr(
  data$Vars$Genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

UniGname = bitr(
  data$Vars$Genes,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = "org.Hs.eg.db"
)
UniGname = UniGname[UniGname$SYMBOL %in% mRNAgenes,]
UniGene = UniGene[UniGene$ENSEMBL %in% UniGname$ENSEMBL,]

### old plot
GRNLayout = function(G = NULL) {
  qcut = sort(degree(G), decreasing = T)[10]
  V(G)$color = "blue"
  V(G)$frame.color = "blue"
  V(G)$color[which(degree(G) >= qcut)] = "red"
  V(G)$frame.color[which(degree(G) >= qcut)] = "darkred"
  V(G)$size = sqrt(degree(G)) * 1.5 + 0.1
  E(G)$color = rgb(0, 0, 0, alpha = .1)
  Vnames = rep(NA, length(V(G)))
  Vnames[which(degree(G) >= qcut)] = DiffGname[names(which(degree(G) >= qcut)), 2]
  plot(
    G,
    vertex.label = Vnames,
    layout = layout.fruchterman.reingold,
    edge.arrow.size = 0.1,
    edge.curve = 0.1,
    vertex.label.font = 1,
    vertex.label.cex = 0.5,
    vertex.label.dist = 1,
    vertex.label.color = "black"
  )
}

GRNLayout(NormalGRN)
GRNLayout(TumorGRN)
GRNLayout(DifferentialGRN)


### GeneSet analysis
### GSEA
getGSEAdb = function(db = NULL) {
  f = readLines(db)
  lst = sapply(f, function(x)
    unlist(strsplit(x, "\t", fixed = TRUE)))
  names(lst) = sapply(lst, function(x)
    x[1])
  lst1 = lapply(lst, function(x)
    x[-(1:2)])
  lst2 = lapply(lst, function(x)
    x[2])
  list(gsea = lst1, url = lst2)
}

GSEAC2 =  getGSEAdb("./script/data2/c2.all.v6.2.entrez.gmt")
EntrezGenes = unique(unlist(GSEAC2$gsea))
names(EntrezGenes) = NULL

testGSEA = function(geneset, diffgene, universe) {
  bgset   = setdiff(universe, diffgene)
  pvalue  = lapply(geneset, function(x) {
                              x1 = length(intersect(diffgene, x))
                              x2 = length(setdiff(diffgene, x))
                              y1 = length(intersect(bgset, x))
                              y2 = length(setdiff(bgset, x))
                              c(fisher.test(matrix(c(x1, x2, y1, y2), nrow = 2), alternative = "two.sided")$p.value,
                                ifelse(x1/x2 >= y1/y2, "+", "-"))
                            })
  pvalue
}

## filter Gene set have lung annotation
GAST_related = lapply(GSEAC2$url, function(url) {
  url = as.character(url)
  mdb = GET(url)
  mdb = readHTMLTable(rawToChar(mdb$content), stringsAsFactors = F)
  i = which(mdb[[1]][,1] == "Full description or abstract")
  description = mdb[[1]][i, 2]
  pattern1 = str_detect(description, "gastric relapse | gastric cancer | gastric tumor | gastric adenocarcinoma")
  pattern2 = all(str_detect(description, c("gastric", "cancer | tumor | adenocarcinoma | carcinoma")))
  pattern1 | pattern2
})
GSEA4GAST = GSEAC2$gsea[unlist(GAST_related)]
GSEA4GAST = GSEA4GAST[sapply(GSEA4GAST, function(x){length(intersect(x, unique(UniGene[,2]))) != 0})]
result = testGSEA(GSEA4GAST, unique(DiffGene[DiffGene$ENSEMBL %in% names(which(degree(DifferentialGRN, mode = "all") >= 1)), 2]),
                  unique(UniGene[,2]))

pvalue_gsea  = as.numeric(sapply(result, `[`, 1))
enrich_gsea  = sapply(result, `[`, 2)
enrich_idx   = which(pvalue_gsea < 0.05 & enrich_gsea == "+")
enrichedGS   = data.frame(gs = names(result)[enrich_idx],
                         description = unlist(GSEAC2$url[names(result)[enrich_idx]]),
                         pval = pvalue_gsea[enrich_idx])
rownames(enrichedGS) = NULL
enrichedGS

DiffGname[names(sort(degree(DifferentialGRN, mode = "all"), decreasing = T))[1:10], 2]

for(n in 1:nrow(enrichedGS)) {
  cat(as.character(enrichedGS[n,1]), " & ", toupper(format(signif(enrichedGS[n,3], 2), scientific = T)), " \\\\ \n")
  cat("\\hline \n")
}

