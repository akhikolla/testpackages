## ----knitr_options, include=FALSE---------------------------------------------
library(knitr)
opts_chunk$set(fig.width=7, fig.height=4.5,
               dev.args=list(pointsize=16))

## ----set.seed, include=FALSE--------------------------------------------------
set.seed(80607574)

## ----create_parent------------------------------------------------------------
library(simcross)
p1 <- create_parent(L=100, allele=1)
p2 <- create_parent(L=100, allele=2)
f1 <- create_parent(L=100, allele=1:2)

## ----sim_f2-------------------------------------------------------------------
(f2 <- cross(f1, f1))

## ----sim_f2_nointerference----------------------------------------------------
f2 <- cross(f1, f1, m=0, obligate_chiasma=TRUE)

## ----AILped_head--------------------------------------------------------------
head(AILped, n=10)

## ----check_AILped-------------------------------------------------------------
check_pedigree(AILped)

## ----sim_ril_ped--------------------------------------------------------------
(ril <- sim_ril_pedigree(ngen=4, parents=1:4))

## ----sim_4way_ped-------------------------------------------------------------
(fourway <- sim_4way_pedigree(ngen=2, nsibs=c(3, 3)))

## ----sim_ail_ped--------------------------------------------------------------
ailped <- sim_ail_pedigree(ngen=12, npairs=100, nkids_per=5)
nrow(ailped)
table(ailped$gen)

## ----sim_do_ped---------------------------------------------------------------
doped <- sim_do_pedigree(ngen=12)
nrow(doped)
table(do=doped$do, gen=doped$gen)

## ----sim_ail_fully------------------------------------------------------------
ailped <- sim_ail_pedigree(ngen=8, npairs=30, nkids_per=5)
xodat <- sim_from_pedigree(ailped, L=100)

## ----last_ail_individual------------------------------------------------------
xodat[[length(xodat)]]

## ----plot_ave_breakpoints-----------------------------------------------------
n_breakpoints <- sapply(xodat, function(a) sum(sapply(a, function(b) length(b$alleles)-1)))
ave_breakpoints <- tapply(n_breakpoints, ailped$gen, mean)
gen <- as.numeric(names(ave_breakpoints))
plot(gen, ave_breakpoints,
     xlab="Generation", ylab="Average no. breakpoints", las=1,
     pch=21, bg="Orchid", main="AIL with 30 breeding pairs")

## ----where_het----------------------------------------------------------------
where_het(xodat[[length(xodat)]])

## ----plot_prop_het------------------------------------------------------------
prop_het <- sapply(lapply(xodat, where_het), function(a) sum(a[,2]-a[,1])/100)
ave_prop_het <- tapply(prop_het, ailped$gen, mean)
plot(gen, ave_prop_het,
     xlab="Generation", ylab="Average proportion heterozygous", las=1,
     pch=21, bg="Orchid", main="AIL with 30 breeding pairs")
abline(h=0.5, lty=2)

## ----reset_seed, include=FALSE------------------------------------------------
set.seed(28998542)

## ----ail_few_pairs------------------------------------------------------------
ailped2 <- sim_ail_pedigree(ngen=8, npairs=3, nkids_per=50)
xodat2 <- sim_from_pedigree(ailped2, L=100)
prop_het2 <- sapply(lapply(xodat2, where_het), function(a) sum(a[,2]-a[,1])/100)
ave_prop_het2 <- tapply(prop_het2, ailped2$gen, mean)
gen2 <- as.numeric(names(ave_prop_het2))
plot(gen2, ave_prop_het2,
     xlab="Generation", ylab="Average proportion heterozygous", las=1,
     pch=21, bg="Orchid", main="AIL with 3 breeding pairs")
abline(h=0.5, lty=2)

## ----get_geno-----------------------------------------------------------------
g30 <- get_geno(xodat, 30)
tail(g30)

## ----construct_map------------------------------------------------------------
map <- seq(0, 100, by=10)
names(map) <- paste0("m", map)
map

## ----convert2geno-------------------------------------------------------------
geno <- convert2geno(xodat, map)
geno[nrow(geno)-(4:0), ]

## ----convert2geno_8wayril-----------------------------------------------------
rilped <- sim_ril_pedigree(ngen=6, parents=1:8)
dat <- sim_from_pedigree(rilped, L=100)
geno <- convert2geno(dat, map)
geno[nrow(geno)-(1:0),,]

## ----sim_founder_alleles------------------------------------------------------
fg <- matrix(sample(1:2, 8*length(map), replace=TRUE), nrow=8)

## ----convert2geno_8wayril_snps------------------------------------------------
snpgeno <- convert2geno(dat, map, fg)
snpgeno[nrow(snpgeno)-(1:0),]

## ----ail_mult_chr-------------------------------------------------------------
ailped3 <- sim_ail_pedigree(ngen=12, npairs=30, nkids_per=3)
xodat3 <- sim_from_pedigree(ailped3, c("1"=100, "2"=75, "X"=100), "X")
xodat3alt <- sim_from_pedigree(ailped3, c("1"=100, "2"=75, "X"=100),
                               c(FALSE, FALSE, TRUE))

## ----ail_mult_chr_no_X--------------------------------------------------------
xodat3_noX <- sim_from_pedigree(ailped3, c("1"=100, "2"=75, "3"=100), "")

## ----construct_marker_map-----------------------------------------------------
map <- list("1"=seq(0, 100, by=5),
            "2"=seq(0, 75, by=5),
            "X"=seq(0, 100, by=5))
for(i in seq(along=map))
    names(map[[i]]) <- paste0("m", names(map)[i], "_", map[[i]])

## ----convert2geno_multchr-----------------------------------------------------
geno <- convert2geno(xodat3, map)

## ----do_founder_geno_mult_chr-------------------------------------------------
fg <- vector("list", length(map))
for(i in seq(along=map))
    fg[[i]] <- matrix(sample(1:2, 8*length(map[[i]]), replace=TRUE), nrow=8)

## ----sim_do_mult_chr----------------------------------------------------------
doped <- sim_do_pedigree(ngen=4, nkids_per=1)
xodat <- sim_from_pedigree(doped, L=c("1"=100, "2"=75, "X"=100), "X")
xodat <- collapse_do_alleles(xodat)
geno <- convert2geno(xodat, map, fg)

