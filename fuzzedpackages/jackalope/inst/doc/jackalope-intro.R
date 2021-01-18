## ----setup, echo = FALSE, message = FALSE-------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
library(jackalope)
library(ape)
set.seed(65456156)

## ----examples-read-assembly-for-show, eval = FALSE----------------------------
#  ref <- read_fasta("dmel-6.27.fasta.gz", cut_names = TRUE)
#  ref$filter_chroms(1e6, method = "size")
#  ref$rm_chroms("Y")
#  ref$merge_chroms(c("2L", "2R"))
#  ref$merge_chroms(c("3L", "3R"))
#  names <- ref$chrom_names()
#  names[grepl("^2", names)] <- "2"
#  names[grepl("^3", names)] <- "3"
#  ref$set_names(names)

## ----examples-create-assembly-------------------------------------------------
ref <- create_genome(n_chroms = 4,
                     len_mean = 1e3 / 4,
                     len_sd = 10)
ref$set_names(c(2:4, "X"))

## ----examples-print-assembly, echo = FALSE------------------------------------
print(ref)

## ----examples-mevo-objects----------------------------------------------------
sub_rate <- evo_rates$subs[evo_rates$species == "Drosophila melanogaster"]
indel_rate <- evo_rates$indels[evo_rates$species == "Drosophila melanogaster"]
# Because both are in units of 10^-10 events per site per generation:
sub_rate <- sub_rate * 1e-10
indel_rate <- indel_rate * 1e-10

sub <- sub_JC69(lambda = sub_rate, mu = NULL)
ins <- indels(rate = indel_rate, max_length = 60,a = 1.60)
del <- indels(rate = indel_rate, max_length = 60, a = 1.51)

theta <- evo_rates[evo_rates$species == "Drosophila melanogaster","theta_s"]

# Originally in units of 1e6 individuals
N0 <- evo_rates[evo_rates$species == "Drosophila melanogaster", "Ne"] * 1e6

## ----examples-reads-for-assembly-pacbio, eval = FALSE-------------------------
#  pacbio(ref, out_prefix = "pacbio", n_reads = 2 * 500e3)

## ----examples-reads-for-assembly-hybrid, eval = FALSE-------------------------
#  pacbio(ref, out_prefix = "pacbio", n_reads = 500e3)
#  illumina(ref, out_prefix = "illumina", n_reads = 500e6, paired = TRUE,
#           read_length = 100)

## ----examples-reads-for-assembly-illumina, eval = FALSE-----------------------
#  illumina(ref, out_prefix = "ill_pe", n_reads = 500e6, paired = TRUE,
#           read_length = 100)
#  illumina(ref, out_prefix = "ill_mp", seq_sys = "MSv3",
#           read_length = 250, n_reads = 50e6, matepair = TRUE,
#           frag_mean = 3000, frag_sd = 500)

## ----examples-assembly-diploid-haplotypes-------------------------------------
haps <- create_haplotypes(ref, haps_theta(theta = theta, n_haps = 2), 
                          sub, ins, del)
haps$set_names(c("A", "B"))

## ----examples-assembly-diploid-haplotypes-print, echo = FALSE-----------------
print(haps)

## ----examples-reads-for-assembly-hybrid-diploid, eval = FALSE-----------------
#  pacbio(haps, out_prefix = "pacbio", n_reads = 500e3)
#  illumina(haps, out_prefix = "illumina", n_reads = 500e6, paired = TRUE,
#           read_length = 100)

## ----examples-reads-for-assembly-diploid-output, eval = FALSE-----------------
#  write_fasta(haps, "haps")
#  write_vcf(haps, "haps", sample_matrix = cbind(1, 2))

## ----examples-divergence-scrm, eval = FALSE-----------------------------------
#  library(scrm)
#  # Function to run scrm for one chromosome and format output
#  one_chrom <- function(.size) {
#      ssites <- scrm(sprintf("10 1 -t %.4f -r 1 %i -I 2 5 5 100", theta * 4 * N0, .size))
#      return(ssites$seg_sites[[1]])
#  }
#  ssites <- list(seg_sites = lapply(ref$sizes(), one_chrom))

## ----examples-divergence-create, eval = FALSE---------------------------------
#  haps <- create_haplotypes(ref, haps_ssites(ssites), sub, ins, del)

## ----examples-divergence-create-do, echo = FALSE------------------------------
# For the purposes of the vignette, I'm just going to use `haps_theta` so I
# don't have to (1) load scrm to build the vignette or (2) keep a relatively
# large internal data file to store the scrm output
# theta of 0.4 gives the same # mutations as using `ssites` (~ 2,500)
set.seed(7809534)
haps <- create_haplotypes(ref, haps_theta(theta = 0.4, n_haps = 10), sub, ins, del)

## ----examples-divergence-create-print, echo = FALSE---------------------------
print(haps)

## ----examples-divergence-write-vcf, eval = FALSE------------------------------
#  write_vcf(haps, "haplotypes")

## ----examples-divergence-illumina-pool, eval = FALSE--------------------------
#  illumina(haps, out_prefix = "haps_illumina", n_reads = 500e6, paired = TRUE,
#           read_length = 100, barcodes = c(rep("AACCGCGG", 5),
#                                           rep("GGTTATAA", 5)))

## ----examples-divergence-illumina-individual, eval = FALSE--------------------
#  illumina(haps, out_prefix = "haps_illumina", n_reads = 500e6, paired = TRUE,
#           read_length = 100, sep_files = TRUE)

## ----examples-phylogeny-tree--------------------------------------------------
tree <- rcoal(10)
tree$edge.length <- 4 * N0 * tree$edge.length / max(node.depth.edgelength(tree))

## ----examples-phylogeny-tree-create-for-show----------------------------------
haps <- create_haplotypes(ref, haps_phylo(tree), sub, ins, del)

## ----examples-phylogeny-tree-haplotypes-print, echo = FALSE-------------------
print(haps)

## ----examples-phylogeny-tree-illumina, eval = FALSE---------------------------
#  haplotype_barcodes <- c("CTAGCTTG", "TCGATCCA", "ATACCAAG", "GCGTTGGA",
#                          "CTTCACGG", "TCCTGTAA", "CCTCGGTA", "TTCTAACG",
#                          "CGCTCGTG", "TATCTACA")
#  illumina(haps, out_prefix = "phylo_tree", seq_sys = "MSv3",
#           paired = TRUE, read_length = 250, n_reads = 50e6,
#           barcodes = haplotype_barcodes)
#  ape::write.tree(tree, "true.tree")

## ----examples-phylogeny-gtrees-scrm, eval = FALSE-----------------------------
#  # Run scrm for one chromosome size:
#  one_chrom <- function(.size) {
#      sims <- scrm(
#          paste("24 1",
#                # Output gene trees:
#                "-T",
#                # Recombination:
#                "-r 1", .size,
#                # 3 species with no ongoing migration:
#                "-I 3 8 8 8 0",
#                # Species 2 derived from 1 at time 1.0:
#                "-ej 1.0 2 1",
#                # Species 3 derived from 2 at time 0.5:
#                "-ej 0.5 3 2"
#          ))
#      trees <- sims$trees[[1]]
#      # scrm outputs branch lengths in units of 4*N0 generations, but we want just
#      # generations:
#      adjust_tree <- function(.p) {
#          # Read to phylo object and adjust branch lengths:
#          .tr <- read.tree(text = .p)
#          .tr$edge.length <- .tr$edge.length * 4 * N0
#          # "prefix" from `.p` showing how large the region this gene tree refers to is
#          prefix <- paste0(strsplit(.p, "\\]")[[1]][1], "]")
#          # Put back together into NEWICK text
#          return(paste0(prefix, write.tree(.tr)))
#      }
#      trees <- sapply(trees, adjust_tree)
#      names(trees) <- NULL
#      return(trees)
#  }
#  # For all chromosomes:
#  gtrees <- list(trees = lapply(ref$sizes(), one_chrom))

## ----examples-phylogeny-gtrees-write-true-gtrees, eval = FALSE----------------
#  write_gtrees(haps_gtrees(gtrees), "gtrees")

## ----examples-phylogeny-gtrees-create-haplotypes, eval = FALSE----------------
#  haps <- create_haplotypes(ref, haps_gtrees(gtrees),
#                            sub, ins, del)

## ----examples-phylogeny-gtrees-create-do, echo = FALSE------------------------
# For the purposes of the vignette, I'm just going to use `haps_theta` so I
# don't have to (1) load scrm to build the vignette or (2) keep a relatively
# large internal data file to store the scrm output
# theta of 0.125 gives the same # mutations as using `gtrees` (~2,600)
set.seed(7809534)
haps <- create_haplotypes(ref, haps_theta(theta = 0.125, n_haps = 24), sub, ins, del)

## ----examples-phylogeny-gtrees-haplotypes-print, echo = FALSE-----------------
print(haps)

## ----examples-phylogeny-write-vcf, eval = FALSE-------------------------------
#  write_vcf(haps, out_prefix = "hap_gtrees",
#            sample_matrix = matrix(1:haps$n_haps(), ncol = 2, byrow = TRUE))

## ----examples-phylogeny-gtrees-illumina, eval = FALSE-------------------------
#  # 2 of each barcode bc it's diploid
#  haplotype_barcodes <- rep(c("TCGCCTTA", "CTAGTACG", "TTCTGCCT", "GCTCAGGA", "AGGAGTCC",
#                              "CATGCCTA", "GTAGAGAG", "CCTCTCTG", "AGCGTAGC", "CAGCCTCG",
#                              "TGCCTCTT", "TCCTCTAC"), each = 2)
#  illumina(haps, out_prefix = "phylo_gtrees", seq_sys = "MSv3",
#           read_length = 250, n_reads = 50e6, paired = TRUE,
#           barcodes = haplotype_barcodes)

