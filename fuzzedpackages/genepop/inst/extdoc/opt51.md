### Sub-option 1: Allele and genotype frequencies
This option provides basic information on the data set. The output file is saved in the file *yourdata*`.INF`. For each locus in each sample, several variables are calculated:

-   allele frequencies.

-   observed and expected genotype proportions.

-   $F_\mathrm{IS}$ estimates for each allele following @WeirC84.

-   global estimate of $F_\mathrm{IS}$ over alleles according to @WeirC84 (W&C) and @RobertsonH84 (R&H).

-   observed and “expected” numbers of homozygotes and heterozygotes. “Expected” here means the expected numbers, conditional on observed allelic counts, under HW equilibrium; the difference from naive products of observed allele frequencies is sometimes called Levene’s correction, after\index{Levene's correction} @Levene49.

-   the genotypic matrix.

A table of allele frequencies for each locus and for each sample is also computed.
