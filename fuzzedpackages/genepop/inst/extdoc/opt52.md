### Sub-option 2: Identity-based gene diversities and $F_\mathrm{IS}$

This option takes the observed frequencies of identical pairs of genes as estimates ($Q$) of corresponding probabilities of identity ($Q$) and then simply computes diversities as $1-Q$: gene diversity within individuals (`1-Qintra`), and among individuals within samples (`1-Qinter`), per locus per sample, and averaged over samples or over loci. One-locus $F_\mathrm{IS}$ estimates are also computed in a way consistent with @WeirC84. No estimate is given when no information is available (e.g. no estimate of diversity between individuals within a sample when only one individual has been genotyped).

For haploid data\index{Haploid data}, only the gene diversity among individuals is computed. Multilocus estimates ignore haploid loci, or on the contrary ignore diploid loci if the setting `EstimationPloidy=Haploid` is used. Single-locus estimates are computed for both haploid and diploid loci irrespective of this setting.

The output is saved in the file *yourdata*`.DIV`.
