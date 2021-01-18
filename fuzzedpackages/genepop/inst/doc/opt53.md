### Sub-option 3: Allele size-based gene diversities and $\rho_{\mathrm{IS}}$

Option 5.3 is analogous to option 5.2. It computes measures of diversity
based on allele size, namely mean squared allele size differences within
individuals (`MSDintra`), and among individuals within samples
(`MSDinter`), per locus per sample, and averaged over samples or over
loci. Corresponding $\rho_\mathrm{IS}$ (the $F_\mathrm{IS}$ analogue, see Section \@ref(rho-stats)) estimates are also computed. Allele size is the allele name unless it has been given through the `AlleleSizes` setting.

For haploid data, only the mean squared difference `MSDinter` among individuals is computed. Multilocus estimates ignore haploid loci, or on the contrary ignore diploid loci if the setting `EstimationPloidy=Haploid` is used. Single-locus estimates are computed for both haploid and diploid loci irrespective of this setting.

The output is saved in the file *yourdata*`.MSD`.
