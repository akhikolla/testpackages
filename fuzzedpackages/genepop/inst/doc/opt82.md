### Sub-option 2: Diploidisation of haploid data
\index{Haploid data!to haploid}
This sub-option “diploidizes” haploid loci. For example, the line\
`popul 1, 01 02 10 00`\
of an haploid dataset with 4 loci, will become\
`popul 1, 0101 0202 1010 0000`.\
Only haploid data are thus modified in a mixed haploid/diploid data file. The new file is named `D`*yourdata*.[^24]

Note that there may no longer be any need for this option for further analyses with Genepop (except perhaps as a preliminary to file conversions, option 7), since Genepop 4.0 now perform analyses on haploid data without such prior “diploidization” (don’t forget the `EstimationPloidy=Haploid` setting).
