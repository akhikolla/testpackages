### Sub-option 1: Tests
For this option the null hypothesis is: “Genotypes at one locus are independent from genotypes at the other locus”. For a pair of diploid loci, no assumption is made about the gametic phase in double heterozygotes. In particular, it is not inferred assuming one-locus HW equilibrium, as such equilibrium is not assumed anywhere in the formulation of the test. The test is thus one of association between diploid genotypes at both loci, sometimes described as a test of the composite linkage disequilibrium\index{Linkage disequilibrium!composite} [@WeirbkII, pp. 126-128]. For a haploid locus and a diploid one, a test of association between the haploid and diploid genotypes is computed (there is no concern about gametic phase in this case). This makes it easy to test for cyto-nuclear disequilibria\index{Linkage disequilibrium!cyto-nuclear}. For a pair of loci with haploid information, a straightforward test of association of alleles at the two loci is computed.

The default test statistic is now the log likelihood ratio statistic ($G$-test). However one can still perform probability tests (as implemented in earlier versions of Genepop) by using the `GameticDiseqTest=Proba`\index{GameticDiseqTest setting} setting.

For a given pair of loci within one sample, the relevant information is represented by a contingency table looking e.g. like

           GOT2
           1.1  1.3  3.3  1.7  3.7
    EST    _________________________
     1.1   1    1    0    0    1      3
     1.2   16   6    1    3    2      28
           _________________________
           17   7    1    3    3      31

for two diploid loci (`1.1`, etc., are the diploid genotypes at each locus). Contingency tables are created for all pairs of loci in each sample, then a $G$ test or a probability test for each table is computed for each table using the Markov chain algorithm of @RaymondR95evol. The number of switches\index{Markov chain algorithms!switches} of the algorithm is given for each table analyzed.[^15]

### Output
Results are stored in the file *yourdata*`.DIS`. Three intractable situations are indicated: empty tables (“No data”), table with one row or one column only (“No contingency table”), and tables for which all rows or all columns marginal sums are 1 (“No information”). For each locus pair within each sample, the unbiased estimate of the P-value is indicated, as well as the standard error. Next, a global test (Fisher’s method) for each pair of loci is performed across samples.\index{Combination of different tests}

See also the next section for analysis of a single table.
