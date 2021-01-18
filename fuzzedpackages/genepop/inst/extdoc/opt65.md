### Sub-option 5: isolation by distance between individuals

This option allows analysis of isolation by distance between pairs of individuals. It provides estimates of “neighborhood size”,\index{Neighborhood size|see{$D\sigma^2$ estimation}} or more precisely of $D\sigma^2$, the product of population density and axial mean square parent-offspring distance, derived from the slope of the regression of pairwise genetic statistics against geographical distance or log(distance) in linear or two-dimensional habitats, respectively. More details are described in @Rousset00 ($\hat{a}$ statistic), @LebloisER03 (bootstrap confidence intervals) and @WattsX07 ($\hat{e}$ statistic). For haploid data, a proxy for the $\hat{a}$ statistic has been introduced in version 4.1.

The position of individuals must be specified as two coordinates standing for their name (i.e. before the comma on the line for each individual), and since each individual is considered as a sample, it must be separated by a `Pop`. An example of such input file is given below: The first individual is located at the point $x = 0.0$, $y = 15.0$ (showing that the decimal separator is a period), the second at the point $x = 0$, $y =30$, etc. This example also shows that *individual identifiers can be added after these coordinates*.

     Title line: A really too small data set
     ADH Locus 1
     ADH #2
     ADH three
     ADH-4
     ADH-5
     Pop
     0.0 15.0,  0201 0303 0102 0302 1011
     Pop
     0 30 Second indiv,  0202 0301 0102 0303 1111
     Pop
     0 45,  0102 0401 0202 0102 1010
     Pop
     0 60,  0103 0202 0101 0202 1011
     Pop
     0 75,  0203 0204 0101 0102 1010
     POP
     15 15,      0102 0202 0201 0405 0807
     Pop
     15 30,      0102 0201 0201 0405 0307
     Pop
     15 45,      0201 0203 0101 0505 0402
     Pop
     15 60,      0201 0303 0301 0303 0603
     Pop
     15 75,      0101 0201 0301 0505 0807

**Missing information** arises when there is no genetic estimate (if a pair of individuals has no genotypes for the same locus, for example), or when geographic distance is zero and log(distance) is used. Genepop will correctly handle such missing information until it comes to the point where regression cannot be computed or there are not several loci to bootstrap over.

Options to be described within option 6.5 are: $\hat{a}$ or $\hat{e}$ pairwise statistics (for diploid data); log transformation for geographic distances; minimal geographic distance; coverage probability of confidence interval; testing a given value of the slope; Mantel test settings; conversion to genetic distance matrix in Phylip format. Allele-size based analogues of $\hat{a}$ or $\hat{e}$ can be defined, but they should perform very poorly [@LebloisER03; @Rousset07w], so such an analysis has been purposely disabled.

**Pairwise statistics for diploid data**: They are selected by the setting `IsolationStatistic=a`\index{IsolationStatistic setting}\index{Dsigma2@$D\sigma^2$
estimation!a@$\hat{a}$ statistic} or `=e`, or at runtime (in batch mode, the default is $\hat{a}$). The $\hat{e}$ statistic is asymptotically biased in contrast to $\hat{a}$, but has lower variance. The bias of the $\hat{e}$-based slope is higher the more limited dispersal is, so it performs less well in the lower range of observed dispersal among various species. Confidence intervals are also biased [@LebloisER03; @WattsX07], being too short in the direction of low $D\sigma^2$ values, and on the contrary conservative in the direction of low $D\sigma^2$ values. Based on the simulation results of @WattsX07, a provisional advice is to run analyses with both statistics, and to derive an upper bound for the $D\sigma^2$ confidence interval (CI), hence the lower bound for the regression slope, from $\hat{e}$ (which has CI shorter than $\hat{a}$, though still conservative) and the other $D\sigma^2$ bound, hence the upper bound for the regression slope, from $\hat{a}$ (which has too short CI, but less biased than the $\hat{e}$ CI). When the $\hat{e}$-based $D\sigma^2$ estimate is below 2500 (linear habitat) or 4 (two-dimensional habitat) it is suggested to derive both bounds from $\hat{a}$.

 Note that $\hat{e}$ is essentially Loiselle’s statistic\index{Dsigma2@$D\sigma^2$ estimation!Loiselle's statistic} [@LoiselleSNG95], which use in this context has previously been advocated by e.g. @VekemansH04.

For **haploid data** (i.e. `EstimationPloidy=Haploid`) the denominators of the $\hat{a}$ and $\hat{e}$ statistics cannot be computed. Ideally the denominator should be the gene diversity among individuals that would compete for the same position, as could be estimated from “group” data. As a reasonable first substitute, Genepop uses a single estimate of gene diversity (from the total sample and for each locus) to compute the denominators for all pairs of individuals. This amount to assume that overall differentiation in the population is weak.

**Log transformation for geographic distances**: This transformation is required for estimation of $D\sigma^2$ when dispersal occurs over a surface rather than over a linear habitat. It is the default option in batch mode. It can be turned on and off by the setting `GeographicScale=Log`\index{GeographicScale setting} or `=Linear` or equivalently by `Geometry=2D` or `=1D`.\index{Geometry setting}

**Coverage probability of confidence interval** This is the target probability that the confidence interval contains the parameter value. The usage is to compute intervals with 95% coverage and equal 2.5% tails, and this is the default coverage in Genepop. This can be changed by the setting `CIcoverage`, e.g. `CIcoverage=0.99` will compute interval with target probabilities 0.5% that either the confidence interval is too low or too high (an unrealistically large number of loci may be necessary to achieve the latter precision).\index{CIcoverage setting}\index{Confidence intervals}

**Minimal and maximal geographic distances:** As discussed in @Rousset97, samples at small geographic distances are not expected to follow the simple theory of the regression method, so the program asks for a minimum geographical distance. Only pairwise comparisons of samples at larger distances are used to estimate the regression coefficient (all pairs are used for the Mantel test). The minimal distance may be specified by the setting `MinimalDistance=`*value*\index{MinimalDistance setting} or at runtime. This being said, it is wise to include all pairs in the estimation as no substantial bias is expected, and this avoids uncontrolled hacking of the data. Thus, the suggested minimal distance here is any distance large enough to exclude only pairs at zero geographical distance. Only non-negative values are accepted, and the default in batch mode is 0.0001.

There is also a setting `MaximalDistance=`*value*.\index{MaximalDistance setting} This should not be abused, and is (therefore) available only through the settings file, not as a runtime option.

**Testing a given value of the slope** The setting `testPoint=0.00123` (say) returns the unidirectional P-value for a specific value of the slope, using the ABC bootstrap method. This is the reciprocal of a confidence interval computation: confidence intervals evaluate parameter values corresponding to given error levels, say the 0.025 and 0.975 unidirectional levels for a 95% bidirectional CI, while this option evaluates the unidirectional P-value associated with a given parameter value.

**Mantel test:**\index{Mantel test} The Mantel test is implemented. See Section \@ref(mantel-test) for limitations of this test. In the present context this is an exact test of the null hypothesis that there in no spatial correlation between genetic samples.

Up to version 4.3 Genepop implemented only a Mantel test based on the rank correlation. It now also implements, and performs by default, Mantel tests based on the regression coefficient for the “genetic distance” statistic used to quantify isolation by distance. The latter tests should generally be more congruent with the confidence intervals based on the same distances than the rank-based tests are. The rank test can now be performed by using the setting `MantelRankTest=` (no *value* needed).\index{MantelRankTest setting}

Ideally the confidence interval for the slope should contain zero if and only if the Mantel test is non-significant. Some exceptions may occur as the bootstrap method is only approximate, but such exceptions appear to be rare. Exceptions may more commonly occur when the bootstrap is based on the regression of genetic “distance” and geographic distance over a selected range of the latter.

The number of permutations may be specified by the setting `MantelPermutations=`*value*,\index{MantelPermutations setting} or else at runtime. In batch mode, if no such value has been given the default behaviour is not to perform the test.

**Export genetic distance matrix in Phylip format**.\index{Phylip package|see PhylipMatrix} This option is activated by the setting `PhylipMatrix=` (no *value* needed).\index{PhylipMatrix setting} It may be useful, if you wish to use Phylip, to draw a tree based on genetic distances. A constant is added to all values if necessary so that all resulting distances are positive. Output is written in the file *yourdata*`.PMA`. No further estimation or testing is done, so the name of the groups/individuals does not need to be their spatial coordinates.

Except for this export option, output files are:

-   the *yourdata*`.ISO` output file, containing (i) a genetic distance ($\hat{a}$ or $\hat{e}$) half-matrix and a geographic (log-)distance half-matrix; missing information is reported as ‘`-`’; (ii) regression estimates and bootstrap confidence intervals; (iii) the result of testing a slope value (using `testPoint`); (iv) results of a Mantel test for evidence of isolation by distance, if requested; (v) a bootstrap interval for the intercept. The order of elements in the half-matrices is:

               1     2     3
         2     x
         3     x     x
         4     x     x     x

-   a *yourdata*`.MIG` output file, containing the same genetic and geographic distances as in the `ISO` file, but with more digits, and without estimation or test results. This file was formerly useful as input for the Isolde program (see “Former option 5 of Genepop”, below), and is a bit redundant now.

-   a *yourdata*`.GRA` output file, where again the genetic and geographic distances are reported, now as $(x,y)$ coordinates for each pair of individuals (one per line). This is useful e.g. for importing the output into programs with good graphics. Pairs with missing values (either $x$ or $y$) are not reported in this file.
