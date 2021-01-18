Version 1.0.2: July 17, 2020
-------------------------------------------------------------------------------

Diversity:

+ Fixed a bug in `plotDiversityTest` that caused all values of `q` to appear on
  the plot rather than just the specified one.

Gene Usage:

+ Fixed a major bug in the single-cell mode of `groupGenes` where the `v_call`  
  column was being used in instead of the `j_call` column for J gene grouping.
+ Added support for TCR genes to `groupGenes`.
+ Changed the `only_igh` argument of `groupGenes` to `only_heavy`.


Version 1.0.1:  May 8, 2020
-------------------------------------------------------------------------------

Backwards Incompatible Changes:

+ Changed default expected data format from the Change-O data format to the
  AIRR Rearrangement standard. For example: where functions used the column 
  name `V_CALL` (Change-O) as the default to identify the field that stored 
  the V gene calls, they now use `v_call` (AIRR). That means, scripts that 
  relied on default values (previously, `v_call="V_CALL"`), will now fail if 
  calls to the functions are not updated to reflect the correct value for the 
  data. If data are in the Change-O format, the current default value 
  `v_call="v_call"` will fail to identify the column with the V gene calls
  as the column `v_call` doesn't exist. In this case, `v_call="V_CALL"` needs 
  to be specified in the function call.
+ `ExampleDb` converted to the AIRR Rearrangement standard and examples updated 
  accordingly. The legacy Change-O version is available as `ExampleDbChangeo`.
+ For consistency with the style of the new data format default, other field 
  names have been updated to use the same capitalization. This change affects:
     - amino acid physicochemical properties (e.g. `GRAVY` to `gravy`); 
     - `countGenes`, `countClones` (e.g., `SEQ_COUNT` to `seq_count`)
     - `estimateAbundance` (e.g., `RANK` to `rank`)
     - `groupGenes` (e.g., `VJ_GROUP` to `vj_group`)
     - `collapseDuplicates` and `makeChangeoClone` (e.g., `SEQUENCE_ID` to 
       `sequence_id`, `COLLAPSE_COUNT` to `collapse_count`)
     - lineage tree functions (`summarizeTrees`, `getPathLengths`, `getMRCA`, 
       `tableEdges`, `testEdges`) also return columns in lower case (e.g., 
       `parent`, `child`, `outdegree`, `steps`, `annotation`, `pvalue`)
+ `IG_COLOR` names converted to official C region identifiers 
   (IGHA, IGHD, IGHE, IGHG, IGHM, IGHK, IGHL).

General:

+ License changed to AGPL-3.
+ `baseTheme` looks is now consistent across `sizing` options.
+ `cpuCount` will now return `1` if the core count cannot be determined.
+ Fixed a bug in `padSeqEnds` wherein the `pad_char` argument was being 
  ignored.

Diversity:

+ Fixed documentation error in diversity vignette for viewing test results.
+ `estimateAbundance` slot `clone_by` now contains the name of the column
   with the clonal group identifier, as specified in the function call. For 
   example, if the function was called with `clone="clone_id"`, 
   then the `clone_by` slot will be `clone_id`.

Lineage:

+ Renamed the `buildPhylipLineage` arguments `vcall`, `jcall` and 
  `dnapars_exec` to `v_call`, `j_call` and `phylip_exec`, respectively.


Version 0.3.0:  July 17, 2019
-------------------------------------------------------------------------------

Deprecated:

+ `rarefyDiversity` is deprecated in favor of `alphaDiversity`, which includes
  the same functionality.
+ `testDiversity` is deprecated. The test calculation have been added to the 
   normal output of `alphaDiversity`.

General:

+ Added `ape` and `tibble` dependencies.

Lineage:

+ Added `readIgphyml` to read in IgPhyML output and `combineIgphyml` to 
  combine parameter estimates across samples.
+ Added `graphToPhylo` and `phyloToGraph` to allow conversion between 
  graph and phylo formats.

Diversity:

+ Fixed a bug in `estimateAbundance` where setting the `clone` column to a 
  non-default value produced an error.
+ Added rarefaction options to `estimateAbundance` through the `min_n`, 
  `max_n`, and `uniform` arguments.
+ Moved the rarefaction calculation for the diversity functions into 
  `estimateAbundance`. `alphaDiversity` will call `estimateAbundance` for 
  bootstrapping if not provided an existing `AbundanceCurve` object.
+ Restructured the `DiversityCurve` and `AbundanceCurve` objects to accomodate
  the new diversity methods.

Gene Usage:

+ `groupGenes` now supports grouping by V gene, J gene, and junction length 
   (`junc_len`) as well, in addition to grouping by V gene and J gene without 
   junction length. Also added support for single-cell input data with the addition
   of new arguments `cell_id`, `locus`, and `only_igh`.


Version 0.2.11:  September 12, 2018
-------------------------------------------------------------------------------

General:

+ Added `nonsquareDist` function to calculate the non-square distance matrix of 
  sequences.
+ Exported some internal utility functions to make them available to dependent 
  packages: `progressBar`, `baseTheme`, `checkColumns` and `cpuCount`.

Diversity:

+ `estimateAbundance`, and `plotAbundanceCurve`, will now allow `group=NULL`
  to be specified to performance abundance calculations on ungrouped data.

Gene Usage:

+ Added `fill` argument to `countGenes`. When set `TRUE` this adds zeroes 
  to the `group` pairs that do not exist in the data.
+ Added new function `groupGenes` to group sequences sharing same V and J gene.
  
Toplogy Analysis:

+ Fixed a bug in tableEdges causing it to fail when no parent/child 
  relationships exist when specifying `indirect=TRUE`.
+ `makeChangeoClone` will now issue an error and terminate, instead of 
  continuing with a warning, when all sequences are not the same length.
  

Version 0.2.10:  March 30, 2018
-------------------------------------------------------------------------------

General:

+ Fixed a bug in `IPUAC_AA` wherein X was not properly matching against Q.
+ Changed behavior in `getAAMatrix` to treat * (stop codon) as a mismatch.


Version 0.2.9:  March 21, 2018
-------------------------------------------------------------------------------

General:

+ Added explicit type casting for known columns to `readChangeoDb`.
+ Added the `padSeqEnds` function which pads sequences with Ns to make
  then equal in length.
+ Added verification of unique sequence IDs to `collapseDuplicates`.

Diversity:

+ Added the `uniform` argument to `rarefyDiversity` allowing users to toggle
  uniform vs non-uniform sampling.
+ Renamed `plotAbundance` to `plotAbundanceCurve`.
+ Changed `estimateAbundance` return object from a data.frame to a new 
  `AbundanceCurve` custom class.
+ Set default `plot` call for `AbundanceCurve` to `plotAbundanceCurve`.
+ Added the `annotate` argument from `plotDiversityCurve` to 
  `plotAbundanceCurve`.
+ Added the `score` argument to  `plotDiversityCurve` to toggle between 
  plotting diversity or evenness.
+ Added the function `plotDiversityTest` to generate a simple plot of
  `DiversityTest` object summaries.

Gene Usage:

+ Added the `omit_nl` argument to `getAllele`, `getGene` and `getFamily` to
  allow optional filtering of non-localized (NL) genes.

Lineage:

+ Fixed a bug in `makeChangeoClone` preventing it from interpreting the `id`
  argument correctly.
+ Added the `pad_end` argument to `makeChangeoClone` to allow automatic 
  padding of ends to make sequences the same length.


Version 0.2.8:  September 21, 2017
-------------------------------------------------------------------------------

General:

+ Updated Rcpp dependency to 0.12.12.
+ Added `dry` argument to `collapseDuplicates` which will annotate duplicate 
  sequences but not remove them when set to `TRUE`.
+ Fixed a bug where `collapseDuplicates` was returning one sequence if all 
  sequences were considered ambiguous.

Lineage:

+ Added ability to change masking character and distance matrix used in 
  `makeChangeoClone` and `buildPhylipLineage` for purposes of (optionally) 
  treating indels as mismatches.
+ Fixed a bug in `buildPhylipLineage` when PHYLIP doesn't generate inferred
  sequences and has only one block.


Version 0.2.7:  June 12, 2017
-------------------------------------------------------------------------------

General:

+ Fixed a bug in `readChangeoDb` causing the `select` argument to do nothing.
+ Added progress package dependency.
+ Internal changes to support Rcpp 0.12.11.

Gene Usage:

+ Renamed the count/frequency columns output by `countGenes` when the `clone` argument
  is specified to `CLONE_COUNT`/`CLONE_FREQ`.
+ Added a vignette describing basic gene usage analysis.


Version 0.2.6:  March 21, 2017
-------------------------------------------------------------------------------

General:

+ License changed to Creative Commons Attribution-ShareAlike 4.0 International
  (CC BY-SA 4.0).
+ Removed data.table dependency and added readr dependency.
+ Performance improvements in `readChangeoDb` and `writeChangeoDb`.


Version 0.2.5:  August 5, 2016
-------------------------------------------------------------------------------

General:

+ Fixed a bug in `seqDist()` wherein distance was not properly calculated in
  some sequences containing gap characters.
+ Added stop and gap characters to `getAAMatrix()` return matrix.


Version 0.2.4:  July 20, 2016
-------------------------------------------------------------------------------

General:

+ Added Rcpp and data.table dependencies.
+ Modified `readChangeoDb()` to wrap `data.table::fread()` instead of 
  `utils::read.table()` if the input file is not compressed.
+ Ported `testSeqEqual()`, `getSeqDistance()` and `getSeqMatrix()` to C++ to 
  improve performance of `collapseDuplicates()` and other dependent functions.
+ Renamed `testSeqEqual()`, `getSeqDistance()` and `getSeqMatrix()` to 
  `seqEqual()`, `seqDist()` and `pairwiseDist()`, respectively.
+ Added `pairwiseEqual()` which creates a logical sequence distance matrix;
  TRUE if sequences are identical, FALSE if not, excluding Ns and gaps.
+ Added translation of ambiguous and gap characters to `X` in 
  `translateDNA()`.
+ Fixed bug in `collapseDuplicates()` wherein the input data type sanity check
  would cause the vignette to fail to build under R 3.3.
+ Replaced the `ExampleDb.gz` file with a larger, more clonal, `ExampleDb` 
  data object.
+ Replaced `ExampleTrees` with a larger set of trees.
+ Renamed `multiggplot()` to `gridPlot()`.

Amino Acid Analysis:

+ Set default to `normalize=FALSE` for charge calculations to be more consistent
  with previously published repertoire sequencing results.
  
Diversity Analysis:

+ Added a `progress` argument to `rarefyDiversity()` and `testDiversity()` to
  enable the (previously default) progress bar.
+ Fixed a bug in `estimateAbundance()` were the function would fail if there 
  was only a single input sequence per group.
+ Changed column names in `data` and `summary` slots of `DiversityTest` to 
  uppercase for consistency with other tools.
+ Added dispatching of `plot` to `plotDiversityCurve` for `DiversityCurve`
  objects.
  
Gene Usage:

+ Added `sortGenes()` function to sort V(D)J genes by name or locus position.
+ Added `clone` argument to `countGenes()` to allow restriction of gene 
  abundance to one gene per clone.

Topology Analysis:

+ Added a set of functions for lineage tree topology analysis.
+ Added a vignette showing basic tree topology analysis.


Version 0.2.3:  February 22, 2016
-------------------------------------------------------------------------------

General:

+ Fixed a bug wherein the package would not build on R < 3.2.0 due to changes
  in `base::nchar()`.
+ Changed R dependency to R >= 3.1.2.


Version 0.2.2:  January 29, 2016
-------------------------------------------------------------------------------

General:

+ Updated license from CC BY-NC-SA 3.0 to CC BY-NC-SA 4.0.
+ Internal changes to conform to CRAN policies.

Amino Acid Analysis:

+ Fixed bug where arguments for the `aliphatic()` function were not being
  passed through the ellipsis argument of `aminoAcidProperties()`.
+ Improved amino acid analysis vignette.
+ Added check for correctness of amino acids sequences to `aminoAcidProperties()`.
+ Renamed `AA_TRANS` to `ABBREV_AA`.

Diversity:

+ Added evenness and bootstrap standard deviation to `rarefyDiversity()` 
  output.

Lineage:

+ Added `ExampleTrees` data with example output from `buildPhylipLineage()`.


Version 0.2.1:  December 18, 2015
-------------------------------------------------------------------------------

General:

+ Removed plyr dependency.
+ Added dplyr, lazyeval and stringi dependencies.
+ Added strict requirement for igraph version >= 1.0.0.
+ Renamed `getDNADistMatrix()` and `getAADistMatrix()` to `getDNAMatrix` and 
  `getAAMatrix()`, respectively.
+ Added `getSeqMatrix()` which calculates a pairwise distance matrix for a set 
  of sequences.
+ Modified default plot sizing to be more appropriate for export to PDF 
  figures with 7-8 inch width.
+ Added `multiggplot()` function for performing multiple panel plots.

Amino Acid Analysis:

+ Migrated amino acid property analysis from Change-O CTL to alakazam. 
  Includes the new functions `gravy()`, `bulk()`, `aliphatic()`, `polar()`, 
  `charge()`, `countPatterns()` and `aminoAcidProperties()`.

Annotation:

+ Added support for unusual TCR gene names, such as 'TRGVA*01'.
+ Added removal of 'D' label (gene duplication) from gene names when parsed 
  with `getSegment()`, `getAllele()`, `getGene()` and `getFamily()`.  May be 
  disabled by providing the argument `strip_d=FALSE`.
+ Added `countGenes()` to tabulate V(D)J allele, gene and family usage.

Diversity:

+ Added several functions related to analysis of clone size distributions, 
  including `countClones()`, `estimateAbundance()` and `plotAbundance()`.
+ Renamed `resampleDiversity()` to `rarefyDiversity()` and changed many of
  the internals. Bootstrapping is now performed on an inferred complete
  relative abundance distribution.
+ Added support for inclusion of copy number in clone size determination
  within `rarefyDiversity()` and `testDiversity()`.
+ Diversity scores and confiderence intervals within `rarefyDiversity()`
  and `testDiversity()` are now calculated using the mean and standard 
  deviation of the bootstrap realizations, rather than the median and
  upper/lower quantiles.
+ Added ability to add counts to the legend in `plotDiversityCurve()`.


Version 0.2.0:  June 15, 2015
-------------------------------------------------------------------------------

Initial public release.

General:

+ Added citations for the `citation("alakazam")` command.


Version 0.2.0.beta-2015-05-30:  May 30, 2015
-------------------------------------------------------------------------------

Lineage:

+ Added more error checking to `buildPhylipLineage()`.


Version 0.2.0.beta-2015-05-26:  May 26, 2015
-------------------------------------------------------------------------------

Lineage:

+ Fixed issue where `buildPhylipLineage()` would hang on R 3.2 due to R change 
  request PR#15508.


Version 0.2.0.beta-2015-05-05:  May 05, 2015
-------------------------------------------------------------------------------

Prerelease for review.
