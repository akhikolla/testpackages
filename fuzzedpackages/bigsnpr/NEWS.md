## bigsnpr 1.5.1

- Add multiple checks in `snp_readBGEN()` to make sure of the expected format.

## bigsnpr 1.5.0

- Add function `snp_fst()` for computing Fst. 

## bigsnpr 1.4.11

- Workaround for error `could not find function "ldpred2_gibbs_auto"`.

## bigsnpr 1.4.9 & bigsparser 0.4.0

- Can now directly do `as_SFBM(corr0)` instead of `bigsparser::as_SFBM(as(corr0, "dgCMatrix"))`. This should also use less memory and be faster.

## bigsnpr 1.4.8

- Add option `sparse` to enable getting also a sparse solution in LDpred2-auto.

## bigsparser 0.3.0

- Faster `bigsparser::as_SFBM()`.

## bigsnpr 1.4.7

- Allow for format `01` or `1` for chromosomes in BGI files.

## bigsnpr 1.4.6

- Fasten `snp_match()`. Also now remove duplicates by default.

## bigsnpr 1.4.3

- Fix a bug when using very large correlation matrices in LDpred2 (although we do not recommend to do so).

## bigsnpr 1.4.2

- All 3 LDpred2 functions now use an SFBM as input format for the correlation matrix.

- Allow for multiple initial values for p in `snp_ldpred2_auto()`.

- Add function `coef_to_liab()` for e.g. converting heritability to the liability scale.

## bigsnpr 1.4.1

- Change default of parameter `alpha` of function `snp_cor()` to `1`.

## bigsnpr 1.4.0

- Add functions `snp_ldpred2_inf()`, `snp_ldpred2_grid()` and `snp_ldpred2_auto()` for running the new LDpred2-inf, LDpred2-grid and LDpred2-auto.

- Add functions `snp_ldsc()` and `snp_ldsc2()` for performing LD score regression.

- Add function `snp_asGeneticPos()` for transforming physical positions to genetic positions.

- Add function `snp_simuPheno()` for simulating phenotypes.

## bigsnpr 1.3.1

- Also use OpenMP for the parallelization of `snp_pcadapt()`, `bed_pcadapt()`, `snp_readBGEN()` and `snp_fastImputeSimple()`.

## bigsnpr 1.3.0

- Parallelization of clumping algorithms has been modified. Before, chromosomes were imputed in parallel. Now, chromosomes are processed sequentially, but computations within each chromosome are performed in parallel thanks to OpenMP. This should prevent major slowdowns for very large samples sizes (due to swapping).

- Use OpenMP to parallelize other functions as well (possibly only sequential until now).

## bigsnpr 1.2.6

- Can now run `snp_cor()` in parallel.

- Parallelization of `snp_fastImpute()` has been modified. Before this version, chromosomes were imputed in parallel. Now, chromosomes are processed sequentially, but computation of correlation between variants and XGBoost models are performed using parallelization.

## bigsnpr 1.2.5

- Add function `snp_subset()` as alias of method `subset()` for subsetting `bigSNP` objects.

## bigsnpr 1.2.4

- Use new class `bed_light` internally to make parallel algorithms faster because they have to transfer less data to clusters. Also define differently functions used in `big_parallelize()` for the same reason.

## bigsnpr 1.2.0

- Use the new implementation of robust OGK Mahalanobis distance in {bigutilsr}.

## bigsnpr 1.1.1

- Fix error `object 'obj.bed' not found` in `snp_readBed2()`.

## bigsnpr 1.1.0

- Cope with new read-only option in {bigstatsr} version >= 1.1.

## bigsnpr 1.0.2

- Add option `backingfile` to `subset.bigSNP()`.

## bigsnpr 1.0.1

- Add option `byrow` to `bed_counts()`.

## bigsnpr 1.0.0

- Add memory-mapping on PLINK (.bed) files with missing values + new functions:
    - `bed()`
    - `bed_MAF()`
    - `bed_autoSVD()`
    - `bed_clumping()`
    - `bed_counts()`
    - `bed_cprodVec()`
    - `bed_pcadapt()`
    - `bed_prodVec()`
    - `bed_projectPCA()`
    - `bed_projectSelfPCA()`
    - `bed_randomSVD()`
    - `bed_scaleBinom()`
    - `bed_tcrossprodSelf()`
    - `download_1000G()`
    - `snp_modifyBuild()`
    - `snp_plinkKINGQC()`
    - `snp_readBed2()`
    - `sub_bed()`
    
- Add 3 parameters to `autoSVD()`: `alpha.tukey`, `min.mac` and `max.iter`.

- Remove option for changing ploidy (that was only partially supported).

- Automatically apply `snp_gc()` to `pcadapt`.

## bigsnpr 0.12.0

- Add `snp_fastImputeSimple()`: fast imputation via mode, mean or sampling according to allele frequencies.

## bigsnpr 0.11.3

- Fix a bug in `snp_readBGEN()` that could not handle duplicated variants or individuals.

## bigsnpr 0.11.1

- When using `snp_grid_PRS()`, it now stores not only the FBM, but also the input parameters as attributes (the whole result basically).

## bigsnpr 0.11.0

- Add 3 SCT functions `snp_grid_*()` to improve from Clumping and Thresholding (preprint coming soon).

- Add `snp_match()` function to match between summary statistics and some SNP information.

## bigsnpr 0.10.2

- Parameter `is.size.in.bp` is deprecated.

## bigsnpr 0.10.1

- Add parameter `read_as` for `snp_readBGEN()`. It is now possible to sample BGEN probabilities as random hard calls using `read_as = "random"`. Default remains reading probabilities as dosages.

## bigsnpr 0.10.0

- For memory-mapping, now use *mio* instead of *boost*.

- `snp_clumping()` (and `snp_autoSVD()`) now has a `size` that is inversely proportional to `thr.r2`.

- `snp_pruning()` is deprecated (and will be removed someday); now always use `snp_clumping()`.

## bigsnpr 0.9.0

- When reading bed files, switch reading of Os and 2s to be consistent with other software.

## bigsnpr 0.8.2

- Add function `snp_assocBGEN()` for computing quick association tests from BGEN files. Could be useful for quick screening of useful SNPs to read in bigSNP format. This function might be improved in the future.

## bigsnpr 0.8.1

- Change url to download PLINK 1.9.

## bigsnpr 0.8.0

- Add function `snp_readBGEN()` to read UK Biobank BGEN files in `bigSNP` format.

## bigsnpr 0.7.0

- Add parameter `is.size.in.bp` to `snp_autoSVD()` for the clumping part.

- Change the threshold of outlier detection in `snp_autoSVD()` (it now detects less outliers). See the documentation details if you don't have any information about SNPs.

## bigsnpr 0.6

- Keep up with {bigstatsr}.

## bigsnpr 0.3.1

- Provide function `snp_gene` (as a gist) to get genes corresponding to 'rs' SNP IDs thanks to package {rsnps} from rOpenSci. See README.

## bigsnpr 0.3.0

- **Package {bigsnpr} is published in [Bioinformatics](http://dx.doi.org/10.1093/bioinformatics/bty185)**.

## bigsnpr 0.2.7

- Faster defaults + possibility to estimate correlations based on a subset of individuals for `snp_fastImpute`. Also store information in an FBM (instead of a data frame) so that imputation can be done by parts (you can stop the imputation by killing the R processes and come back to it later). Note that the defaults used in the *Bioinformatics* paper were `alpha = 0.02` and `size = 500` (instead of `1e-4` and `200` now, respectively). These new defaults are more stringent on the SNPs that are used, which makes the imputation faster (30 min instead of 42-48 min), without impacting accuracy (still 4.7-4.8% of errors).

## bigsnpr 0.2.5

- **This package won't be on CRAN**.

## bigsnpr 0.2.4

- No longer download PLINK automatically (because it is a CRAN policy violation).
