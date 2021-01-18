BinaryDosage: Creates, Merges, and Reads Binary Dosage Files
================

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/USCbiostats/BinaryDosage?branch=master&svg=true)](https://ci.appveyor.com/project/USCbiostats/BinaryDosage)
[![Travis build
status](https://travis-ci.org/USCbiostats/BinaryDosage.svg?branch=master)](https://travis-ci.org/USCbiostats/BinaryDosage)
[![Codecov test
coverage](https://codecov.io/gh/USCbiostats/BinaryDosage/branch/master/graph/badge.svg)](https://codecov.io/gh/USCbiostats/BinaryDosage?branch=master)
<!-- badges: end -->

# Binary Dosage Files

### Introduction

Genotype imputation is an essential tool in genomics, enabling
association testing with markers not directly genotyped, increasing
statistical power, and facilitating data pooling between studies that
employ different genotyping platforms. Two commonly used software
packages for imputation are
[minimac](https://genome.sph.umich.edu/wiki/Minimac) and
[Impute2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html).
Furthermore, services such as the [Michigan Imputation
Server](https://imputationserver.sph.umich.edu/index.html) have made
genotype imputation much more accessible and streamlined.

While a number of software options are available for analyses of imputed
data (e.g. PLINK, EPACTS), fewer are available for Genomewide Gene x
Environment Interaction Scan (GWIS). Furthermore, data management tasks
such as parsing, subsetting, and merging, while manageable in smaller
studies, quickly become unwieldy and prohibitively slow with very large
samples sizes. We aim to address these limitations by converting
imputation outputs into a binary dosage file. The benefits of a binary
format are two fold - decreased hard drive storage requirements
(compared to a VCF file), and speed of parsing/analyses. The
BinaryDosage package contains functions to convert VCF and Impute2
formatted files into binary dosage files, along with functions to merge
samples.

For GWAS/GWIS analysis of BinaryDosage files, please refer to the
[**GxEScanR**](https://github.com/USCbiostats/GxEScanR) package.

### Description

##### Binary dosage data sets contain the following information:

  - Sample information
      - Family ID
      - Subject ID
  - SNP information
      - Chromosome number  
      - SNP ID  
      - Location in base pairs  
      - Reference allele  
      - Alternate allele  
  - Genetic information
      - Dosage values
      - Genotype probabilities, Pr(*g=0*), Pr(*g=1*), Pr(*g=2*)

There are 4 formats for a binary dosage data set. Data sets in formats
1, 2, and 3 have 3 files, a sample information file, a SNP information
file, and a genetic information file. Data sets in format 4 have just 1
file. This file contains all the information listed above and may
contain the following information.

**Note:** Format 4 is recommended and is the default value for all
functions.

  - Additional SNP information
      - Alternate allele frequency
      - Minor allele frequency
      - Average call rate
      - Imputation r squared
  - Merging information
      - Number of data sets merged
      - Sample size of each data set merged

### Functions

  - **vcftobd** - Converts a VCF file to a binary dosage data set
  - **gentobd** - Converts a GEN (impute2) file to a binary dosage data
    set
  - **bdmerge** - Merges multiple binary dosage data sets into a single
    data set
  - **getbdinfo** - Creates an R List containing information about a
    binary dosage data set (required for **getsnp** and **bdapply**)
  - **getvcfinfo** - Creates an R List containing information about a
    VCF file (required for **vcfapply**)
  - **getgeninfo** - Creates an R List containing information about a
    GEN file (required for **genapply**)
  - **bdapply** - Applies a function to the data for each SNP in a
    binary dosage file (requires list returned by **getbdinfo**)
  - **vcfapply** - Applies a function to the data for each SNP in a VCF
    file (requires list returned by **getvcfinfo**)
  - **genapply** - Applies a function to the data for each SNP in a GEN
    file (requires list returned by **getgeninfo**)
  - **getsnp** - Obtain genotype Dosages/Genotype Probabilities from a
    binary dosage file, outputs results to an R list

# Installation

1.  Install the [devtools](https://github.com/hadley/devtools) package
2.  Install the
    [BinaryDosage](https://github.com/USCbiostats/BinaryDosage) package
    directly from the USCbiostats repository on GitHub:

<!-- end list -->

``` r
remove.packages("BinaryDosage")
devtools::install_github("https://github.com/USCbiostats/BinaryDosage")

library(BinaryDosage)
```

# Usage

#### General Workflow

The general workflow for using binary dosage data sets is as follows:

  - Convert VCF or GEN files to a binary dosage data set
      - Note: When converting a VCF file to a binary dosage data set,
        the information file associated with the vcf can be used to add
        additional imputation information to the binary dosage data set
      - Note: When converting a GEN file to a binary dosage data set,
        the subject IDs can either be on the first line of the GEN file
        or in a separate sample file
  - Merge binary dosage datasets into a single data set
  - Apply a function to each SNP in the data set using bdapply
  - Extract SNPs for further analysis

#### Examples

The examples below use the default values for the functions. More
information about the functions and their options can be found using the
help files and the vignettes.

##### Example files

In the examples below the input files are included with the binary
dosage package and the output files are written to R temporary files. In
normal use, the user would provide the names of the input and output
files.

###### Input files

Example datasets *set1a.vcf* and *set1b.vcf* are representative of VCF
output files obtained from the Michigan Imputation Server. An
information file is also included for each set, *set1a.info* and
*set1b.info*.

Example datasets *set3a.imp* and *set3b.imp* are representative of files
return by the Impute imputation software. For GEN files the subject IDs
are contained in separated files. For this example these are
*set3a.sample* and *set3b.sample*.

The VCF and GEN files contain the same data. These files are in the
extdata directory of the BinaryDosage package. These sets contain the
following:

| Set   | Number of subjects | Number of SNPS |
| ----- | ------------------ | -------------- |
| 1a,3a | 60                 | 10             |
| 1b,3b | 40                 | 10             |

Since these files are distributed with the Binary Dosage package, it is
necessary to get the complete file name and path for use in the
following examples. The following code gets all the file names needed
for the examples.

``` r
library(BinaryDosage)

# Get the file names for the VCF and information files
vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
vcf1ainfo <- system.file("extdata", "set1a.info", package = "BinaryDosage")
vcf1bfile <- system.file("extdata", "set1b.vcf", package = "BinaryDosage")
vcf1binfo <- system.file("extdata", "set1b.info", package = "BinaryDosage")

# Get the file names for the GEN and sample files
gen3afile <- system.file("extdata", "set3a.imp", package = "BinaryDosage")
gen3asample <- system.file("extdata", "set3a.sample", package = "BinaryDosage")
gen3bfile <- system.file("extdata", "set3b.imp", package = "BinaryDosage")
gen3bsample <- system.file("extdata", "set3b.sample", package = "BinaryDosage")
```

###### Output files

The binary dosage output files will be written to temporary files. There
needs to be only one output file per data set because the examples use
the default format value of 4. The following code creates these
temporary output files.

``` r

# The output files for set 1
bdfile1a <- tempfile()
bdfile1b <- tempfile()
mergebd1 <- tempfile()

# The output files for set 3
bdfile3a <- tempfile()
bdfile3b <- tempfile()
mergebd3 <- tempfile()
```

##### Converting VCF files to a binary dosage data set

Converting a VCF file into a binary dosage file is simple. The user
passes the names of the VCF and information files along with the name
for the binary dosage file to the
<span style="font-family:Courier">vcftobd</span> function. There are
some options available for the
<span style="font-family:Courier">vcftobd</span> functions such as using
gz compressed files vcf files. More information about these options can
be found using the help files or reading the vignette
<span style="font-family:Courier">usingvcffiles</span>.

The following commands convert VCF data sets 1a and 1b into the binary
dosage format.

``` r

vcftobd(vcffiles = c(vcf1afile, vcf1ainfo), bdfiles = bdfile1a)
vcftobd(vcffiles = c(vcf1bfile, vcf1binfo), bdfiles = bdfile1b)
```

##### Converting GEN files to the Binary Dosage Format and Merging into one data set

Converting GEN files to binary dosage files is a little more difficult
than converting VCF files. This is because GEN files aren’t as strictly
formatted as VCF files. The user needs to have knowledge of how the GEN
file is formatted. More information on this can be found in the help
files and the vignette
<span style="font-family:Courier">usinggenfiles</span>.

In the example GEN file, the first column contains “--” for each SNP and
the second column contains the SNP ID in the format

<span style="font-family:Courier">\<chromosome\>:\<location\>\_\<reference
allele\>\_\<alternate allele\></span>

Because of this formatting, the function
<span style="font-family:Courier">gentobd</span> requires the
<span style="font-family:Courier">snpcolumns</span> parameter to have
the value <span style="font-family:Courier">c(0L, 2L:5L)</span>. To
convert the GEN data sets to binary dosage data sets, the names of the
input and output files are passed to
<span style="font-family:Courier">gentobd</span> along with the needed
value for snpcolumns.

The following commands convert the two GEN files into binary dosage
files.

``` r

gentobd(genfiles = c(gen3afile, gen3asample), snpcolumns = c(0L, 2L:5L), bdfiles = bdfile3a)
gentobd(genfiles = c(gen3bfile, gen3bsample), snpcolumns = c(0L, 2L:5L), bdfiles = bdfile3b)
```

##### Merging binary dosage files

Merging binary dosage files is done by SNP ID. The files to merge cannot
have the same subject IDs. See the vignette
<span style="font-family:Courier">usingbdfiles</span> for more
information. In this example we are assuming two separate groups of
subjects were imputed separately to the same reference panel.

To merge files, the user calls the
<span style="font-family:Courier">bdmerge</span> function and passes the
names of the files to merge along with a file name for the merged data
set. Other options exist for bdmerge and can be found in the help files
and the vignette <span style="font-family:Courier">mergingfiles</span>.

The following code first merges the binary files bdfile1a and bdfile1b
created from the VCF files into a single file, mergedbd1, and then does
the analogous action for the binary dosage files created from the GEN
files.

``` r
bdmerge(mergefiles = mergebd1, bdfiles = c(bdfile1a, bdfile1b))
bdmerge(mergefiles = mergebd3, bdfiles = c(bdfile3a, bdfile3b))
```

##### Applying a function to all the SNPs in a data set

Once binary dosage files have been created, a function can be applied to
all the SNPs in a file.

###### Defining the function

The function applied to the SNPs in a binary dosage file must have the
following four parameters, dosage, p0, p1, and p2. These are the dosage,
Pr(*g=0*), Pr(*g=1*), and Pr(*g=2*), respectively. Other parameters can
also be passed. For more information on defining the function see the
vignette <span style="font-family:Courier">usingbdfiles</span>.

The following code defines a function to calculate the alternate allele
frequency.

``` r
calculateaaf <- function(dosage, p0, p1, p2) {
    return(mean(dosage, na.rm = TRUE)/2)
}
```

###### Applying the function

To apply the function the user needs to call bdapply and pass
information about the binary dosage file and the function. The
information about the dosage file is obtained by calling the function
getbdinfo. If the user is going to call bdapply multiple times, the user
may wish to save the results of getbdinfo.

``` r
mergebd1info <- getbdinfo(mergebd1)
aaf1 <- bdapply(mergebd1info, calculateaaf)

mergebd3info <- getbdinfo(mergebd3)
aaf3 <- bdapply(mergebd3info, calculateaaf)
```

###### Checking the results

Since the VCF and GEN files contain the same information, the alternate
allele frequencies should be the same. The following code creates a data
frame with the SNP IDs and the alternate allele frequencies for both
data sets.

``` r
aaf <- cbind(mergebd1info$snps, aaf_set1 = unlist(aaf1), aaf_set3 = unlist(aaf3))
```

Here is a table showing the results.

| chromosome | location | snpid       | reference | alternate | aaf\_set1 | aaf\_set3 |
| :--------- | -------: | :---------- | :-------- | :-------- | --------: | --------: |
| 1          |    10000 | 1:10000:C:A | C         | A         |    0.3527 |    0.3527 |
| 1          |    11000 | 1:11000:T:C | T         | C         |    0.0135 |    0.0135 |
| 1          |    12000 | 1:12000:T:C | T         | C         |    0.2400 |    0.2400 |
| 1          |    13000 | 1:13000:T:C | T         | C         |    0.3375 |    0.3375 |
| 1          |    14000 | 1:14000:G:C | G         | C         |    0.1901 |    0.1901 |
| 1          |    15000 | 1:15000:A:C | A         | C         |    0.5627 |    0.5627 |
| 1          |    16000 | 1:16000:G:A | G         | A         |    0.4569 |    0.4569 |
| 1          |    17000 | 1:17000:C:A | C         | A         |    0.4578 |    0.4578 |
| 1          |    18000 | 1:18000:C:G | C         | G         |    0.2591 |    0.2591 |
| 1          |    19000 | 1:19000:T:G | T         | G         |    0.2431 |    0.2431 |

##### Extracting a SNP from the data set

After doing an analysis, the user may want to extract a SNP from the
data set for further analysis. This can be done using the getsnp
function. By default the function returns a list with the dosage values
for all the subjects. The genotype probabilities can be added to the
list by setting the dosageonly option to FALSE. See the help files or
the vignette <span style="font-family:Courier">usingbdfiles</span> for
more information.

The following code extracts the 6th SNP from both the binary dosage data
sets generated above.

``` r
# Get the dosage values for the 6th SNP
set1snp6 <- getsnp(mergebd1info, 6)
# Get the dosage values for the 6th SNP
set3snp6 <- getsnp(mergebd3info, 6)
```

The results from the above lines were merged into a data frame with the
subject IDs. Here are the first 10 lines of the data frame.

| subjectid | set1snp6 | set3snp6 |
| :-------- | -------: | -------: |
| I1        |    1.000 |    1.000 |
| I2        |    1.849 |    1.849 |
| I3        |    1.000 |    1.000 |
| I4        |    2.000 |    2.000 |
| I5        |    1.046 |    1.046 |
| I6        |    1.915 |    1.915 |
| I7        |    2.000 |    2.000 |
| I8        |    2.000 |    2.000 |
| I9        |    1.000 |    1.000 |
| I10       |    1.000 |    1.000 |
