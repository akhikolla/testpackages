# RstoxData

![R-CMD-check](https://github.com/StoXProject/RstoxData/workflows/R-CMD-check/badge.svg)

## Introduction

RstoxData is a package to parse and manipulate various data formats for fisheries.

The package contains functions for reading, filtering and writing fisheries' trawl survey samples (biotic), acoustic trawl survey and commercial catch samples (landings) data in XML formats.

Loaded data can also be filtered by using any supported R conditional syntax such as `longitude > 10`, or by user pre-defined functions such as `inside()`.

### Supported formats:

#### A. Norwegian Institute of Marine Research (IMR) ([data definitions](https://www.imr.no/formats/)):

  1. Biotic (fish trawl survey) XML format (up to version 3.1).
  2. Acoustic (acoustic/echosounder trawl surveys) XML format (version 1).
  3. Landings (commercial catch samples) XML format (version 2).

#### B. Norwegian Directorate of Fisheries:
  1. Sales Notes data in the LSS format.
  2. Electronic logbooks (ERS) in tabular format.

#### C. International Council for the Exploration of the Sea (ICES) ([data definitions](https://ices.dk/data/data-portals/Pages/acoustic.aspx)):

  1. Biotic (fish trawl survey) XML format.
  2. Acoustic (acoustic/echosounder trawl surveys) XML format.

## Installation

1. Install from CRAN:

    ```r
    install.packages("RstoxData")
    ```

2. Install the latest release from our local repository:
    ```r
    install.packages("RstoxData", repos = c("https://stoxproject.github.io/repo/", "https://cloud.r-project.org/"))
    ```

3. Install the latest version from GitHub:
    ```r
    devtools::install_github("https://github.com/StoXProject/RstoxData")
    ```

On computers that return errors when trying to run the Rtools through RStudio (most institutional Windows machines), install the binary directly from https://github.com/StoXProject/RstoxData/releases.
Download the newest RstoxData zip file, click the `Packages` tab -> `Install` -> `Install from:` `Package Archive File` -> `Install`. If the installer does not complain, the package is installed correctly.

## License

LGPL-3 © Norwegian Institute of Marine research (IMR) ([homepage](https://www.hi.no/en)).

The development of RstoxData package is mainly supported by IMR's [REDUS](http://www.redus.no) and SEA2DATA projects.

RstoxData is part of the bigger [StoX ecosystem](https://stoxproject.github.io).

---

### For historical release notes, see: https://github.com/StoXProject/RstoxData/blob/master/NEWS.md