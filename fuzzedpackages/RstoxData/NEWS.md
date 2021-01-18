# RstoxData v1.0.17 (2020-11-23)

* In `writeICESDatras()`: Change country code 'NOR' to 'NO'.
* Change the default NumberOfCores to 1L in all parallel-able functions.
* Refactor `WriteICESAcoustic()`, `WriteICESBiotic()`, `WriteICESDatras()` into
`prepareICESAcoustic()`, `prepareICESBiotic()`, `prepareICESDatras()`.
* Tests: do not attempt to copy file outside `tempdir()` in `test-readXmlFile.R` file.
* Github actions: Update `check-full.yaml` file.
* `prepareICESAcoustic()` to use `data.table` and not use `format()`.
* Fix parallel `lapplyOnCores()` behavior in Windows platform.
* Remove `StoxAcousticStartMiddleStopDateTime()` function.

# RstoxData v1.0.16 (2020-11-11)

* Prepare for CRAN submission: Remove all attempts to modify the global environment.
* Prepare for CRAN submission: Shorten the package title.
* Update variables in `processDataSchema.json`.
* Delete `zzz.R` and `pkgname.R` and moved the contents to `RstoxData-package.R`.

# RstoxData v1.0.15 (2020-11-05)

* Added a `NEWS.md` file to track changes to the package.
* Added a `cran-comments.md` file for CRAN submission.
* `StoxExport` functions are now writing output files in `tempdir()`.
* Various small fixes for adhering to CRAN policies.

# RstoxData v1.0.14 (2020-10-29)

* New minor release.

# RstoxData v1.0.14 (2020-10-01)

* Changed SampleKey to catchsampleid.

# RstoxData v1.0.12 (2020-09-30)

Changes:
* StoxExport: Adopt the latest unified DATRAS format.
* StoxExport: Fix function for getting ICES ship code.
* stoxBioticObject: Fix DateTime column from ICES Biotic has a different class compared to the one from NMD Biotic.
* Documentation fix.

# RstoxData v1.0.11 (2020-09-18)

* New minor release.

# RstoxData v1.0.10 (2020-09-08)

* New `getStoxKeys()`.

# RstoxData v1.0.9 (2020-09-07)

Added `mergeByStoxKeys()`.

# RstoxData v1.0.8 (2020-09-01)

Added support for undefined `signifDigits`.

# RstoxData v1.0.7 (2020-09-01)

Modified `setPrecisionLevelOneDT()` to keep at least 6 significant digits.

# RstoxData v1.0.6 (2020-08-28)

Removed unit in variable and parameter names.

# RstoxData v1.0.5 (2020-08-25)

Added IDs separated by '-'.

# RstoxData v1.0.4 (2020-08-21)

New minor release.

# RstoxData v1.0.3 (2020-08-15)

Fix small bug in `DefineDataTranslation()`.

# RstoxData v1.0.2 (2020-08-14)

Added `Cruise` in `StoxAcoustic`.

# RstoxData v1.0.1.9001 (2020-07-17)

Added `Cruise` in the Cruise table of `StoxBiotic`.

# RstoxData v1.0.1.9000 (2020-07-17)

Added Redefine, Translate and Convert.

# RstoxData v1.0.1 (2020-07-11)

Changed `AcousticCategory` to character.

# RstoxData v0.8.15 (2020-06-21)

Added warning for non-unique `LogKey`.

# RstoxData v0.8.14 (2020-06-16)

Added `AddStoxBioticVariables()` and ConversionType in ConvertStoxBiotic.

# RstoxData v0.8.13 (2020-06-11)

Add NMD Biotic format v3.1 support.

# RstoxData v0.8.12 (2020-06-10)

Renamed Cores to NumberOfCores and LengthInterval to LengthIntervalCentimeters.

# RstoxData v0.8.11.9000 (2020-06-08)

Filter: Filtering upwards now able to skip empty tables.

# RstoxData v0.8.11 (2020-06-06)

Use Github actions.

# RstoxData v0.8.10.9002 (2020-06-02)

Appveyor: fix logic.

# RstoxData v0.8.10.9001 (2020-06-02)

Appveyor: test logic.

# stoxData v0.8.10 (2020-06-02)

New minor release.

# RstoxData v0.8.9 (2020-05-28)

Minor update to stox acoustic.

# RstoxData v0.8.8 (2020-05-19)

Minor update to stox acoustic.

# RstoxData v0.8.7 (2020-05-08)

Minor update to stox acoustic.

# RstoxData v0.8.6.9001 (2020-05-05)

Use strict lss format.

# RstoxData v0.8.6 (2020-05-04)

New features:
* Add multiple data export feature:
  1. ICES acoustic XML format to ICES acoustic CSV format
  2. NMD biotic v3 XML format to ICES biotic CSV format
  3. NMD biotic V3 XML format to NS-IBTS ICES Datras (CSV) format

Bug fixes:
* Fixes ICES Acoustic/Biotic vocabulary generation when parsing ICES's XML files

# RstoxData v0.8.5 (2020-05-01)

New features:
* Add filtering with upward propagation 

Bugfixes:
* Propagate down now able to skip empty tables (e.g., NMD v3 prey)

# RstoxData v0.8.4 (2020-04-28)

New minor release.

# RstoxData v0.8.3 (2020-04-16)

Add define and update variables.

# RstoxData v0.8.2 (2020-04-07)

More fixes in stox acoustic.

# RstoxData v0.8.1 (2020-04-02)

Fixes in stox acoustic.

# RstoxData v0.8.0.9002 (2020-04-01)

Minor bug fix.

# RstoxData v0.8.0 (2020-04-01)

New stox acoustic.

# RstoxData v0.7.0 (2020-03-27)

Fixed bug with non-unique tables from ICESAcoustic.

# RstoxData v0.7 (2020-03-26)

New stox acoustic.

# RstoxData v0.6.7 (2020-03-19)

New minor release. Notable new features are the support for NMD Biotic v1.4 and ICES Biotic data input.

Changelog:
* Updated README
* Increment version
* Add test for StoxBiotic
* Add ICES Biotic to StoXBiotic format conversion
* readXmlFile: Fix typo
* Add support for converting NMD Biotic v1.4 format to StoxBiotic

# RstoxData v0.6.6 (2020-02-11)

New minor release.

# RstoxData v0.6.5 (2019-12-19)

New minor release.

Changelog:
* Update README(.md)
* Increment version
* Refresh xsdObjects data
* Fix documentation
* Add verbose output switch for reading XMLs
* Add tests for StoxBiotic and StoxAcoustic functions
* Fix merging proses in StoxAcoustic() after correction of the column types
* Correct the data types in result tables
* Fix StoxBiotic function

# RstoxData v0.6.4 (2019-12-18)

Updated XML reader, RstoxBiotic and RstoxAcoustic functions.

Changes:
* Update DESCRIPTION
* Appveyor needs the whole R version digit
* Increment minor version
* Update CIs to use only R 3.5
* Merge branch 'biotic_acoustic' of github.com:StoXProject/RstoxData into biotic_acoustic
* Filter draft
* XML stream read speed improvements
* Result tables is now ordered
* Suppress warning on utils data
* metadata from readXmlFile is now a data.table
* Update ICES XSDs raw and compiled data objects
* Update stoxBioticObject data
* Fix processBioticData length resolution converter
* icesBiotic is now supported
* Fix NMDBiotic v.1x xsdObjects
* Added documentation for DataTypes
* Merge branch 'biotic_acoustic' of github.com:StoXProject/RstoxData into biotic_acoustic
* Added documentation for BioticData
* StoxBiotic: Remove duplicate rows from SpeciesCategory table
* Update data for conversion
* Latitude and Latitude2
* Latest StoxBiotic updates
* Merge branch 'biotic_acoustic' of github.com:StoXProject/RstoxData into biotic_acoustic
* Minor fix
* Fix (another) encoding problem
* Fix documents
* Update StoxBiotic conversion data and scripts
* Clean up data.table warnings and inefficiencies
* Removed PDFs
* Add StoxAcoustic and StoxBiotic

# RstoxData v0.6.3 (2019-12-06)

New minor release, various bug fixes.

# RstoxData v0.6.2 (2019-11-26)

New minor release, mainly for various bug fixes.

# RstoxData v0.6.1 (2019-06-21)

New Feature: Reading ICES acoustic XML files.

# RstoxData v0.6.0 (2019-06-12)

Changes:
* Rebrand to RstoxData.
* Add support for all biotic formats.

# RNMDAPI v0.5.1 (2019-05-27)

Changes:
* BUG FIX: Correctly handle decimal type in the resulting tables.

# RNMDAPI v0.5.0 (2019-04-09)

Changes:
* Supports reading ZIP compressed XML files in the push XML parser. This function is built on top of MINIZ compression library.

# RNMDAPI v0.4.0 (2019-04-06)

Changes:
* Faster and more memory efficient in reading XML files.
* New XML engine, all XML processing are now done inside the C++ code.
* A new streaming pull parser is implemented to avoid reading the whole big XML files into memory.
* Now supports landings v2 XML format.
