Version 1.1.0: August 10, 2020
-------------------------------------------------------------------------------

+ Fixed a bug in the clonal clustering methods causing TCR data to fail.
+ Added support for single-cell data in the clonal clustering methods, which
  are enabled by defining the optional `cell_id` column.


Version 1.0.1: May 24, 2020
-------------------------------------------------------------------------------

+ Fixed a fatal error in `identicalClones`, `hierarchicalClones` and 
  `spectralClones` when specifying `nproc` > 1.
  
  
Version 1.0.0: May 15, 2020
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
  accordingly.
+ Split `defineClonesScoper` function to three functions: `identicalClones`, 
  `hierarchicalClones`, and `spectralClones`.
  
General:

+ License changed to AGPL-3.

Cloning:

+ Fixed a platform precision incompatibility bug which caused spectral cloning
  results to be non-reproducible across platforms.
+ Added largest distance-to-nearest filter to clustering process.


Version 0.2.0:  August 5, 2019
-------------------------------------------------------------------------------

Deprecated:

+ function `analyzeClones` is deprecated. The clonal analysis has been added 
  to the main function `defineClonesScoper` as an argument `analyze_clones`. 
+ the out put class is deprecated. Results would be reported as a list if 
  argument `analyze_clones` set to be true, otherwise a single dataframe is
  returned.
+ `plot_neighborhoods` from clonal analysis has been deprecated.
+ `neighborhoods` from clonal analysis has been deprecated.

General:

+ New models, `hierarchical` for hierarchical-clustering based, and `identical` 
  for clustering among identical junction sequences are added. 
+ New method for spectral-clustering based model has been added through the 
  `vj` in argument `method`.

Clonal analysis:

+ Switched the meaning of the "inter" and "intra" labels in 
  `calculateInterVsIntra` function. Now, "inter" is the label used to form 
  distances that mean between clones, and "intra" is the label used to form 
  distances that mean on the inside, within each clone.
+ Changed the `plotInterVsIntra` output from a density plot to a histogram.
+ Changed the way to calculate the effective threshold. Now, desntiy approach 
  is used.
    

Version 0.1.0:  October 4, 2018
-------------------------------------------------------------------------------

Initial public release.
