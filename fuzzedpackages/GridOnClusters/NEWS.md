---
output:
  pdf_document: default
  html_document: default
---
# NEWS

## Version 0.0.8
  
  2020-09-13
  
  1. Added function plotGOCpatterns to plot the continuous data along with the
  cluster preserving grid.
  2. Created a manual for the plotGOCpatterns() function.
  3. Updated the code for the Examples vignette to use plotGOCpatterns.

  2020-08-10
  
  1. Created version 0.0.8 from 0.0.7.
  2. Added an additional parameter 'min_level' to denote the minimum number of 
  discretization levels required for each dimension. 
  3. Updated the manual of discrete.jointly() function.
  4. Added an entry in reference and citation.
  5. Updated README.md with badges.

## Version 0.0.7
  
  2020-04-03
  
  1. Tidied up code for the Examples vignette.
  2. Updated the manual of discrete.jointly() function.
  3. Made minor editorial changes in DESCRIPTION and README.md.
  4. Resolved signed/unsigned mismatches.

  2020-03-31
  
  1. Created version 0.0.7 from 0.0.6.
  2. Fixed memory leak in Clusters.cpp when calculating median.


## Version 0.0.6

  2020-03-26
  
  1. Created version 0.0.6 from 0.0.5.
  2. Rewrote Prep_Index() to work in between two consecutive points, 
  rather than on top of a single point.
  3. Using distance() in Prep_Index() to calculate the distance for two 
  iterators.
  4. Using "ceil" in Binary_Index_Searching() to consider even/odd cases
  when determining grid lines.
  5. Fixed potential overflow issues.

## Version 0.0.5 (not released to the public)

  2020-03-25
  
  1. Fixed a bug in the prep_index() function.
  2. Fixed prep_index() (lines 120 and 125) such that grid lines are
  put at the midpoint between two conseuctive points, instead of on
  one of the points.
  3. Updated vignette. Example 2 seems always correct now.
  
  2020-03-24

  1. Created version 0.0.5 from 0.0.4.
  2. Function discretize.jointly() now returns cluster labels of each
  observation and a similarity score (ARI) between the joint
  discretization and the cluster labels of each observation. 
  3. The class Cluster has a new constructor that takes cluster
  labels and the input data to compute median for each dimension.
  4. Find_grid() is now based on median.
  5. Using 'dqrng' in test cases to avoid RNG issue in testing.
  6. Rewrote multiple functions in Joint_Grid.cpp to avoid push_back().
  7. New visualization code in vignette now shows cluster labels for
  each observation.

## Version 0.0.4 (not released to the public)

  2020-03-20
  
  1. Created version 0.0.4 from 0.0.3.
  2. Fixed typos in DESCRIPTION and README files.

## Version 0.0.3

  2020-03-17
  
  1. Created version 0.0.3 from 0.0.2. Package renamed to GridOnClusters
  2. Function joint.grid.discretize.R() renamed to discretize.jointly()
  3. Return values of function discretize.jointly() changed to include
  both the discretized data and the grid
  4. Manual for discretize.jointly() updated.
  5. Line 104, 105 in Joint_Grid.cpp commented out
  6. Rewrote Find_Grid() to avoid push_back() in Joint_Grid.cpp
  7. Created a vignette to include examples.
  
## Version 0.0.2 (not released to the public)

  2020-03-14
  
  1. Created the initial version 0.0.1. Package renamed to QNJGD

## Version 0.0.1 (not released to the public)

  2020-03-09
  
  1. Created the initial version 0.0.1. Package named JointGridDiscr
