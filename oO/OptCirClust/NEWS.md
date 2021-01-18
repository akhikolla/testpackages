# NEWS

## Version 0.0.3

    2020-12-16
  1. Fixed inconsistent statements of parameters between the text
     and the code in the two tutorial vignettes.
  2. Included a missing header "cmath" that caused the package 0.0.2 fail 
     to build for Flavor r-oldrel-macos-x86_64


    2020-12-15
  1. Updated DESCRIPTION.
  2. Replaced a hardcoded infinity by numeric_limits<double>::infinity() in
     OptFramedClust.cpp.
  
  
    2020-12-13
  1. Version 0.0.3 was created from version 0.0.2


## Version 0.0.2

    2020-12-12
  
  1. Fixed an issue of missing the C++ header file <limits> in the OptFramedClust.h.
  2. Removed unused C++ header files.
  3. Add three more parameters (fill, border, border.lty) to give the user more control over visualizing circular data clustering results by function plot.CirClust().
  4. Updated vignettes and manuals.
  
  
    2020-12-09
  
  1. Version 0.0.2 was created from version 0.0.1
  2. Return value of CirClust are now fully documented in the manual according to the code.
  3. Checked the return value of FramedClust. 
  

## Version 0.0.1
    
    2020-12-07
    
  1. Updated plot.CirClust() function to avoid calling `par()` function. 
     This will give the user full control of graphic parameters.
     Updated the function's manual, examples, and related vignettes. 
  2. Updated package DESCRIPTION and README.md.
  3. Changed package title.
  
    
    2020-12-04
    
  1. Included with the package a FASTA file downloaded from
     GenBank accession number CP019943.1, representing the
     circular genome of Candidatus Carsonella ruddii, used
     in the vignette "Circular genome clustering".
  
  
    2020-12-03
    
  1. Fixed issues with example(FramedClust).
  2. Optimized plot.CirClust output  
    
    
    2020-12-01
  1. Fixed an issue where FramedClust crashed if the first.frame is
     not 1. E.g., 
        FramedClust(1:10, 2, 4, 2, 7)
        FramedClust(1:10, 2, 4, 2, 6)

    
    2020-11-30
    
  1. The orange arrows in the "Circular genome clustering" vignette
     are visible now.
  2. Changed the frame.width argument in FramedClust function to
     frame.size, to reduce potential confusion.
  
    
    2020-11-29
    
  1. In FramedClust, replaced the default 
        last.frame = length(X)-frame.width
     by 
        last.frame = length(X)-frame.width+1,
  
  2. Replaced the circular data examples in plot.FramedClust
  3. Defined the meaning of frame.width in the user document and the
  details section.
  4. Updated the FramedClust() manual. 
     + explained the NA in the 'cluster' component
     + removed the 'ID' component. 
     + removed "best.frame -> ID", and "Clearly define Border".
      
        
    2020-11-24

  1. Updated the package title and description. 
  
  
    2020-10-06
  
   1. Created a unified function for framed clustering 
  FramedClust(
    X, K, frame.width,     
    first.frame = 1,
    last.frame = length(X)-frame.width+1, 
    method = c("linear.polylog", "kmeans", "Ckmeans.1d.dp")
  )
  
<!--  2. create an R interface to lin_polylog_framed_clust 
  
  (do we need this seperate function, we are handeling all the framed clusters in the FramedClust function)
  
  lin.polylog.framed.clust(X, K, frame.width, first.frame, last.frame), inside which you call lin_polylog_framed_clust()
  
  For lin.polylog.framed.clust, the indices must be 1-based to be consistent with R;
  For lin_polylog_framed_clust, the indices must be 0-based to be consistent with C/C++.
-->

<!--  
    2020-09-29
  
  1. Do not export lin_polylog_framed_clust(), 
     kmeans.framed.clust(), quad.framed.clust() functions. Thus, we do not have to maintain three manuals that look almost identical.
     
  2. Reorder/rename the arguments to
  
  + lin_polylog_framed_clust(X, K, frame_width, first_frame, last_frame, prev, next)
  + kmeans.framed.clust(X, K, frame.width, first.frame, last.frame)
  + quad.framed.clust(X, K, frame.width, first.frame, last.frame) 
-->


    2020-09-13
  
  1. Completed the manual for function lin_polylog_framed_clust.
  2. Wrote the description of the package.
  3. Created a README.md file.
  4. Tested the kmeans.framed.clust() function.
    
    
    2020-08-11
    
  1. Package name changed to OptCirClust to avoid confusion with a
     similar R package.
  
    2020-08-10
  
  1. Added more tests for function CirClust(), such as irregular input data.
  2. Changed the for-loops in CirClust() to vectorized operations.
  3. Corrected the mean in OCC.
  
  
    2020-08-08
  
  1. Created CirClust version 0.0.1 from optCclust version 0.1.0
  2. Renamed function Cir.Clust() to CirClust() to reduce confusion.
     This function returns an object of class "CirClust".
  3. Renamed C/C++ function FOCC_Interface() to 
     lin_polylog_framed_clust().
  4. Renamed file "FOCC_Interface.cpp" to 
     "lin_polylog_framed_clust.cpp"
  5. Exported lin_polylog_framed_clust as a function for the package.

<!--  
    2020-08-08
  
  1. Shall we remove the following two arguments  (Yes)
       prev_k_f = -1,
       next_k_f = -1
       from FramedClust? I would think the user will never
       use them. Is so, will they cause confusion to the user?
-->  
  
