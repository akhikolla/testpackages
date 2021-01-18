# CHANGELOG

Version 0.1.2 implements new functions to adapt the package to the high-performance computer.

## NEW FUNCTIONS

+ Version 0.1.2 adds a new method to get\_profile. It is called get\_profile\_all and it reads the info of all pixels that matched a substance/cluster. It only works for 'clipper' objects (actually, x should be a matrix). It will replace get\_profile\_sinfo some day. Is fast, it does not plot, it does its job.
+ Version 0.1.2 adds to new functions to load and write *SAM* objects. They are called sam\_write and sam\_load.
+ Version 0.1.2 adds a matrix\_sam function to perform a spectral angle mapper between to matrices (x can even be a vector). It is connected with the C++ sam\_internal which now is exported. I wrote it only for debugging purposes.
+ Version 0.1.2 adds a mosak_sam function to cut off pixles that have a poor match.

## FUNCTIONS THAT CHANGED

+ mosaic\_chunk was rewritten to add three extra argument to the function. The function can now read a chunk without pointing a SpectralInfo object. This change has the following implications:
	+ Arguments should be called explicitly, as some of them can be missing.
	+ If info is omitted you should provide the fpa size, the wavenumbers and the file path. You can omit the path if the connection to the \*.dmd file includes its full path.

## OTHER CHANGES

+ I added a configure file. Previous versions install well in Linux and Windows machines. However, they only work in macOS whenever the user have GDAL installed using homebrew. Compilation failed for users that isntall KingChaos's GDAL. The configure file solves the problem. It comes from rgdal package.
