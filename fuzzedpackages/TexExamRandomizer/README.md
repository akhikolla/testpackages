# TexExamRandomizer
[![Build Status](https://travis-ci.org/alexrecuenco/TexExamRandomizer.svg?branch=master)](https://travis-ci.org/alexrecuenco/TexExamRandomizer)
R package, randomizing 'LaTeX' exams and grading them. 

Randomizes 'LaTeX' exams with a flexible set of options that can be provided directly in the document, using JSON-format. 

Look at the vignettes of the package to see all the options provided in the package.

Look as well inside the inst/extdata folder to see examples of different formats that the software understands

## Installation.

To install, 

    devtools::install_github("alexrecuenco/TexExamRandomizer")
    
## Using the package with TexShop (In MAC OS)


The package already comes with the intention of being used directly on a tipical TeX typewriter, such as TexShop. 

You will find inside the exec/ folder a couple of `.engine` files. 
By placing those files in your "engines" folder in your TexShop distribution, 
and by placing the corresponding executable scripts found on that same folder somewhere that can be found on your path, 
you will be able to use this package to create random exams without even having to move out of your TeX environment once it is set up.


