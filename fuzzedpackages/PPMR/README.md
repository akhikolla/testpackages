# PPMR

PPMR (Probalistic polygenic two sample mendelian randomization), is an efficient R package for two sample MR analysis, accounting for the correlatded instruments and the horizontal pleiotropy. PPMR can provide the estimate of causal effect, the estimates of horizontal pleitropy, and the two corresponding p values.


# Installation
It is easy to install the development version of PPMR package using the 'devtools' package. The typical install time on a "normal" desktop computer is less than one minute.

```
# install.packages("devtools")
library(devtools)
install_github("yuanzhongshang/PPMR")
```


# Usage
The main functions is *PMR_individual* for individual level data, and *PMR_summary* for summary data. You can find the instructions by '?PMR_individual' and '?PMR_summary'. 

library(PPMR)

?PMR_individual

?PMR_summary

# Example

One simple example to use the package can be found at https://github.com/yuanzhongshang/PPMR/tree/master/example

# Results reproduced 

All results from all methods used in the PPMR paper can be reproduced at https://github.com/yuanzhongshang/PPMRreproduce

# Development
This R package is developed by Zhongshang Yuan and Xiang Zhou.

