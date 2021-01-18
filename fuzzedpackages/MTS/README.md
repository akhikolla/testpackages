MTS
===
Installation

```
git clone git://github.com/d-/MTS.git
R CMD build MTS/
```

This will create a file named "MTS_VERSION.tar.gz".

Then move the file into your working directory in R and type:
```
install.packages("MTS_VERSION.tar.gz",repos=NULL,type="source")
library(MTS)
```

Alternatively, a simpler solution is to use the 'devtools' package.

```
install.packages("devtools")
library(devtools)
install_github('MTS','d-')
```
