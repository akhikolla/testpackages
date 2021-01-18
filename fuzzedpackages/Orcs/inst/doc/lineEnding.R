## ----pkg, echo=FALSE, message=FALSE-------------------------------------------
library(Orcs)

## ----lineEnding, eval=FALSE---------------------------------------------------
#  ## input file
#  infile <- file.path(system.file(package = "Orcs"), "DESCRIPTION")
#  
#  system(paste("file", infile))
#  # > C:/Users/.../R/win-library/3.5/Orcs/DESCRIPTION: ASCII English text, with CRLF line terminators
#  
#  ## convert to dos line endings and write to output file
#  outfile = file.path(tempdir(), "DESCRIPTION4wd")
#  lineEnding(infile, outfile = outfile, to = "unix")
#  
#  system(paste("file", outfile))
#  # > C:\Users\...\AppData\Local\Temp\RtmpMX3o1b/DESCRIPTION4wd: ASCII English text

