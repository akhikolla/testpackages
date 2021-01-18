.onAttach <- function(libname, pkgname) {
  f <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"),
                c("Version", "Date"))
  packageStartupMessage("Successfully loaded stepR package version ", f[1,1],".\n",
                        "Several new functions are added in version 2.0-0.",
                        " Some older functions are deprecated (still working) ",
                        "and may be defunct in a later version.",
                        " Please read the documentation for more details.",
                        domain = NULL, appendLF = TRUE)
}
