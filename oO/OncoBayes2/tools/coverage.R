# Helper script to assess the test coverage of your package
# Change into the root directory of your package to run this script
#

library(OncoBayes2)
library(covr)

pkg_cov <- package_coverage()
print(pkg_cov)

# Run shiny application to explore in more detail which code parts are not yet
# testd well:
shine(pkg_cov)
