# Documentation
docOrg <- lapply(c(dir("man", pattern="*.md", recursive=TRUE, full.names = TRUE), "NAMESPACE"), readLines)
devtools::document()
docNew <- lapply(c(dir("man", pattern="*.md", recursive=TRUE, full.names = TRUE), "NAMESPACE"), readLines)
if (!identical(docOrg, docNew)) {
     stop("Documentation was not updated, try running devtools::document().")
}
# Pkgdown
# if (!require(pkgdown, quietly = TRUE)) {
#     install.packages("pkgdown")
# }
# 
# list_files = c(dir("docs", recursive=TRUE, full.names = TRUE), "NAMESPACE")
# list_files = list_files[!list_files=="docs/pkgdown.yml"]  # contains datetime
# list_files = list_files[!list_files=="docs/reference/index.html"]  # contains datetime
# list_files = list_files[!list_files=="docs/articles/glmdisc.html"]  # contains random stuff
# print(list_files)
# docOrg <- lapply(list_files, readLines)
# pkgdown::build_site()
# docNew <- lapply(list_files, readLines)
# print(length(docOrg))
# print(length(docNew))
# for (j in 1:length(docOrg)) {
#      if (!identical(docOrg[[j]], docNew[[j]])) {
#           print(paste0("File ", j, " is not OK."))
#      } else {
#           print(paste0("File ", j, " OK."))
#      }
# }
# if (!identical(docOrg, docNew)) {
#      stop("Pkgdown website was not updated, try running pkgdown::build_site().")
# }
