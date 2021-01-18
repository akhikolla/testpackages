
.onAttach <- function(...) {
  if (!interactive()){return()}
  
  # Create a list of helpful tips
  pkg_hints = c(
    "To see the user guides use `browseVignettes('optiSel')`"
  )
  
  # Randomly pick one hint
  startup_hint = sample(pkg_hints, 1)
   
  # Display hint
  packageStartupMessage(paste(strwrap(startup_hint), collapse = "\n"))
}