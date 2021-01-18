tabHome = fluidPage(

  column(
    box(
      div(
        h1(paste("Genepop-shiny: version", config$version)),
        p("This project combines the 'shiny' and 'genepop' R packages to bring a user interface to the Genepop software."),
        br(), br(),
        h1(paste("Genepop: version", config$version_genepop)),
        p("Genepop implements a mixture of traditional methods and some more focused developments:"),
        tags$ul(
          tags$li("It computes exact tests for Hardy-Weinberg equilibrium, for population differentiation and for genotypic disequilibrium among pairs of loci;"),
          tags$li("It computes estimates of F-statistics, null allele frequencies, allele size-based statistics for microsatellites, etc., and of number of immigrants by Barton and Slatkin’s 1986 private allele method;"),
          tags$li("It performs analyses of isolation by distance from pairwise comparisons of individuals or population samples, including confidence intervals for neighborhood size.")
        )
      ),
      width = NULL
    ),
    width = 12)
)
