MenuGauche = sidebarMenu(id = "sidebarmenu",

    menuItem("Home", tabName = "Home",  icon = icon("home", lib="font-awesome")),

    menuItem("Load Data", tabName = "Load",  icon = icon("file-text-o", lib="font-awesome")),

    #menuItem("Documentation genepop", tabName = "Docu",  icon = icon("book", lib="font-awesome")),

    menuItem("Hardy-Weinberg exact tests", tabName = "NULL", icon = icon("pencil", lib = "font-awesome"),
       menuSubItem("H1= Heterozygote deficiency", tabName = "Opt11", href = NULL, newtab = TRUE,
         icon = shiny::icon("angle-double-right"), selected = NULL  ),
       menuSubItem("H1= Heterozygote excess", tabName = "Opt12", href = NULL, newtab = TRUE,
         icon = shiny::icon("angle-double-right"), selected = NULL  ),
       menuSubItem("Probability test ", tabName = "Opt13", href = NULL, newtab = TRUE,
         icon = shiny::icon("angle-double-right"), selected = NULL  ),
       menuSubItem("Global Heterozygote deficiency", tabName = "Opt14", href = NULL, newtab = TRUE,
         icon = shiny::icon("angle-double-right"), selected = NULL  ),
       menuSubItem("Global Heterozygote excess", tabName = "Opt15", href = NULL, newtab = TRUE,
           icon = shiny::icon("angle-double-right"), selected = NULL  )
    ),

    menuItem("haploid and genotypic disequilibrium", tabName = "NULL", icon = icon("pencil", lib = "font-awesome"),
      menuSubItem("Test each pair of loci in each population", tabName = "Opt21", href = NULL, newtab = TRUE,
         icon = shiny::icon("angle-double-right"), selected = NULL),

      menuSubItem("Create genotypic contingency tables", tabName = "Opt22", href = NULL, newtab = TRUE,
         icon = shiny::icon("angle-double-right"), selected = NULL)
    ),

  menuItem("Population differentiation", tabName = "NULL", icon = icon("pencil", lib = "font-awesome"),
     menuSubItem("Genic differentiation: all pop", tabName = "Opt31", href = NULL, newtab = TRUE,
       icon = shiny::icon("angle-double-right"), selected = NULL  ),

     menuSubItem("Genic differentiation: all pairs of population", tabName = "Opt32", href = NULL, newtab = TRUE,
       icon = shiny::icon("angle-double-right"), selected = NULL  ),

     menuSubItem("Genotypic differentiation: all pop", tabName = "Opt33", href = NULL, newtab = TRUE,
       icon = shiny::icon("angle-double-right"), selected = NULL  ),

     menuSubItem("Genotypic differentiation: all pairs of population", tabName = "Opt34", href = NULL, newtab = TRUE,
       icon = shiny::icon("angle-double-right"), selected = NULL  )
  ),

  menuItem("Nm estimates (private allele method)", tabName = "Opt41",  icon = icon("pencil", lib="font-awesome")),

  menuItem("Allele frequencies, various Fis and gene diversities", tabName = "NULL", icon = icon("pencil", lib = "font-awesome"),
           menuSubItem("Allele and genotype frequencies per locus and per sample", tabName = "Opt51", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Gene diversities and Fis using allele identity", tabName = "Opt52", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Gene diversities and Fis using allele size", tabName = "Opt53", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  )
  ),

  menuItem("Fst and other correlations, isolation by distance", tabName = "NULL", icon = icon("pencil", lib = "font-awesome"),
           menuSubItem("Allele identity (F-statistics) for all populations", tabName = "Opt61", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Allele identity (F-statistics) for all population pairs", tabName = "Opt62", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Allele size (Rho-statistics) for all populations", tabName = "Opt63", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Allele size (Rho-statistics) for all population pairs", tabName = "Opt64", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Isolation by distance between individuals", tabName = "Opt65", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Isolation by distance between groups", tabName = "Opt66", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  )
  ),

  menuItem("Ecumenicism: file conversion", tabName = "NULL", icon = icon("pencil", lib = "font-awesome"),
           menuSubItem("GENEPOP to FSTAT (F statistics)", tabName = "Opt71", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("GENEPOP to BIOSYS (letter code)", tabName = "Opt72", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("GENEPOP to BIOSYS (number code)", tabName = "Opt73", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("GENEPOP to LINKDOS (D statistics)", tabName = "Opt74", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  )
  ),

  menuItem("Null alleles and miscellaneous input file utilities", tabName = "NULL", icon = icon("pencil", lib = "font-awesome"),
           menuSubItem("Null allele: estimates of allele frequencies", tabName = "Opt81", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Diploidisation of haploid data", tabName = "Opt82", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Relabeling alleles", tabName = "Opt83", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Conversion to individual data with population names", tabName = "Opt84", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Conversion to individual data with individual names", tabName = "Opt85", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  ),
           menuSubItem("Random sampling of haploid genotypes from diploid ones", tabName = "Opt86", href = NULL, newtab = TRUE,
                       icon = shiny::icon("angle-double-right"), selected = NULL  )
  ),

  tags$br(), tags$br(), tags$br(),
  actionButton("close", "Close current session", icon("ban"), style=config$color$stopButton),
  #tags$br(), tags$br(), tags$br(),
  #actionButton("interrupt", "Interrupt genepop", icon("ban"), style=config$color$stopButton),
  tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),tags$br(),

  menuItem("Team", icon = icon("book", lib="font-awesome"),
        menuItem("François Rousset",  href = "http://www.isem.univ-montp2.fr/recherche/equipes/genetique-evolutive/personnel/rousset-francois/", newtab = TRUE,     icon = shiny::icon("male"), selected = NULL  ),
        menuItem("Khalid Belkhir",  href = "http://www.isem.univ-montp2.fr/recherche/les-plate-formes/bioinformatique-labex/personnel/belkhir-khalid/", newtab = TRUE,     icon = shiny::icon("male"), selected = NULL  ),
        menuItem("Jimmy Lopez",  href = "http://www.isem.univ-montp2.fr/recherche/les-plate-formes/bioinformatique-labex/personnel/", newtab = TRUE,   icon = shiny::icon("male"), selected = NULL  )
      )
  )
