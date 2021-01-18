#' Population constructor.
#'
#' The constructor for the \code{Population} object.
#'
#' \code{Population()} creates a new \code{Population} object based
#' on arguments which optionally modify a previously defined
#' \code{Population} object. If no \code{Population} object is given,
#' the new \code{Population} is created using only the arguments
#' given.
#'
#' The arguments \code{vcf}, \code{map}, \code{genotypes},
#' \code{literal} and \code{alleleFrequencies} all work together in a
#' specific way.
#'
#' If a VCF file is supplied via the \code{vcf} argument, the
#' \code{map}, \code{genotypes} and \code{alleleFrequencies}
#' arguments are not needed, since a map and set of genotypes are
#' given within the VCF file. If the number of individuals' genotypes
#' given by the VCF file does not match the number of individuals
#' specified by the \code{popSize} argument, the supplied genotypes
#' within the VCF file are used to suggest allele frequencies only:
#' this behaviour can also be forced by setting the \code{literal}
#' argument to \code{FALSE}.
#'
#' The \code{genotypes} argument supplies genotypes directly. In this
#' case, the user should supply a phased, individual-major genotypes matrix:
#' one individual per row and two columns per single nucleotide
#' polymorphism (SNP). Odd columns are assumed to be the haplotypes
#' inherited from the sires, while even columns are assumed to be the
#' haplotypes inherited from the dams. As with genotypes supplied via
#' a VCF file, if the number of individuals' genotypes given by the
#' matrix does not match the number of individuals specified by the
#' \code{popSize} argument, the supplied genotypes are used to
#' suggest allele frequencies only; again, this behaviour can also be
#' forced by setting the \code{literal} argument to \code{FALSE}.
#'
#' When supplying genotypes either directly or via a VCF file, all
#' SNPs should be biallelic and phased, with no missing values.
#' Genotypes supplied directly should have variants coded as either
#' 0 or 1.
#'
#' Any map (either supplied directly or via a VCF file) will be sorted,
#' such that all SNPs along the first chromosome listed will appear at the
#' start of the map, sorted in terms of base-pair distance; the second
#' chromosome to appear will then be treated similarly, and so on. SNPs
#' will be referenced within the population in this order.
#'
#' The \code{alleleFrequencies} argument is used when genotypes are
#' not given directly. In this case, the \code{literal} argument has
#' no meaning.
#'
#' An example \code{map data.frame} has been included in the epinetr
#' package as \code{map100snp}. Note that all chromosomes must be
#' autosomal, whether given via the \code{map} parameter or via a VCF
#' file.
#'
#' When supplying an existing \code{Population} object, any additive
#' effects and epistatic network will be carried over from the
#' previous population unless new QTLs are supplied.
#'
#' Note that if \code{broadH2} is equal to \code{narrowh2}, no
#' epistatic effects will be present; if \code{narrowh2} is 0, no
#' additive effects will be present; if \code{broadH2} is 1, no
#' environmental effects will be present.
#'
#' @param pop an optional \code{Population} object to use for default
#' values
#' @param popSize an optional number of individuals in the new
#' \code{Population} to be created
#' @param vcf an optional VCF file which will provide the map and
#' genotypes
#' @param map an optional \code{data.frame} specifying the name,
#' chromosome and position (in base-pairs) of each SNP
#' @param QTL an optional argument giving either a single number
#' which specifies the number of SNPs to randomly select as QTLs, or a
#' vector of SNP IDs (from \code{map}) to select as QTLs
#' @param genotypes an optional matrix of genotypes to use for the
#' population; see below for details
#' @param literal an optional \code{logical} value specifying whether
#' to use the genotypes directly or to generate new genotypes based
#' on the allele frequencies of the given genotypes
#' @param alleleFrequencies an optional vector of allele frequencies for
#' generating genotypes
#' @param broadH2 initial broad-sense heritability within the new
#' \code{Population}
#' @param narrowh2 initial narrow-sense heritability within the new
#' \code{Population}
#' @param traitVar initial phenotypic variance within the new
#' \code{Population}
#'
#' @return The constructor creates a new object of class
#' \code{'Population'}.
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @export
#'
#' @examples
#' # Construct a new population of size 500, with random allele
#' # frequencies and 20 QTLs chosen at random, broad-sense
#' # heritability set to 0.9, narrow-sense heritability set to 0.75
#' # and overall trait variance set to 40.
#'
#' pop <- Population(
#'   popSize = 500, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100), broadH2 = 0.9,
#'   narrowh2 = 0.75, traitVar = 40
#' )
#'
#' # Construct a new population of size 500 using directly supplied
#' # genotypes and 20 QTLs chosen at random, broad-sense heritability
#' # set to 0.7, narrow-sense heritability set to 0.3 and overall
#' # trait variance set to 10.
#'
#' pop2 <- Population(
#'   map = map100snp, genotypes = geno100snp,
#'   literal = TRUE, QTL = 20,
#'   broadH2 = 0.7, narrowh2 = 0.3, traitVar = 10
#' )
#'
#' # Modify the previous population to have narrow-sense heritabilty
#' # set to 0.45 and overall trait variance set to 20.
#'
#' pop2 <- Population(pop2, narrowh2 = 0.45, traitVar = 20)
#' @seealso \code{\link{addEffects}}, \code{\link{attachEpiNet}},
#' \code{\link{print.Population}}
Population <- function(pop = NULL, popSize = NULL, vcf = NULL, map = NULL, QTL = NULL,
                       genotypes = NULL, literal = TRUE, alleleFrequencies = NULL, broadH2 = NULL,
                       narrowh2 = NULL, traitVar = NULL) {

  # Create new Population object
  pop2 <- list()
  class(pop2) <- "Population"

  if (!is.null(pop)) {
    testPop(pop)
  }

  # Get population size
  pop2 <- addPopSize(pop, pop2, popSize, genotypes, literal)

  # Read VCF file
  if (!is.null(vcf)) {
    vcflist <- readVCF(vcf)

    if (is.null(map)) {
      map <- vcflist$map
    }

    if (is.null(genotypes)) {
      genotypes <- vcflist$genotypes
    }
  }

  # Add map to population
  pop2 <- addMap(pop, pop2, map)

  # Add QTL
  pop2 <- addQTL(pop, pop2, QTL)

  # Add haplotypes to the population
  pop2 <- addHaplotypes(pop, pop2, genotypes, literal, alleleFrequencies)

  # Add variances
  pop2 <- addVariances(pop, pop2, traitVar, broadH2, narrowh2)

  # Add sexes to the population
  pop2 <- addSex(pop, pop2)

  # Add pedigree data frame and animal IDs
  pop2$ID <- 1:pop2$popSize
  zeroes <- rep(0, pop2$popSize)
  pop2$ped <- data.frame(
    ID = pop2$ID, Sire = zeroes, Dam = zeroes, Additive = zeroes,
    Epistatic = zeroes, Environmental = zeroes, Phenotype = zeroes
  )

  # Add additive effects if QTL have not changed and narrow-sense
  # heritability is greater than zero
  if (is.null(QTL) && !is.null(pop$additive) && pop2$h2 > 0) {
    invisible(utils::capture.output(pop2 <- addEffects(pop2, pop$additive)))
  }

  # Add epistatic network if QTL have not changed and narrow-sense
  # heritability is less than broad-sense heritability
  if (is.null(QTL) && !is.null(pop$epiNet) && pop2$h2 < pop2$H2) {
    pop2$epiNet <- pop$epiNet
    pop2 <- calcEpiScale(pop2)
  }

  # Update pedigree data frame if we have all necessary effects in place
  if ((!is.null(pop2$additive) || pop2$h2 == 0) && (!is.null(pop2$epiNet) ||
    pop2$h2 == pop2$H2) || pop2$h2 == 0 && pop2$H2 == 0) {
    pop2 <- updatePedigree(pop2)
  }

  if (is.null(pop2$additive) && pop2$h2 > 0) {
    message("Run addEffects() to attach additive effects to population.")
  } else if (is.null(pop2$epiNet) && pop2$H2 > pop2$h2) {
    message("Run attachEpiNet() to attach epistatic effects to population.")
  }

  return(pop2)
}
