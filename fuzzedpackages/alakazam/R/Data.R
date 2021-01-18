# Documentation and definitions for data and constants

#### Sysdata ####

# 1x20 vector of default amino acid hydropathy scores
# HYDROPATHY_KYTJ82

# 1x20 vector of default amino acid bulkiness scores
# BULKINESS_ZIMJ68

# 1x20 vector of default amino acid polarity scores
# POLARITY_GRAR74

# 1x7 vector of default amino acid pK values
# PK_EMBOSS

#### Data ####

#' Example AIRR database
#'
#' A small example database subset from Laserson and Vigneault et al, 2014.
#'
#' @format   A data.frame with the following AIRR style columns:
#'   \itemize{
#'     \item  \code{sequence_id}:                Sequence identifier
#'     \item  \code{sequence_alignment}:         IMGT-gapped observed sequence.
#'     \item  \code{germline_alignment}:         IMGT-gapped germline sequence.
#'     \item  \code{germline_alignment_d_mask}:  IMGT-gapped germline sequence with N, P and 
#'                                               D regions masked.
#'     \item  \code{v_call}:                     V region allele assignments.
#'     \item  \code{v_call_genotyped}:           TIgGER corrected V region allele assignment.
#'     \item  \code{d_call}:                     D region allele assignments.
#'     \item  \code{j_call}:                     J region allele assignments.
#'     \item  \code{c_call}:                     Isotype (C region) assignment.
#'     \item  \code{junction}:                   Junction region sequence.
#'     \item  \code{junction_length}:            Length of the junction region in nucleotides.
#'     \item  \code{np1_length}:                 Combined length of the N and P regions proximal
#'                                               to the V region.
#'     \item  \code{np2_length}:                 Combined length of the N and P regions proximal
#'                                               to the J region.
#'     \item  \code{duplicate_count}:            Copy count (number of duplicates) of the sequence.
#'     \item  \code{clone_id}:                   Change-O assignment clonal group identifier.
#'     \item  \code{sample_id}:                  Sample identifier. Time in relation to vaccination.
#' }
#' 
#' @seealso \link{ExampleDbChangeo} \link{ExampleTrees}
#' 
#' @references
#' \enumerate{
#'   \item  Laserson U and Vigneault F, et al. High-resolution antibody dynamics of 
#'            vaccine-induced immune responses. 
#'            Proc Natl Acad Sci USA. 2014 111:4928-33.
#' }
"ExampleDb"

#' Example Change-O database
#'
#' A small example database subset from Laserson and Vigneault et al, 2014.
#'
#' @format   A data.frame with the following Change-O style columns:
#'   \itemize{
#'     \item  \code{SEQUENCE_ID}:           Sequence identifier
#'     \item  \code{SEQUENCE_IMGT}:         IMGT-gapped observed sequence.
#'     \item  \code{GERMLINE_IMGT_D_MASK}:  IMGT-gapped germline sequence with N, P and 
#'                                          D regions masked.
#'     \item  \code{V_CALL}:                V region allele assignments.
#'     \item  \code{V_CALL_GENOTYPED}:      TIgGER corrected V region allele assignment.
#'     \item  \code{D_CALL}:                D region allele assignments.
#'     \item  \code{J_CALL}:                J region allele assignments.
#'     \item  \code{JUNCTION}:              Junction region sequence.
#'     \item  \code{JUNCTION_LENGTH}:       Length of the junction region in nucleotides.
#'     \item  \code{NP1_LENGTH}:            Combined length of the N and P regions proximal
#'                                          to the V region.
#'     \item  \code{NP2_LENGTH}:            Combined length of the N and P regions proximal
#'                                          to the J region.
#'     \item  \code{SAMPLE}:                Sample identifier. Time in relation to vaccination.
#'     \item  \code{ISOTYPE}:               Isotype assignment.
#'     \item  \code{DUPCOUNT}:              Copy count (number of duplicates) of the sequence.
#'     \item  \code{CLONE}:                 Change-O assignment clonal group identifier.
#' }
#' 
#' @seealso \link{ExampleDb} \link{ExampleTrees}
#' 
#' @references
#' \enumerate{
#'   \item  Laserson U and Vigneault F, et al. High-resolution antibody dynamics of 
#'            vaccine-induced immune responses. 
#'            Proc Natl Acad Sci USA. 2014 111:4928-33.
#' }
"ExampleDbChangeo"


#' Example Ig lineage trees
#'
#' A set of Ig lineage trees generated from the \code{ExampleDb} file, subset to
#' only those trees with at least four nodes.
#'
#' @format   A list of igraph objects output by \link{buildPhylipLineage}.
#'           Each node of each tree has the following annotations (vertex attributes):
#'   \itemize{
#'     \item  \code{sample_id}:          Sample identifier(s). Time in relation to vaccination.
#'     \item  \code{c_call}:             Isotype assignment(s). 
#'     \item  \code{duplication_count}:  Copy count (number of duplicates) of the sequence.
#'   }
#'   
#' @seealso \link{ExampleTrees}
"ExampleTrees"


#### Constants ####

#' Default colors
#' 
#' Default color palettes for DNA characters, Ig isotypes, and TCR chains.
#' 
#' @format  Named character vectors with hexcode colors as values.
#' \itemize{
#'   \item  \code{DNA_COLORS}:  DNA character colors 
#'                              \code{c("A", "C", "G", "T")}.
#'   \item  \code{IG_COLORS}:   Ig isotype colors 
#'                              \code{c("IGHA", "IGHD", "IGHE", "IGHG", "IGHM", "IGHK", "IGHL")}.
#'   \item  \code{TR_COLORS}:   TCR chain colors 
#'                              \code{c("TRA", "TRB", "TRD", "TRG")}.
#' }
#' 
#' @examples 
#' # IG_COLORS as an isotype color set for ggplot
#' isotype <- c("IGHG", "IGHM", "IGHM", "IGHA")
#' db <- data.frame(x=1:4, y=1:4, iso=isotype)
#' g1 <- ggplot(db, aes(x=x, y=y, color=iso)) + 
#'     scale_color_manual(name="Isotype", values=IG_COLORS) +
#'     geom_point(size=10)
#' plot(g1)
#' 
#' # DNA_COLORS to translate nucleotide values to a vector of colors 
#' # for use in base graphics plots
#' seq <- c("A", "T", "T", "C")
#' colors <- translateStrings(seq, setNames(names(DNA_COLORS), DNA_COLORS))
#' plot(1:4, 1:4, col=colors, pch=16, cex=6)
#' 
#' @name   DEFAULT_COLORS
NULL

#' @rdname   DEFAULT_COLORS
#' @export
DNA_COLORS <- c("A"="#64F73F", 
                "C"="#FFB340", 
                "G"="#EB413C", 
                "T"="#3C88EE")

#' @rdname DEFAULT_COLORS
#' @export
IG_COLORS <- c("IGHA"="#377EB8", 
               "IGHD"="#FF7F00", 
               "IGHE"="#E41A1C", 
               "IGHG"="#4DAF4A", 
               "IGHM"="#984EA3",
               "IGHK"="#E5C494",
               "IGHL"="#FFD92F")

#' @rdname DEFAULT_COLORS
#' @export
TR_COLORS <- c("TRA"="#CBD5E8", 
               "TRB"="#F4CAE4", 
               "TRD"="#FDCDAC", 
               "TRG"="#E6F5C9")

#' IUPAC ambiguous characters
#'
#' A translation list mapping IUPAC ambiguous characters code to corresponding nucleotide
#' amino acid characters.
#' 
#' @format  A list with single character codes as names and values containing character 
#'          vectors that define the set of standard characters that match to each each 
#'          ambiguous character.
#' \itemize{
#'   \item  \code{IUPAC_DNA}:  DNA ambiguous character translations.
#'   \item  \code{IUPAC_AA}:   Amino acid ambiguous character translations.
#' }
#' 
#' @name    IUPAC_CODES
NULL

#' @rdname  IUPAC_CODES
#' @export
IUPAC_DNA <- list("A"="A", 
                  "C"="C", 
                  "G"="G", 
                  "T"="T",
                  "M"=c("A","C"), 
                  "R"=c("A","G"), 
                  "W"=c("A","T"), 
                  "S"=c("C","G"), 
                  "Y"=c("C","T"), 
                  "K"=c("G","T"), 
                  "V"=c("A","C","G"), 
                  "H"=c("A","C","T"), 
                  "D"=c("A","G","T"), 
                  "B"=c("C","G","T"),
                  "N"=c("A","C","G","T"))

#' @rdname    IUPAC_CODES
#' @export
IUPAC_AA <-  list("A"="A", 
                  "B"=c("N","R"),
                  "C"="C", 
                  "D"="D",
                  "E"="E",
                  "F"="F",
                  "G"="G",
                  "H"="H",
                  "I"="I",
                  "J"=c("I","L"),
                  "K"="K",
                  "L"="L",
                  "M"="M",
                  "N"="N",
                  "P"="P",
                  "Q"="Q",
                  "R"="R",
                  "S"="S",
                  "T"="T",
                  "V"="V",
                  "W"="W",
                  "X"=c("A","B","C","D","E","F","G","H",
                        "I","J","K","L","M","N","P","Q",
                        "R","S","T","V","W","X","Y","Z",
                        "*"),
                  "Y"="Y",
                  "Z"=c("E","Q"),
                  "*"="*")


#' Amino acid abbreviation translations
#' 
#' Mappings of amino acid abbreviations.
#' 
#' @format  Named character vector defining single-letter character codes to 
#'          three-letter abbreviation mappings.
#' 
#' @name   ABBREV_AA
#' 
#' @examples 
#' aa <- c("Ala", "Ile", "Trp")
#' translateStrings(aa, ABBREV_AA)
#' 
#' @export
ABBREV_AA <- c("A"="Ala",
               "R"="Arg",
               "N"="Asn",
               "D"="Asp",
               "C"="Cys",
               "Q"="Gln",
               "E"="Glu",
               "G"="Gly",
               "H"="His",
               "I"="Ile",
               "L"="Leu",
               "K"="Lys",
               "M"="Met",
               "F"="Phe",
               "P"="Pro",
               "S"="Ser",
               "T"="Thr",
               "W"="Trp",
               "Y"="Tyr",
               "V"="Val")


#' IMGT V-segment regions
#'
#' A list defining the boundaries of V-segment framework regions (FWRs) and complementarity 
#' determining regions (CDRs) for IMGT-gapped immunoglobulin (Ig) nucleotide sequences 
#' according to the IMGT numbering scheme.
#' 
#' @format  A list with regions named one of \code{c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")} 
#'          with values containing a numeric vector of length two defining the 
#'          \code{c(start, end)} positions of the named region.
#'          
#' @references
#'   \url{http://imgt.org}
#' 
#' @export
IMGT_REGIONS <- list("fwr1"=c(1, 78),
                     "cdr1"=c(79, 114),
                     "fwr2"=c(115, 165),
                     "cdr2"=c(166, 195),
                     "fwr3"=c(196, 312))