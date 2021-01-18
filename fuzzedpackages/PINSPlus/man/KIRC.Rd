\name{KIRC}
\docType{data}
\alias{KIRC}
\title{Kidney renal clear cell carcinoma dataset}
\description{
    The Cancer Genome Atlas Kidney Renal Clear Cell Carcinoma (TCGA-KIRC) data collection is part of a larger effort to build a research community focused on connecting cancer phenotypes to genotypes by providing clinical images matched to subjects from The Cancer Genome Atlas (TCGA). Clinical, genetic, and pathological data resides in the Genomic Data Commons (GDC) Data Portal while the radiological data is stored on The Cancer Imaging Archive (TCIA).
    
    This embed version of KIRC in PINPlus package is the reduced version of KIRC using Principle Component Analysis.
}
\format{
    A list containing properties:
    \tabular{rlll}{
        Name          \tab Type       \tab Description\cr
        GE            \tab data.frame \tab mRNA expression  data\cr
        ME            \tab data.frame \tab DNA Methylation data\cr
        MI            \tab data.frame \tab miRNA expression  data\cr
        survival      \tab data.frame \tab Clinical survival data\cr
    }
}
\source{\url{https://wiki.cancerimagingarchive.net/display/Public/TCGA-KIRC}}
\references{
    Akin, O., Elnajjar, P., Heller, M., Jarosz, R., Erickson, B. J., Kirk, S., ... Filippini, J. (2016). Radiology Data from The Cancer Genome Atlas Kidney Renal Clear Cell Carcinoma [TCGA-KIRC] collection. The Cancer Imaging Archive.
}
\keyword{datasets}