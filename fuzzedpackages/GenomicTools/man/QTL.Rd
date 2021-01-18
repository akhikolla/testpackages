\name{QTL}
\alias{QTL}

\title{
  Perform a QTL Analysis
}

\description{
  This function performs a QTL analysis.
}
\usage{
  QTL(pheno, phenoSamples=NULL, geno=NULL, genoSamples=NULL, 
      method="LM", mc=1, sig=NULL, testType="asymptotic", 
      nper=2000, which=NULL, verbose=TRUE)
}

\arguments{
  \item{pheno}{Matrix or Vector with phenotype values.}
  \item{phenoSamples}{Sample names for the phenotype values, see details (optional).}
  \item{geno}{Genotype data.}
  \item{genoSamples}{Sample names for the genotype values, see details (optional).}
  \item{method}{Method of choice for the QTL, see details.}
  \item{mc}{Amount of cores for parallel computing.}
  \item{sig}{Significance level for the QTL testing, see details.}
  \item{testType}{Type of significance test, see details.}
  \item{nper}{Sets the amount of permutations, if permuation tests are used.}
  \item{which}{Names of phenotypes for that the QTL should be performed.}
  \item{verbose}{Logical, if the method should report intermediate results.}
}

\details{
This function performs a QTL analysis and offers different types of tests. The type of test
can be specified with the \code{method} option and possible options are \code{"LM"} and \code{"directional"}. 
The option \code{"LM"} fits for each provided SNP a linear model for the genotype information and the corresponding
phenotype values. The null hypothesis for each test is then that the slope is equal to zero and the alternative is that it is not zero.

The \code{"directional"} option applies a new directional test based on probabilistic indices for triples as described 
in Fischer, Oja, et al. (2014). Being \eqn{\mathbf{x}_0=(x_{01},x_{02},\dots,x_{0N_0})', \mathbf{x}_1=(x_{11},x_{12},\dots,x_{1N_1})'}
and \eqn{\mathbf{x}_2=(x_{21},x_{22},\dots,x_{2N_2})'} the phenotype values that are linked to the three genotype
groups \eqn{0,1} and \eqn{2} with underlying distributions \eqn{F_0, F_1} and \eqn{F_2}.  We first calculate the probabilistic
indices \eqn{P_{0,1,2} = \frac{1}{N_0 N_1 N_2} \sum_i \sum_j \sum_k I(x_{0i} < x_{1j} < x_{2k})}
and \eqn{P_{2,1,0} = \frac{1}{N_0 N_1 N_2} \sum_i \sum_j \sum_k I(x_{2i} < x_{1j} < x_{0k})}. These are the probabilities that the phenotype
values of the three groups follow a certain order what we would expect for possible QTLs. The null hypothesis that we have then
in mind is that the phenotype values from these three groups have the same distribution \eqn{H_0: F_0 = F_1 = F_2} and the
two alternatives are that the distributions have a certain stochastical order \eqn{H_1: F_0 < F_1 < F_2} and \eqn{H_2: F_2 < F_1 < F_0}.

 The test is applied for the two probabilistic indices \eqn{P_{0,1,2}} and
\eqn{P_{2,1,0}} and combines the two resulting p-values \eqn{p_{012}=p_1} and \eqn{p_{210}=p_2} from previous tests then as overall p-value
\eqn{\min(2 \min(p_1 , p_2 ), 1)}. In the two-group case (this means only two different genotypes are present for a certain SNP) a two-sided
Wilcoxon rank-sum test is applied.

The phenotype values used in the QTL are specified in \code{pheno}. If several phenotypes should be tested, then \code{pheno} is a matrix and each
column refers to a phenotype and each row to an individuum. Sample names can either be given as row names in the matrix or as separate
vector in \code{phenoSamples}. If only values of one phenotype should be tested then \code{pheno} can be a vector.

The genotype information is provided in the \code{geno} object. Here one can i.a. specify the file name of a ped/map file pair.
The function then imports the genotype information using the function \code{importPED}. In case the genotype information has been imported
already earlier using \code{importPED} the resulting \code{PedMap} object can also be given as a parameter for \code{geno}.

The option \code{genoSamples} is used in case that the sample names in the ped/map file (or SnpMatrix)
do not match with \code{rownames(pheno)} given in the expression matrix. The vector \code{genoSamples} is as
long as the \code{geno} object has samples, but gives then for each row in \code{geno} the corresponding
name in the \code{pheno} object. The function finds then also the smallest union between the two data objects.
If there are repeated measurements per individual for the genotypes we take by default only the first
appearance in the data and neglect all successive values. Currently this cannot be changed. In case
this behavior is not desired, the user has to remove the corresponding rows from \code{geno} before starting the 
calculation.

If the code is executed on a Linux OS the user can specify with the \code{mc} option the amount of CPU cores used for the
calculation. 

If the \code{sig} option is set to a certain significance level, then the method only reports those SNPs that are tested to be
significant. This can reduce the required memory drastically, but shouldn't be used in case the results will be later plotted as Manhattan plot. 
}

\value{
 A list of class \code{qtl} containing the values 
 \item{pheno}{The \code{pheno} object from the function call.}
 \item{geno}{The \code{geno} object from the function call.}
 \item{genoSamples}{The \code{genoSamples} object from the function call.}
 \item{windowSize}{The \code{windowSize} object from the function call.}

 and an incapsulated list \code{qtl} where each list item is a tested phenotype and contains the items

 \item{ProbeLoc}{Used position of that gene. (Only different from 1 if multiple locations are considered.)}
 \item{TestedSNP}{Details about the considered SNPs.}
 \item{p.values}{P values of the test.}
 \item{GeneInfo}{Details about the center gene.}
}

\references{
  \cite{Fischer, D., Oja, H., Sen, P.K., Schleutker, J., Wahlfors, T. (2014): Generalized Mann-Whitney Type Tests for Microarray Experiments, Scandinavian Journal of Statistics, 3, pages 672-692, 
           doi: 10.1111/sjos.12055.}

  \cite{Fischer, D., Oja, H. (2013): Mann-Whitney Type Tests for Microarray Experiments: The R Package gMWT, submitted article.}
}


\author{
  Daniel Fischer
}

\keyword{
  methods
}
 
