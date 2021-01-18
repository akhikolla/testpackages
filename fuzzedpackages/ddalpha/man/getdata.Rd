\name{getdata}
\alias{getdata}
\alias{data}

\alias{baby}
\alias{banknoten}
\alias{biomed}
\alias{bloodtransfusion}
\alias{breast_cancer_wisconsin}
\alias{bupa}
\alias{chemdiab_1vs2}
\alias{chemdiab_1vs3}
\alias{chemdiab_2vs3}
\alias{cloud}
\alias{crabB_MvsF}
\alias{crabF_BvsO}
\alias{crabM_BvsO}
\alias{crabO_MvsF}
\alias{crab_BvsO}
\alias{crab_MvsF}
\alias{cricket_CvsP}
\alias{diabetes}
\alias{ecoli_cpvsim}
\alias{ecoli_cpvspp}
\alias{ecoli_imvspp}
\alias{gemsen_MvsF}
\alias{glass}
\alias{groessen_MvsF}
\alias{haberman}
\alias{heart}
\alias{hemophilia}
\alias{indian_liver_patient_1vs2}
\alias{indian_liver_patient_FvsM}
\alias{iris_setosavsversicolor}
\alias{iris_setosavsvirginica}
\alias{iris_versicolorvsvirginica}
\alias{irish_ed_MvsF}
\alias{kidney}
\alias{pima}
\alias{plasma_retinol_MvsF}
\alias{segmentation}
\alias{socmob_IvsNI}
\alias{socmob_WvsB}
\alias{tae}
\alias{tennis_MvsF}
\alias{tips_DvsN}
\alias{tips_MvsF}
\alias{uscrime_SvsN}
\alias{vertebral_column}
\alias{veteran_lung_cancer}
\alias{vowel_MvsF}
\alias{wine_1vs2}
\alias{wine_1vs3}
\alias{wine_2vs3}

\title{
Data for Classification
}
\description{
50 multivariate data sets for binary classification. For more details refer {\url{http://www.wisostat.uni-koeln.de/de/forschung/software-und-daten/data-for-classification/}}

The \code{getdata} function gets the data set from the package, and returns it. The dataset itself does not appear in the global environment and the existing variables with the same name remain unchanged.
}
\usage{
# load the data set
# data(name)

# load the data set by name
# data(list = "name")

# load the data set by name to a variable
# getdata("name")

}
\arguments{
  \item{name}{
the data set name. 
}
}
\format{
  A data frame with \code{n} observations on the \code{d} variables. The last \code{d+1} column is the class label.
  \describe{
    \item{\code{x[,1:d]}}{numeric values}
    \item{\code{x[,d+1]}}{the numeric class label (0 or 1) or (1 or 2)}
  }
}
\details{
The package contains data sets used in the joint project of the University of Cologne and the Hochschule Merseburg "Classifying real-world data with the DDalpha-procedure". Comprehensive description of the methodology, and experimental settings and results of the study are presented in the work:

Mozharovskyi, P., Mosler, K., and Lange, T. (2015). Classifying real-world data with the DD\eqn{\alpha}-procedure. \emph{Advances in Data Analysis and Classification} \bold{9} 287--314.

For a more complete explanation of the technique and further experiments see:
Lange, T., Mosler, K., and Mozharovskyi, P. (2014). Fast nonparametric classification based on data depth. \emph{Statistical Papers} \bold{55} 49--69.

50 binary classification tasks have been obtained from partitioning 33 freely accessible data sets. Multiclass problems were reasonably split into binary classification problems, some of the data set were slightly processed by removing objects or attributes and selecting prevailing classes. Each data set is provided with a (short) description and brief descriptive statistics. The name reflects the origination of the data. A letter after the name is a property filter, letters (also their combinations) in brackets separated by "vs" are the classes opposed. The letters (combinations or words) stand for labels of classes (names of properties) and are intuitive. Each description contains a link to the original data.

The data have been collected as open source data in January 2013. Owners of the package decline any responsibility regarding their correctness or consequences of their usage. If you publish material based on these data, please quote the original source. Special requests regarding citations are found on data set's web page.

}
\references{
Lange, T., Mosler, K., and Mozharovskyi, P. (2014). Fast nonparametric classification based on data depth. \emph{Statistical Papers} \bold{55} 49--69.

Mozharovskyi, P., Mosler, K., and Lange, T. (2015). Classifying real-world data with the DD\eqn{\alpha}-procedure. \emph{Advances in Data Analysis and Classification} \bold{9} 287--314.

The general list of sources consists of:

UCI Machine Learning Repository {\url{http://archive.ics.uci.edu/ml}}\cr
R-packages {\url{https://CRAN.R-project.org/}} \cr
{\url{http://stat.cmu.edu}} \cr
{\url{http://stat.ethz.ch/Teaching/Datasets}} \cr
{\url{http://www.stats.ox.ac.uk/pub/PRNN}} 
}

\seealso{
\code{\link[utils:data]{utils:data}}
}

\note{
List of the datasets:

baby\cr
banknoten\cr
biomed\cr
bloodtransfusion\cr
breast_cancer_wisconsin\cr
bupa\cr
chemdiab_1vs2\cr
chemdiab_1vs3\cr
chemdiab_2vs3\cr
cloud\cr
crabB_MvsF\cr
crabF_BvsO\cr
crabM_BvsO\cr
crabO_MvsF\cr
crab_BvsO\cr
crab_MvsF\cr
cricket_CvsP\cr
diabetes\cr
ecoli_cpvsim\cr
ecoli_cpvspp\cr
ecoli_imvspp\cr
gemsen_MvsF\cr
glass\cr
groessen_MvsF\cr
haberman\cr
heart\cr
hemophilia\cr
indian_liver_patient_1vs2\cr
indian_liver_patient_FvsM\cr
iris_setosavsversicolor\cr
iris_setosavsvirginica\cr
iris_versicolorvsvirginica\cr
irish_ed_MvsF\cr
kidney\cr
pima\cr
plasma_retinol_MvsF\cr
segmentation\cr
socmob_IvsNI\cr
socmob_WvsB\cr
tae\cr
tennis_MvsF\cr
tips_DvsN\cr
tips_MvsF\cr
uscrime_SvsN\cr
vertebral_column\cr
veteran_lung_cancer\cr
vowel_MvsF\cr
wine_1vs2\cr
wine_1vs3\cr
wine_2vs3\cr


Also functional data sets can be loaded:

geneexp\cr
growth\cr
medflies\cr
population\cr
population2010\cr
tecator

}


\examples{
# load a dataset using data()
data(hemophilia)
data(list = "hemophilia")
rm(hemophilia)

# load data set using getdata()
hemophilia = "This is some existing object called 'hemophilia'. It remains unchanged"
d = getdata("hemophilia")
head(d)
print(hemophilia)

#get the list of all data sets
names = data(package = "ddalpha")$results[,3]

}


\keyword{datasets}