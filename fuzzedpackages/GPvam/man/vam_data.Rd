\name{vam_data}
\alias{vam_data}
\docType{data}
\title{
Simulated Data
}
\description{
A simulated data set used to illustrate the functionality of the package. The data are simulated according to the VP model, and demonstrate the stability of the program in the presence of perfectly correlated future year effects.
}
\usage{data(vam_data)}
\format{
  A data frame with 3750 observations on 1250 students over 3 years, with 50 teachers in each year. The data set contains the following 5 variables.
  \describe{
    \item{\code{y}}{a numeric vector representing the student score}
    \item{\code{student}}{a numeric vector}
    \item{\code{year}}{a numeric vector}
    \item{\code{teacher}}{a numeric vector}
    \item{\code{cont_var}}{a numeric vector representing a continuous covariate}
  }
}
\details{
The data set may be reproduced with the following code.

set.seed(0) \cr
years<-3 \cr
#teacher in each year \cr
teachers<-50 \cr
#students in each class \cr
students<-25 \cr
alpha<-.4 \cr
eta.stu<-rnorm(students*teachers,0,5) \cr
z1<-rep(1:teachers,each=students) \cr
z2<-sample(rep(1:teachers,each=students)) \cr
z3<-sample(rep(1:teachers,each=students)) \cr
cont_var1<-rnorm(students*teachers,0,4) \cr
cont_var2<-rnorm(students*teachers,0,4) \cr
cont_var3<-rnorm(students*teachers,0,4) \cr
gam1<- rnorm(teachers,0,5) \cr
gam2<- rnorm(teachers,0,5) \cr
gam3<- rnorm(teachers,0,5) \cr
eps1<- rnorm(students*teachers,0,5) \cr
eps2<- rnorm(students*teachers,0,5) \cr
eps3<- rnorm(students*teachers,0,5) \cr
y1<-eta.stu+gam1[z1]+cont_var1+eps1 \cr
y2<-eta.stu+gam1[z1]*alpha+gam2[z2]+cont_var2+eps2 \cr
y3<-eta.stu+gam1[z1]*alpha+gam2[z2]*alpha+gam3[z3]+cont_var3+eps3 \cr
student<-1:(students*teachers) \cr
teacher<-c(z1,z2,z3) \cr
cont_var<-c(cont_var1,cont_var2,cont_var3) \cr
year<-c(rep(1:3,each=students*teachers)) \cr
y<-c(y1,y2,y3) \cr
vam_data<-as.data.frame(cbind(student,teacher,year,y,cont_var)) \cr

}

\examples{
data(vam_data)
print(vam_data[1,])
}
\keyword{datasets}
