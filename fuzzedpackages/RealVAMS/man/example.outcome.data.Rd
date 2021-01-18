\name{example.outcome.data}
\alias{example.outcome.data}
\docType{data}
\title{
Simulated Data
}
\description{
A simulated data set used to illustrate the functionality of the package. This data set represents binary outcome measurements on 625 students (with one missing).
}
\usage{data(example.outcome.data)}
\format{
  A data frame with 624 observations. The data set contains the following 2 variables.
  \describe{
    \item{\code{r}}{a numeric vector composed of 0's and 1's representing a binary outcome measured on students.}
    \item{\code{student}}{a numeric vector}
  }
}
\details{
The data set may be reproduced with the following code.
\preformatted{
set.seed(0)                                                                                                                                                         
library(MASS)
#number of years: fixed at 3 for this code block
years<-3   
#teacher in each year
teachers<-25
#students in each class
students<-25
alpha<-.5
eta.stu.j <- 
mvrnorm(n=teachers*students,mu=c(0,0),Sigma=cbind(c(5,.2),c(.2,.1)))
eta.stu<-eta.stu.j[,1]
eta.stu.r<-  eta.stu.j[,2]
z1<-rep(1:teachers,each=students)
z2<-sample(rep(1:teachers,each=students))
z3<-sample(rep(1:teachers,each=students))
cont_var1<-rnorm(students*teachers,0,.5)
cont_var2<-rnorm(students*teachers,0,.5)
cont_var3<-rnorm(students*teachers,0,.5)
gam <- mvrnorm(n=teachers*years,mu=c(0,0),Sigma=cbind(c(5,.6),c(.6,.6)))
eps1<- rnorm(students*teachers,0,sqrt(5))
eps2<- rnorm(students*teachers,0,sqrt(5))
eps3<- rnorm(students*teachers,0,sqrt(5))
gam1<-gam[seq(1,teachers),1]
gam2<-gam[seq((teachers+1),(2*teachers)),1]
gam3<-gam[seq((2*teachers+1),(3*teachers)),1]
gam1.r<-gam[seq(1,teachers),2]
gam2.r<-gam[seq((teachers+1),(2*teachers)),2]
gam3.r<-gam[seq((2*teachers+1),(3*teachers)),2]
y1<-50+ eta.stu+gam1[z1]+cont_var1+eps1
y2<-eta.stu+gam1[z1]*alpha+gam2[z2]+cont_var2+eps2
y3<-100+eta.stu+gam1[z1]*alpha+gam2[z2]*alpha+gam3[z3]+cont_var3+eps3
r1<-rbinom(students*teachers,1,
           pnorm(.1+eta.stu.r+gam1.r[z1]+gam2.r[z2]+gam3.r[z3]))
student<-1:(students*teachers)
teacher<-c(z1,z2,z3)
cont_var<-c(cont_var1,cont_var2,cont_var3)
year<-c(rep(1:3,each=students*teachers))
y<-c(y1,y2,y3)
vam_data2<-as.data.frame(cbind(student,teacher,year,y,cont_var))
vam_data2<-vam_data2[order(vam_data2$student,vam_data2$year),]
example.score.data<-vam_data2
vam_data2.r<-as.data.frame(cbind(student,r=r1))
example.outcome.data<-vam_data2.r
}
}
\examples{
data(example.outcome.data)
print(example.outcome.data[1,])
}
\keyword{datasets}
