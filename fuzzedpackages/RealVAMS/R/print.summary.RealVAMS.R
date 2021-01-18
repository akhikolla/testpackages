print.summary.RealVAMS <-
function (x, ...) 
{
  count.table<-table(table(x$joined.table[x$joined.table$response=="score","student"]))
    cat("Number of test scores:", sum(x$y.response=="score"), "\n")
    cat("Number of students (bottom row) with a given number of test scores recorded (top row)")
    print(count.table)
  cat("Number of outcome indicators recorded:", sum(x$y.response=="outcome"), "\n")
    cat("Number of years:", x$num.year, "\n")
    cat("Number of students:", x$num.student, "\n")
    for (i in 1:x$num.year) {
        cat("Number of teachers in year", i, ":", x$num.teach[i], 
            "\n")
    }
    cat("\n")
  #  cat("Number of EM iterations: ", x$iter, "\n")
  #  cat("-2 log-likelihood", -2 * x$loglik, "\n")
  #  cat("AIC", x$AIC, "\n")
  #  cat("AICc", x$AICc, "\n")
    cat("Outcome family:\n")
    print(x$outcome.family)
    cat("Persistence: ",x$persistence,"\n")
    cat("Persistence Parameters (of year [column index] teacher on year [row index] score): \n")
    print(x$persistence_parameters)
    
    cat("\n")
    for (i in 1:x$num.year) {
        cat("Covariance matrix for year", i, "teachers.\n")
        print(as.matrix(x$teach.cov[[i]]))
        cat("\n")
        cat("with correlation matrix\n")
        print(round(cov2cor(as.matrix(x$teach.cov[[i]])), 4))
        cat("\n")
    }
    if (!any(is.na(x$R_i))) {
        cat("Block of error covariance matrix (R):\n")
        print(as.matrix(x$R_i))
        cat("with correlation matrix\n")
        print(round(cov2cor(as.matrix(x$R_i)), 4))
        cat("\n")
    }
   # if (!is.na(x$stu.cov)) {
  #      cat("Student variance component:\n")
  #      print(x$stu.cov)
  #      cat("\n")
  #  }
    cat("Parameter estimates:\n")
    print(x$parameters)
    cat("\n")
    cat("Distribution of raw marginal residuals for test scores\n")
    print(summary(x$mresid))
    cat("\n")
    cat("Distribution of raw conditional residuals for test scores\n")
    print(summary(x$cresid))
    cat("\n")
  
    
   # cat("Distribution of scaled conditional residuals\n")
  #  print(summary(x$sresid))
  #  cat("\n")

   
}
