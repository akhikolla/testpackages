rococo <- function(x, y, similarity=c("linear", "exp", "gauss", "epstol",
                                      "classical"),
                   tnorm="min", r=0, noVarReturnZero=TRUE)
{
     if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y))
          stop("'x' and 'y' need to be numeric vectors of the same length")

     if (!is.numeric(r) || any(r < 0))
          stop("'r' must be a (vector of) non-negative real number(s)")

     if (missing(similarity))
         similarity <- "linear"

     similarity <- match.arg(similarity, several.ok=TRUE)

     if (length(r) == 1)
         r[2] <- r[1]

     tnormSel <- -1

     if (is.character(tnorm))
     {
         tnorm <- match.arg(tnorm, c("min", "prod", "lukasiewicz"))

         tnormSel <- switch(tnorm, min=1, prod=2, lukasiewicz=3)
     }
     else if (is.function(tnorm))
     {
          if (length(formals(tnorm)) != 2)
               stop("'tnorm' should be a function of two arguments, ",
                    "e.g. 'tnorm=function(a, b) a * b' for the product t-norm")

          if (abs(tnorm(1, 0.5) - 0.5) > .Machine$double.eps ||
              abs(tnorm(0, 0.25)) > .Machine$double.eps)
               stop("supplied function does not appear to be a valid t-norm")

          if (requireNamespace("compiler", quietly=TRUE))
               tnorm <- compiler::cmpfun(tnorm)
     }
     else
         stop("'tnorm' should be valid string or a function of two arguments, ",
              "e.g. 'tnorm=function(a, b) a * b' for the product t-norm")

     if (length(similarity) > 1 && similarity[1] != similarity[2])
     {
         if (similarity[1] != "classical" && r[1] == 0)
         {
             r[1] <- 0.1 * IQR(x)

             if (r[1] == 0)
             {
                 warning("IQR(x) = 0 => ",
                         "could not determine tolerance; ",
                         "reverting to classical crisp similarity")
                 similarity[1] <- "classical"
             }
         }

         if (similarity[2] != "classical" && r[2] == 0)
         {
             r[2] <- 0.1 * IQR(y)

             if (r[2] == 0)
             {
                 warning("IQR(y) = 0 => ",
                         "could not determine tolerance; ",
                         "reverting to classical crisp similarity")
                 similarity[2] <- "classical"
             }
         }

         Rx <- switch(similarity[1],
                      linear=.Call("rcor_matrix_linear",
                                   vx=as.double(x), as.double(r[1])),
                      exp=.Call("rcor_matrix_exp",
                                vx=as.double(x), as.double(r[1])),
                      gauss=.Call("rcor_matrix_gauss",
                                  vx=as.double(x), as.double(r[1])),
                      epstol=.Call("rcor_matrix_epstol",
                                   vx=as.double(x), as.double(r[1])),
                      classical=.Call("rcor_matrix_classical",
                                   vx=as.double(x), as.double(r[1])))
         Ry <- switch(similarity[2],
                      linear=.Call("rcor_matrix_linear",
                                   vx=as.double(y), as.double(r[2])),
                      exp=.Call("rcor_matrix_exp",
                                vx=as.double(y), as.double(r[2])),
                      gauss=.Call("rcor_matrix_gauss",
                                  vx=as.double(y), as.double(r[2])),
                      epstol=.Call("rcor_matrix_epstol",
                                   vx=as.double(y), as.double(r[2])),
                      classical=.Call("rcor_matrix_classical",
                                   vx=as.double(y), as.double(r[2])))
     }
     else
     {
         if (similarity[1] != "classical")
         {
             if (r[1] == 0)
             {
                 r[1] <- 0.1 * IQR(x)

                 if (r[1] == 0)
                 {
                     warning("IQR(x) = 0 => ",
                             "could not determine tolerance; ",
                             "reverting to classical crisp similarity")
                     similarity[1] <- "classical"
                 }
             }

             if (r[2] == 0)
             {
                 r[2] <- 0.1 * IQR(y)

                 if (r[2] == 0)
                 {
                     warning("IQR(y) = 0 => ",
                             "could not determine tolerance; ",
                             "reverting to classical crisp similarity")
                     similarity[2] <- "classical"
                 }
             }
         }

         if (length(similarity) > 1 && similarity[1] != similarity[2])
         {
             Rx <- switch(similarity[1],
                          linear=.Call("rcor_matrix_linear",
                                       vx=as.double(x), as.double(r[1])),
                          exp=.Call("rcor_matrix_exp",
                                    vx=as.double(x), as.double(r[1])),
                          gauss=.Call("rcor_matrix_gauss",
                                      vx=as.double(x), as.double(r[1])),
                          epstol=.Call("rcor_matrix_epstol",
                                       vx=as.double(x), as.double(r[1])),
                          classical=.Call("rcor_matrix_classical",
                                          vx=as.double(x), as.double(r[1])))
             Ry <- switch(similarity[2],
                          linear=.Call("rcor_matrix_linear",
                                       vx=as.double(y), as.double(r[2])),
                          exp=.Call("rcor_matrix_exp",
                                    vx=as.double(y), as.double(r[2])),
                          gauss=.Call("rcor_matrix_gauss",
                                      vx=as.double(y), as.double(r[2])),
                          epstol=.Call("rcor_matrix_epstol",
                                       vx=as.double(y), as.double(r[2])),
                          classical=.Call("rcor_matrix_classical",
                                          vx=as.double(y), as.double(r[2])))
         }
         else
         {
             matrices <- switch(similarity[1],
                                linear=.Call("rcor_matrices_linear",
                                             vx=as.double(x),
                                             vy=as.double(y),
                                             as.double(r[1]),
                                             as.double(r[2])),
                                exp=.Call("rcor_matrices_exp",
                                          vx=as.double(x),
                                          vy=as.double(y),
                                          as.double(r[1]),
                                          as.double(r[2])),
                                gauss=.Call("rcor_matrices_gauss",
                                            vx=as.double(x),
                                            vy=as.double(y),
                                            as.double(r[1]),
                                            as.double(r[2])),
                                epstol=.Call("rcor_matrices_epstol",
                                             vx=as.double(x),
                                             vy=as.double(y),
                                             as.double(r[1]),
                                             as.double(r[2])),
                                classical=.Call("rcor_matrices_classical",
                                                vx=as.double(x),
                                                vy=as.double(y),
                                                as.double(r[1]),
                                                as.double(r[2])))

             Rx <- matrices$Rx
             Ry <- matrices$Ry
         }
     }

     if (tnormSel > 0)
     {
          res <- .Call("rcor", Rx, Ry, as.integer(tnormSel))
          c <- res$c
          d <- res$d
     }
     else
     {
          c <- sum(mapply(tnorm, Rx, Ry))
          d <- sum(mapply(tnorm, Rx, t(Ry)))
     }

     if (!is.na(c) && !is.na(d) && c + d > 0)
         gamma <- (c - d) / (c + d)
     else if (identical(noVarReturnZero, TRUE))
         gamma <- 0
     else
     {
         gamma <- NA
         warning("no variation in at least one of the two observables")
     }

     gamma
}
