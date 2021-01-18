stat2 <- function(x){ # Anderson-Darling
    x <- sort(x)
    n <- length(x)
    logp1 <- pnorm((x - mean(x)) / sd(x), log.p = TRUE)
    logp2 <- pnorm(-(x - mean(x)) / sd(x), log.p = TRUE)
    h <- (2 * seq(1:n) - 1) * (logp1 + rev(logp2))
    A <- -n - mean(h)
    AA <- (1 + 0.75 / n + 2.25 / n ^ 2) * A
    if (AA < 0.2) {
        pval <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA ^ 2)
    }
    else if (AA < 0.34) {
        pval <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA ^ 2)
    }
    else if (AA < 0.6) {
        pval <- exp(0.9177 - 4.279 * AA - 1.38 * AA ^ 2)
    }
    else if (AA < 10) {
        pval <- exp(1.2937 - 5.709 * AA + 0.0186 * AA ^ 2)
    }
    else pval <- 3.7e-24
    res <- list(stat = A, p.value = pval) 
    return(res)
}

stat3 <- function(x){ # Zhang Zc
    n <- length(x)
    phi <- pnorm(sort((x - mean(x)) / sd(x)))
    res <- sum((log((1 / phi - 1) / ((n - 0.5) / ((1:n) - 0.75) - 1))) ^ 2)
    return(res)
}

stat4 <- function(x) { # Zhang Za
    n <- length(x)
    phi <- pnorm(sort((x - mean(x)) / sd(x)))
    res <- -sum(((log(phi)) / (n - (1:n) + 0.5)) + ((log(1.0 - phi)) / ((1:n) - 0.5)))
    return(10.0 * res - 32.0) # See page 711 in their paper
}

stat6 <- function(x) { # D'Agostino-Pearson. Credits to: fBasics package, dagoTest()
    n <- length(x)
    meanX <- mean(x)
    s <- sqrt(mean((x - meanX) ^ 2))
    a3 <- mean((x - meanX) ^ 3)/s ^ 3
    a4 <- mean((x - meanX) ^ 4) / s ^ 4
    SD3 <- sqrt(6 * (n - 2) / ((n + 1) * (n + 3)))
    SD4 <- sqrt(24 * (n - 2) * (n - 3) * n / ((n + 1) ^ 2 * (n + 3) * 
                                              (n + 5)))
    U3 <- a3 / SD3
    U4 <- (a4 - 3 + 6 / (n + 1)) / SD4
    b <- (3 * (n ^ 2 + 27 * n - 70) * (n + 1) * (n + 3)) / ((n - 2) * 
                                                            (n + 5) * (n + 7) * (n + 9))
    W2 <- sqrt(2 * (b - 1)) - 1
    delta <- 1 / sqrt(log(sqrt(W2)))
    a <- sqrt(2 / (W2 - 1))
    Z3 <- delta * log((U3 / a) + sqrt((U3 / a) ^ 2 + 1))
    B <- (6 * (n * n - 5 * n + 2) / ((n + 7) * (n + 9))) * sqrt((6 * 
                                                                 (n + 3) * (n + 5)) / (n * (n - 2) * (n - 3)))
    A <- 6 + (8 / B) * ((2 / B) + sqrt(1 + 4 / (B ^ 2)))
    jm <- sqrt(2 / (9 * A))
    pos <- ((1 - 2 / A) / (1 + U4 * sqrt(2 / (A - 4)))) ^ (1 / 3)
    Z4 <- (1 - 2 / (9 * A) - pos) / jm
    omni <- Z3 ^ 2 + Z4 ^ 2
    pomni <- 1 - pchisq(omni, 2)
    res <- list(stat = omni, p.value = pomni) 
    return(res)
}


stat7 <- function(x){ # Jarque-Bera. Credits to: tseries package, jarque.bera.test()
    n <- length(x)
    m1 <- sum(x) / n
    m2 <- sum((x - m1) ^ 2) / n
    m3 <- sum((x - m1) ^ 3) / n
    m4 <- sum((x - m1) ^ 4) / n
    b1 <- (m3 / m2 ^ (3 / 2)) ^ 2
    b2 <- (m4 / m2 ^ 2)
    STATISTIC <- n * b1 / 6 + n * (b2 - 3) ^ 2 / 24
    PVAL <- 1 - pchisq(STATISTIC, df = 2)
    res <- list(stat = STATISTIC, p.value = PVAL) 
    return(res)
}


stat8 <- function(x){ # Doornik-Hansen. Credits to: normwhn.test package, normality.test1()
    n <- length(x)
    m1 <- sum(x) / n
    v.matrix <- 1 / sqrt(cov.wt(as.matrix(x))$cov)
    xhatprime.matrix <- t(x - m1)
    trprime <- t(v.matrix %*% xhatprime.matrix)
    v1 <- mean(trprime)
    v2 <- (n ^ (-1)) * sum((trprime - v1) ^ 2)
    v3 <- (n ^ (-1)) * sum((trprime - v1) ^ 3)
    v4 <- (n ^ (-1)) * sum((trprime - v1) ^ 4)
    rtb1 <- v3 / (v2 ^ (3 / 2))
    b2 <- v4 / v2 ^ 2
    beta <- (3 * (n ^ 2 + 27 * n - 70) * (n + 1) * (n + 3)) / ((n - 
        2) * (n + 5) * (n + 7) * (n + 9))
    w2 <- (-1) + sqrt(2 * (beta - 1))
    delta <- 1 / sqrt(log(sqrt(w2)))
    f <- (w2 - 1) / 2
    g <- (n + 1) * (n + 3) / (6 * (n - 2))
    h <- sqrt(f * g)
    y <- rtb1 * h
    z1 <- delta * log(y + sqrt(y ^ 2 + 1))
    del <- ((n - 3) * (n + 1) * (n ^ 2 + 15 * n - 4))
    aye <- ((n - 2) * (n + 5) * (n + 7) * (n ^ 2 + 27 * n - 70)) / (6 * 
        del)
    cee <- ((n - 7) * (n + 5) * (n + 7) * (n ^ 2 + 2 * n - 5)) / (6 * 
        del)
    alp <- aye + ((rtb1 ^ 2) * cee)
    kap <- ((n + 5) * (n + 7) * (n ^ 3 + 37 * n ^ 2 + 11 * n - 313)) / (12 * 
        del)
    chi <- (b2 - 1 - rtb1 ^ 2) * (2 * kap)
    chi <- abs(chi)
    z2 <- (((chi / (2 * alp)) ^ (1 / 3)) - 1 + (1 / ((9 * alp)))) * sqrt(9 * 
        alp)
    stat <- z1 ^ 2 + z2 ^ 2
    pval <- 1 - pchisq(stat, 2)
    res <- list(stat = stat, p.value = pval) 
    return(res)
}


stat17 <-function (x) { # Bonett-Seier. Credits to moments package, bonett.test()
    x <- sort(x)
    n <- length(x)
    rho <- sqrt(sum((x - mean(x)) ^ 2) / n)
    tau <- sum(abs(x - mean(x))) / n
    omega <- 13.29 * (log(rho) - log(tau))
    z <- sqrt(n + 2) * (omega - 3) / 3.54
    pval <- 2* pnorm(z, lower.tail = FALSE)
        if (pval > 1) 
            pval <- 2 - pval
    res <- list(stat = z, p.value = pval) 
    return(res)
}


stat26 <- function(x){ # Chen-Shapiro
  n <- length(x)
  m <- qnorm(((1:n) - 0.375) / (n + 0.25))
  xo <- sort(x)
  QH <- sum(diff(xo) / diff(m)) / (n - 1) / sd(x); 
  res <- sqrt(n) * (1 - QH)
  return(res)
}

stat29 <- function(x){ # BCMR
  n <- length(x)
  xo <- sort(x)
  dt <- 0.000001
  k <- 1:n
  sigma <- 0
  for (i in k){
    u <- seq((k[i] - 1) / n + dt / 2, k[i] / n - dt / 2, dt)
    sigma <- sigma +  sum(qnorm(u)) * xo[i] * dt
  }
  stat <- 1 - sigma ^ 2 / (var(x) * (n - 1) / n)
  return(stat)
}

stat30 <- function(x) { # Coin. Credits to http://digilander.libero.it/polynomial/B3test.R
    normorder <- function(n){
        dt <- 0.01
        t <- seq(-10, 10, dt)
        int <- function(r) {
            rinv <- n - r + 1
            logh <- lgamma(n + 1) - lgamma(rinv) - lgamma(n - rinv + 1) +
                (rinv - 1) * pnorm(t, low = FALSE, log.p = TRUE) + (n - rinv) * pnorm(t, low = TRUE, log.p = TRUE) +
                dnorm(t, log = TRUE) 
            h <- t * exp(logh)
            return(sum(h) * dt)
        }
        res <- apply(cbind(1:n), 1, int)
        return(res)
    }
    n <- length(x)
    a <- normorder(n)
    x <- sort(x)
    z <- (x - mean(x)) / sd(x)
    mod <- lm(z ~ 0 + a + I(a ^ 3))
    stat <- mod$coef[2] ^ 2
    return(stat = stat)
}

stat36 <- function(x) { # Xapd
    n <- length(x) 
    euler <- - digamma(1)
    z <- (x - mean(x)) / sqrt(var(x) * (n - 1) / n) 
    z <- z[! (z == 0)] 
    B2 <- sum(z ^ 2 * sign(z)) / n 
    K2 <- sum(z ^ 2 * log(abs(z))) / n 
    var.B2 <- (1 / n) * (3.0 - 8.0 / pi) * (1.0 - 1.9 / n) 
    Z.B2 <-  B2 / sqrt(var.B2)
    esp.net.K2 <- ((2.0 - log(2.0) - euler) / 2.0) ^ (1.0 / 3.0) *
        (1 - 1.026 / n) 
    var.net.K2 <- (1.0 / n) * (1.0 / 72.0) *
        ((2.0 - log(2.0) - euler) / 2.0) ^ (-4.0 / 3.0) *
        (3.0 * pi ^ 2 - 28.0) * (1.0 - 2.25 / n ^ 0.8) ;
    Z.net.K2 <- ((K2 - B2 ^ 2) ^ (1.0 / 3.0) - esp.net.K2) /
        sqrt(var.net.K2)
    Xapd.stat <- Z.B2 ^ 2 + Z.net.K2 ^ 2 ;
    p.value <- pchisq(Xapd.stat, 2, lower.tail = FALSE) ;
    res <- list(test = "2nd-power skewness and kurtosis omnibus normality test Xapd",
                B2 = B2, K2 = K2, Z.B2 = Z.B2, Z.net.K2 = Z.net.K2,
                Xapd.stat = Xapd.stat, p.value = p.value)
    return(res) 
}

stat37 <- function(x) { # Zepd
    n <- length(x) 
    euler <- -digamma(1)
    z <- (x - mean(x)) / sqrt(var(x) * (n - 1) / n) 
    z <- z[! (z == 0)] 
    K2 <- sum(z ^ 2 * log(abs(z))) / n 
    alpha.n <- -0.06 + 2.1 / n ^ 0.67
    K2.box.cox <- ((2.0 * K2) ^ alpha.n - 1.0) / alpha.n
    esp.K2.box.cox <- - ((2.0 - log(2.0) - euler) ^ (-0.06) - 1.0) / 0.06 -
        1.32 / n ^ 0.95
    var.K2.box.cox <- (1.0 / n) * ((2.0 - log(2.0) - euler) ^ (-2.12) *
                                   (3.0 * pi ^ 2 - 28.0) / 2.0 - 3.78 / n ^ 0.733)
    stat.Z.K2 <- (K2.box.cox - esp.K2.box.cox) / sqrt(var.K2.box.cox)
    p.value <- 2.0 * pnorm(abs(stat.Z.K2), lower.tail = FALSE)
    res <- list(test = "2nd-power kurtosis directional normality test Zepd",
                stat.Z.K2 = stat.Z.K2, p.value = p.value)
    return(res)
}
