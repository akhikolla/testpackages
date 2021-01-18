## -----------------------------------------------------------------------------
jackknife_error <- function (samples, na.rm = FALSE) {
  ## Number of jackknife samples.
  N <- length(samples)

  if (na.rm) {
    selection <- !is.na(samples)
    samples <- samples[selection]

    ## Number of non-NA samples.
    m <- sum(selection)
    factor <- N / m
  } else {
    factor <- 1.0
  }

  sqrt(factor * (N - 1) / N * sum((samples - mean(samples))^2))
}

## -----------------------------------------------------------------------------
N = 10
data = rnorm(N)

be = sd(data)
je = jackknife_error(data)
expected_factor = sqrt((N-1)^2 / N)
actual_factor = je / be

actual_factor / expected_factor

## -----------------------------------------------------------------------------
jackknife_cov <- function (x, y = NULL, na.rm = FALSE, ...) {
    factor <- 1.0
    
    if (is.null(y)) {
        N <- nrow(x)
        if (na.rm) {
            na_values <- apply(x, 2, function (row) any(is.na(row)))
            m <- sum(na_values)
            x <- x[!na_values, ]
            factor <- N / m
        }
    } else {
        N <- length(x)
        if (na.rm) {
        na_values <- is.na(x) | is.na(y)
            m <- sum(na_values)
            x <- x[!na_values]
            y <- y[!na_values]
            factor <- N / m
        }
    }
    
    (N-1)^2 / N * factor * cov(x, y, ...)
}

## -----------------------------------------------------------------------------
x <- rnorm(10)
y <- rnorm(10)

cov(x, y)
jackknife_cov(x, y)

cov(cbind(x, y))
jackknife_cov(cbind(x, y))

## -----------------------------------------------------------------------------
x <- rnorm(1000)

jackknife_samples_1 <- sapply(1:length(x), function (i) mean(x[-i]))
jackknife_samples_2 <- sapply(2:length(x), function (i) mean(x[-c(i-1, i)]))

jse1 <- jackknife_error(jackknife_samples_1)
jse2 <- jackknife_error(jackknife_samples_2)

jse2/jse1

