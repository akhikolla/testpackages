
## -------------------------------------------------------------------
fz <- function (m, delta) {

    ## Argument
    idx <- 0:(m - 1)
    z <- complex(argument = -delta*c(idx, idx - m)^2)

    stats::fft(z)
}

fft_array <- function (mv_FFT_vals_array, FFT_vect_list) {
    list(mv_FFT_vals_array, FFT_vect_list)
}
