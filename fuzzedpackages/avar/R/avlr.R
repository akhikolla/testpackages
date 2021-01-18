#' @title Computes the Allan Variance Linear Regression estimator
#'
#' @description
#' Estimate the parameters of time series models based on the Allan Variance Linear Regression (AVLR) approach
#'
#' @param x     A \code{vec} of time series observations or an \code{imu} object.
#' @param qn    A \code{vec} specifying on which scales the parameters of a Quantization Noise (QN) should be computed.
#' @param wn    A \code{vec} specifying on which scales the parameters of a White Noise (WN) should be computed.
#' @param rw    A \code{vec} specifying on which scales the parameters of a Random Wakk (RW) should be computed.
#' @param dr    A \code{vec} specifying on which scales the parameters of a Drift (DR) should be computed.
#' @param qn_gyro    A \code{vec} specifying on which scales the parameters of a Quantization Noise (QN) should be computed for the gyroscope component.
#' @param wn_gyro    A \code{vec} specifying on which scales the parameters of a White Noise (WN) should be computed for the gyroscope component.
#' @param rw_gyro    A \code{vec} specifying on which scales the parameters of a Random Wakk (RW) should be computed for the gyroscope component.
#' @param dr_gyro    A \code{vec} specifying on which scales the parameters of a Drift (DR) should be computed for the gyroscope component.
#' @param qn_acc     A \code{vec} specifying on which scales the parameters of a Quantization Noise (QN) should be computed for the accelerometer component.
#' @param wn_acc     A \code{vec} specifying on which scales the parameters of a White Noise (WN) should be computed for the accelerometer component.
#' @param rw_acc     A \code{vec} specifying on which scales the parameters of a Random Wakk (RW) should be computed for the accelerometer component.
#' @param dr_acc     A \code{vec} specifying on which scales the parameters of a Drift (DR) should be computed for the accelerometer component.
#' @param ci    A \code{boolean} to compute parameter confidence intervals.
#' @param B     A \code{double} for the number of bootstrap replicates to compute the parameter confidence intervals.
#' @param alpha A \code{double} defining the level of the confidence interval (1 - `alpha`).
#' @param ...   Further arguments passed to other methods.
#'
#' @return If the input \code{x} is a \code{vec}, then the function returns a \code{list} that contains:
#' \itemize{
#'  \item "estimates": The estimated value of the parameters.
#'  \item "implied_ad": The Allan deviation implied by the estimated parameters.
#'  \item "implied_ad_decomp": The Allan deviation implied by the estimated parameters for each individual model (if more than one is specified).
#'  \item "av": The \code{avar} object computed from the provided data.
#' }
#' If the input \code{x} is of the class \code{imu_avar}, then the function returns a \code{list} that contains:
#' \itemize{
#'  \item "gyro": The estimation results correseponding to the gyroscope component.
#'  \item "acc": The estimation results correseponding to the accelerometer component.
#'  \item "imu_av": The \code{imu_avar} object computed based on the IMU data.
#' }
#' @importFrom stats dnorm
#' @rdname avlr
#' @examples
#' \donttest{
#' set.seed(999)
#'
#' N = 100000
#' Xt = rnorm(N) + cumsum(rnorm(N, 0, 3e-3))
#'
#' av = avar(Xt)
#' plot(av)
#'
#' # Input time series
#' fit = avlr(Xt, wn = 1:8, rw = 11:15)
#' fit
#'
#' # Input directly Allan variance
#' fit = avlr(av, wn = 1:8, rw = 11:15)
#' fit
#'
#' # Plot functions
#' plot(fit)
#' plot(fit, decomp = TRUE)
#' plot(fit, decomp = TRUE, show_scales = TRUE)
#' }
avlr = function(x, ...){
  UseMethod("avlr")
}


#' @rdname avlr
#' @export
avlr.default = function(x, qn = NULL, wn = NULL, rw = NULL, dr = NULL,
                        ci = FALSE, B = 100, alpha = 0.05, ...){

  if(is.null(x) | length(x) <=1){
    stop("Please provide a time series vector or an 'imu' object")
  }else if(class(x)[1] != "avar"){
    if(dim(as.matrix(x))[2] >1){
      stop("Please provide a time series vector or an 'imu' object")
    } else {
      x = avar(x)
    }
  }

  if(sum(sapply(list(qn,wn,rw,dr), is.null)) == 4){
    stop("Please specify a least one process.")
  }

  # Fit parameter
  fit = fit_avlr(qn = qn, wn = wn, rw = rw, dr = dr, ad = sqrt(x$allan), scales = x$levels)

  implied_ad = apply(fit$implied, 1, sum)
  estimates = t(t(fit$param))
  rownames(estimates) = fit$process
  colnames(estimates) = "Value"

  # Bootstrap parameters
  if (ci == TRUE){
    out_boot = boostrap_ci_avlr(model = fit$model_estimated,
                                B = B, n = x$n, qn = qn,
                                wn = wn, rw = rw, dr = dr,
                                alpha = alpha)
    fit$param = 2*fit$param - out_boot$mu
    out_boot$ci = cbind(fit$param - dnorm(1-alpha/2)*out_boot$sd, fit$param + dnorm(1-alpha/2)*out_boot$sd)
    print("Parameter estimates corrected for bias via bootstrap")
  }else{
    out_boot = NULL
  }

  # Used scales
  scales_used = matrix(NA, 4, 2)
  if (!is.null(qn))
    scales_used[1, ] = range(qn)

  if (!is.null(wn))
    scales_used[2, ] = range(wn)

  if (!is.null(rw))
    scales_used[3, ] = range(rw)

  if (!is.null(dr))
    scales_used[4, ] = range(dr)

  x = structure(list(estimates = fit$param,
                     process_desc = fit$process,
                     implied_ad = implied_ad,
                     implied_ad_decomp = fit$implied,
                     av = x,
                     model = fit$model_estimated,
                     ci = out_boot, scales_used = scales_used), class = "avlr")
  invisible(x)
}

#' @title Internal function to the Allan Variance Linear Regression estimator
#'
#' @description
#' Estimate the parameters of time series models based on the Allan Variance Linear Regression (AVLR) approach
#' @param qn     A \code{vec} specifying on which scales the parameters of a Quantization Noise (QN) should be computed.
#' @param wn     A \code{vec} specifying on which scales the parameters of a White Noise (WN) should be computed.
#' @param rw     A \code{vec} specifying on which scales the parameters of a Random Wakk (RW) should be computed.
#' @param dr     A \code{vec} specifying on which scales the parameters of a Drift (DR) should be computed.
#' @param ad     A \code{vec} of the Allan variance.
#' @param scales A \code{vec} of the scales.
#' @return       A \code{list} with the estimated parameters.
fit_avlr = function(qn, wn, rw, dr, ad, scales){
  # Number of processes needed
  n_processes = 4 - sum(sapply(list(qn,wn,rw,dr), is.null))

  # Initialisation
  process = rep(NA,n_processes)
  param = rep(NA,n_processes)
  implied = matrix(NA,length(scales),n_processes)
  counter = 0

  # Fit WN (if any)
  if(!is.null(wn)){
    if(length(wn) < 1 || !is.whole(wn) || min(wn) < 1 || max(wn) > length(ad)){
      stop("wn incorrectly formatted.")
    }
    counter = counter + 1
    process[counter] = "WN"
    param[counter] = exp(mean(log(ad[wn]) + log(scales[wn])/2))
    implied[,counter] = param[counter]/sqrt(scales)

    if (counter == 1){
      model_estimated = WN(sigma2 = (param[counter])^2)
    }else{
      model_estimated = model_estimated + WN(sigma2 = (param[counter])^2)
    }
  }

  # Fit QN (if any)
  if(!is.null(qn)){
    if(length(qn) < 1 || !is.whole(qn) || min(qn) < 1 || max(qn) > length(ad)){
      stop("qn incorrectely formatted.")
    }
    counter = counter + 1
    process[counter] = "QN"
    param[counter] = (1/sqrt(3))*exp(mean(log(ad[qn]) + log(scales[qn])))
    implied[,counter] = sqrt(3)*param[counter]/(scales)

    if (counter == 1){
      model_estimated = QN(q2 = (param[counter])^2)
    }else{
      model_estimated = model_estimated + QN(q2 = (param[counter])^2)
    }
  }

  # Fit RW (if any)
  if(!is.null(rw)){
    if(length(rw) < 1 || !is.whole(rw) || min(rw) < 1 || max(rw) > length(ad)){
      stop("rw incorrectely formatted.")
    }
    counter = counter + 1
    process[counter] = "RW"
    param[counter] = sqrt(3)*exp(mean(log(ad[rw]) - log(scales[rw])/2))
    implied[,counter] = param[counter]*sqrt(scales/3)

    if (counter == 1){
      model_estimated = RW(gamma2 = (param[counter])^2)
    }else{
      model_estimated = model_estimated + RW(gamma2 = (param[counter])^2)
    }
  }

  # Fit drift (if any)
  if(!is.null(dr)){
    if(length(dr) < 1 || !is.whole(dr) || min(dr) < 1 || max(dr) > length(ad)){
      stop("dr incorrectely formatted.")
    }
    counter = counter + 1
    process[counter] = "DR"
    param[counter] = sqrt(2)*exp(mean(log(ad[dr]) - log(scales[dr])))
    implied[,counter] = param[counter]*scales/2

    if (counter == 1){
      model_estimated = DR(omega = param[counter])
    }else{
      model_estimated = model_estimated + DR(omega = param[counter])
    }
  }

  return(list(implied = implied, model_estimated = model_estimated, param = param,
              process = process))
}


#' @rdname avlr
#' @export
avlr.imu_avar = function(x, qn_gyro = NULL, wn_gyro = NULL, rw_gyro = NULL, dr_gyro = NULL,
                         qn_acc = NULL, wn_acc = NULL, rw_acc = NULL, dr_acc = NULL,
                         B = 100, alpha = 0.05, ...){

  if(sum(sapply(list(qn_gyro,wn_gyro,rw_gyro,dr_gyro), is.null)) == 4 && "Gyroscope" %in% x$type){
    stop("Please specify a least one process (Gyro).")
  }

  if(sum(sapply(list(qn_acc,wn_acc,rw_acc,dr_acc), is.null)) == 4 && "Accelerometer" %in% x$type){
    stop("Please specify a least one process (Accel).")
  }

  # Retrive scales
  scales = x$avar[[1]]$levels
  J = length(scales)

  # Fit parameter of Gyro
  if ("Gyroscope" %in% x$type){
    # Number of axes
    m = sum(x$type %in% "Gyroscope")

    # Adjust selected scales
    if (!is.null(qn_gyro))
      qn_gyro = rep(qn_gyro,  m) + J*rep(0:(m-1), each = length(qn_gyro))

    if (!is.null(wn_gyro))
      wn_gyro = rep(wn_gyro,  m) + J*rep(0:(m-1), each = length(wn_gyro))

    if (!is.null(rw_gyro))
      rw_gyro = rep(rw_gyro,  m) + J*rep(0:(m-1), each = length(rw_gyro))

    if (!is.null(dr_gyro))
      dr_gyro = rep(dr_gyro,  m) + J*rep(0:(m-1), each = length(dr_gyro))

    # Index gyro
    ind_gyro = 1:(length(x$type)) * (x$type %in% "Gyroscope")
    ind_gyro = ind_gyro[ind_gyro > 0]
    ad_gyro = NULL
    for (i in 1:m){
      ad_gyro = c(ad_gyro, sqrt(x$avar[[ind_gyro[i]]]$allan))
    }
    scales_gyro = rep(scales, m)

    # Fit parameter
    fit_gyro = fit_avlr(qn = qn_gyro, wn = wn_gyro, rw = rw_gyro, dr = dr_gyro,
                        ad = ad_gyro, scales = scales_gyro)

    implied_ad_gyro = apply(fit_gyro$implied, 1, sum)
    estimates_gyro = t(t(fit_gyro$param))
    rownames(estimates_gyro) = fit_gyro$process
    colnames(estimates_gyro) = "Value"

    # Used scales
    scales_used_gyro = matrix(NA, 4, 2)
    if (!is.null(qn_gyro))
      scales_used_gyro[1, ] = c(qn_gyro[1], qn_gyro[length(qn_gyro)/m])

    if (!is.null(wn_gyro))
      scales_used_gyro[2, ] = c(wn_gyro[1], wn_gyro[length(wn_gyro)/m])

    if (!is.null(rw_gyro))
      scales_used_gyro[3, ] = c(rw_gyro[1], rw_gyro[length(rw_gyro)/m])

    if (!is.null(dr_gyro))
      scales_used_gyro[4, ] = c(dr_gyro[1], dr_gyro[length(dr_gyro)/m])

    gyro_out = list(estimates = fit_gyro$param,
                    process_desc = fit_gyro$process,
                    implied_ad = implied_ad_gyro,
                    implied_ad_decomp = fit_gyro$implied,
                    model = fit_gyro$model_estimated,
                    scales_used = scales_used_gyro)
  }else{
    gyro_out = NULL
  }


  # Fit parameter of Acce
  if ("Accelerometer" %in% x$type){
    # Number of axes
    m = sum(x$type %in% "Accelerometer")

    # Adjust selected scales
    if (!is.null(qn_acc))
      qn_acc = rep(qn_acc,  m) + J*rep(0:(m-1), each = length(qn_acc))

    if (!is.null(wn_acc))
      wn_acc = rep(wn_acc,  m) + J*rep(0:(m-1), each = length(wn_acc))

    if (!is.null(rw_acc))
      rw_acc = rep(rw_acc,  m) + J*rep(0:(m-1), each = length(rw_acc))

    if (!is.null(dr_acc))
      dr_acc = rep(dr_acc,  m) + J*rep(0:(m-1), each = length(dr_acc))

    # Find index accel
    ind_acc = 1:(length(x$type)) * (x$type %in% "Accelerometer")
    ind_acc = ind_acc[ind_acc > 0]
    ad_acc = NULL
    for (i in 1:m){
      ad_acc = c(ad_acc, sqrt(x$avar[[ind_acc[i]]]$allan))
    }
    scales_acc = rep(scales, m)

    # Fit parameter
    fit_acc = fit_avlr(qn = qn_acc, wn = wn_acc, rw = rw_acc, dr = dr_acc,
                       ad = ad_acc, scales = scales_acc)

    implied_ad_acc = apply(fit_acc$implied, 1, sum)
    estimates_acc = t(t(fit_acc$param))
    rownames(estimates_acc) = fit_acc$process
    colnames(estimates_acc) = "Value"

    # Used scales
    scales_used_acc = matrix(NA, 4, 2)
    if (!is.null(qn_acc))
      scales_used_acc[1, ] = c(qn_acc[1], qn_acc[length(qn_acc)/m])

    if (!is.null(wn_acc))
      scales_used_acc[2, ] = c(wn_acc[1], wn_acc[length(wn_acc)/m])

    if (!is.null(rw_acc))
      scales_used_acc[3, ] = c(rw_acc[1], rw_acc[length(rw_acc)/m])

    if (!is.null(dr_acc))
      scales_used_acc[4, ] = c(dr_acc[1], dr_acc[length(dr_acc)/m])

    acc_out = list(estimates = fit_acc$param,
                   process_desc = fit_acc$process,
                   implied_ad = implied_ad_acc,
                   implied_ad_decomp = fit_acc$implied,
                   model = fit_acc$model_estimated,
                   scales_used = scales_used_acc)
  }else{
    acc_out = NULL
  }

  out = list(gyro = gyro_out, acc = acc_out, imu_av = x)
  class(out) = "imu_avlr"
  invisible(out)
}



#' Print avlr object
#'
#' Displays information about the avlr object
#' @method print avlr
#' @export
#' @keywords internal
#' @param x   A \code{avlr} object
#' @param ... Other arguments passed to specific methods
#' @return Text output via print
#' @examples
#' \donttest{
#' set.seed(999)
#'
#' N = 100000
#' Xt = rnorm(N) + cumsum(rnorm(N, 0, 3e-3))
#'
#' fit = avlr(Xt, wn = 1:7, rw = 12:15)
#' print(fit)
#' }
print.avlr = function(x, ...) {
  if(is.null(x$ci)){
    cat("\n Estimates: \n")
    estimates = t(t(x$estimates))
    rownames(estimates) = x$process_desc
    colnames(estimates) = "Value"
    print(estimates)
  }else{
    cat("\n Estimates: \n")
    estimates = t(t(x$estimates))
    mat = cbind(t(t(x$estimates)),x$ci$ci, x$ci$sd )
    rownames(mat) = x$process_desc
    colnames(mat) = c("Value", "CI Low", "CI High", "SD")
    print(mat)
  }
}


#' Print imu_avlr object
#'
#' Displays information about the avlr object
#' @method print imu_avlr
#' @export
#' @keywords internal
#' @param x   A \code{avlr} object
#' @param ... Other arguments passed to specific methods
#' @return Text output via print
#' @examples
#' \donttest{
#' data(navchip_av)
#' navchip_avlr = avlr(navchip_av, wn_gyro = 1:20, rw_gyro = 1:20, wn_acc = 1:20, rw_acc = 1:20)
#' print(navchip_avlr)
#' }
print.imu_avlr = function(x, ...) {
  if("Gyroscope" %in% x$imu_av$type){
    cat("\n Estimates for gyroscopes: \n")
    estimates = t(t(x$gyro$estimates))
    rownames(estimates) = x$gyro$process_desc
    colnames(estimates) = "Value"
    print(estimates)
  }

  if("Accelerometer" %in% x$imu_av$type){
    cat("\n Estimates for accelerometers: \n")
    estimates = t(t(x$acc$estimates))
    rownames(estimates) = x$acc$process_desc
    colnames(estimates) = "Value"
    print(estimates)
  }
}


#' @title Plot the AVLR with the Allan Deviation for IMU
#'
#' @description
#' Displays a plot of the Allan deviation (AD) with the CI values and the AD implied by the estimated parameters for the IMU.
#' @method plot imu_avlr
#' @param x                An \code{avlr} object.
#' @param xlab             A \code{string} that gives a title for the x axis.
#' @param ylab             A \code{string} that gives a title for the y axis.
#' @param main             A \code{string} that gives an overall title for the plot.
#' @param col_ad           A \code{string} that specifies the color of the line allan deviation line.
#' @param col_ci           A \code{string} that specifies the color of the shaded area covered by the confidence intervals.
#' @param nb_ticks_x       An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y       An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param ci_ad            A \code{boolean} that determines whether to plot the confidence interval shaded area.
#' @param point_pch        A \code{double} that specifies the symbol type to be plotted.
#' @param point_cex        A \code{double} that specifies the size of each symbol to be plotted.
#' @param ...              Additional arguments affecting the plot.
#' @return Plot of Allan deviation and relative confidence intervals for each scale.
#' @author Stephane Guerrier and Justin Lee
#' @export
#' @examples
#' \donttest{
#' data(navchip_av)
#' navchip_avlr = avlr(navchip_av, wn_gyro = 1:20, rw_gyro = 1:20, wn_acc = 1:20, rw_acc = 1:20)
#' plot(navchip_avlr)
#' }
plot.imu_avlr = function(x, xlab = NULL, ylab = NULL, main = NULL,
                         col_ad = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                         ci_ad = NULL, point_pch = NULL, point_cex = NULL, ...){
  type = unique(x$imu_av$type)

  if ("Gyroscope" %in% type){
    gyro_index = which(x$imu_av$type == "Gyroscope")
  }else{
    gyro_index = NULL
  }

  if ("Accelerometer" %in% type){
    accel_index = which(x$imu_av$type == "Accelerometer")
  }else{
    accel_index = NULL
  }

  ncol = length(unique(x$imu_av$axis))
  nrow = length(type)

  m = length(x$imu_av$avar)
  J = length(x$imu_av$avar[[1]]$allan)

  # remove negative CI values
  index_to_remove = c()
  for (i in 1:m) {
    if(length(which(x$imu_av$avar[[i]]$lci<0)) > 0){
      index_to_remove = c(index_to_remove, which(x$imu_av$avar[[i]]$lci<0))
    }
  }
  index_to_remove = unique(index_to_remove)
  index_to_keep = which(seq(1:J) != index_to_remove)

  J = length(index_to_keep)
  scales = x$imu_av$avar[[1]]$levels[index_to_keep]

  ci_up = ci_lw = av = matrix(NA, J, m)

  for (i in 1:m){
    ci_up[,i] = x$imu_av$avar[[i]]$uci[index_to_keep]
    ci_lw[,i] = x$imu_av$avar[[i]]$lci[index_to_keep]
    av[,i] = x$imu_av$avar[[i]]$allan[index_to_keep]
  }

  # Axes
  if (is.null(nb_ticks_x)){
    nb_ticks_x = 6
  }

  if (is.null(nb_ticks_y)){
    nb_ticks_y = 5
  }

  # Range
  x_range = range(scales)
  x_low = floor(log10(x_range[1]))
  x_high = ceiling(log10(x_range[2]))

  x_ticks = seq(x_low, x_high, by = 1)
  if (length(x_ticks) > nb_ticks_x){
    x_ticks = x_low + ceiling((x_high - x_low)/(nb_ticks_x + 1))*(0:nb_ticks_x)
  }
  x_labels = sapply(x_ticks, function(i) as.expression(bquote(10^ .(i))))


  # Line and CI colors
  if (is.null(col_ad)){
    col_ad = "darkblue"
  }

  if (is.null(col_ci)){
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }

  if (is.null(point_pch)){
    point_pch = 16
  }

  if (is.null(point_cex)){
    point_cex = 1.25
  }

  # Main Title
  if (is.null(main)){
    main = paste("Allan Variance Representation - ", x$imu_av$sensor, " @ ", x$imu_av$freq, " Hz", sep="")
  }

  # Labels
  if (is.null(xlab)){
    xlab = bquote(paste("Averaging time ", tau, " [sec]", sep = " "))
  }

  if (is.null(ylab)){
    ylab = expression(paste("Allan Deviation ", phi, sep = ""))
  }


  # Main plot
  par(omi=rep(1.0, 4), mar=c(0,0,0,0), mfrow=c(nrow,ncol))

  # Gyro
  if (!is.null(gyro_index)){
    y_range = c(min(ci_lw[,gyro_index]), max(c(ci_up[,gyro_index], x$gyro$implied_ad[1:length(scales)])))
    y_low = floor(log10(y_range[1]))
    y_high = ceiling(log10(y_range[2]))

    y_ticks <- seq(y_low, y_high, by = 1)
    if (length(y_ticks) > nb_ticks_y){
      y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
    }
    y_labels <- sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))


    for (i in seq_along(gyro_index)){
      plot(NA, xlim = range(scales), ylim = y_range, xaxt="n", yaxt="n", log = "xy", bty = "n")
      box(col = "grey")

      mtext(paste("Axis - ", x$imu_av$axis[gyro_index][i], sep = ""), 3, line = 0.5)

      if (i == 1){
        axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2, cex = 1.25)
      }

      if (i == 1){
        mtext("Gyroscope", 2, line = 4.5)
        mtext(ylab, 2, line = 2.5)
      }

      abline(h = 10^y_ticks, col = "grey85")
      abline(v = 10^x_ticks, col = "grey85")

      # CI for AD
      if(ci_ad == TRUE || is.null(ci_ad)){
        polygon(c(scales, rev(scales)), c(ci_lw[,gyro_index[i]], rev(ci_up[,gyro_index[i]])),
                border = NA, col = col_ci)
      }

      # Add AD
      lines(scales, sqrt(av[,gyro_index[i]]), type = "l", col = col_ad, pch = 16)
      lines(scales, sqrt(av[,gyro_index[i]]), type = "p", col = col_ad, pch = point_pch, cex = point_cex)

      if (is.null(accel_index)){
        axis(1, at = 10^x_ticks, labels = x_labels, padj = -0.2, cex = 1.25)
      }

      lines(scales,x$gyro$implied_ad[1:length(scales)], type = "b", lwd = 1.75, col = "#F47F24", pch = 1, cex = 1.5)
    }
  }

  # Accel
  if (!is.null(accel_index)){
    y_range = c(min(ci_lw[,accel_index]), max(c(ci_up[,accel_index], x$acc$implied_ad[1:length(scales)])))
    y_low = floor(log10(y_range[1]))
    y_high = ceiling(log10(y_range[2]))

    y_ticks <- seq(y_low, y_high, by = 1)
    if (length(y_ticks) > nb_ticks_y){
      y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
    }
    y_labels <- sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))


    for (i in seq_along(accel_index)){
      plot(NA, xlim = range(scales), ylim = y_range, xaxt="n", yaxt="n", log = "xy", bty = "n")
      box(col = "grey")
      if (i == 1){
        axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2, cex = 1.25)
      }


      if (i == 1){
        mtext("Accelerometer", 2, line = 4.5)
        mtext(ylab, 2, line = 2.5)
      }

      if (length(accel_index) == 3 && i == 2){
        mtext(xlab, 1, line = 3.5)
      }

      if (is.null(gyro_index)){
        mtext(paste("Axis - ", x$imu_av$axis[gyro_index][i], sep = ""), 3, line = 0.5)
      }

      abline(h = 10^y_ticks, col = "grey85")
      abline(v = 10^x_ticks, col = "grey85")

      # CI for AD
      if(ci_ad == TRUE || is.null(ci_ad)){
        polygon(c(scales, rev(scales)), c(ci_lw[,accel_index[i]], rev(ci_up[,accel_index[i]])),
                border = NA, col = col_ci)
      }

      lines(scales,x$acc$implied_ad[1:length(scales)], type = "b", lwd = 1.75, col = "#F47F24", pch = 1, cex = 1.5)

      # Add AD
      lines(scales, sqrt(av[,accel_index[i]]), type = "l", col = col_ad, pch = 16)
      lines(scales, sqrt(av[,accel_index[i]]), type = "p", col = col_ad, pch = point_pch, cex = point_cex)
      axis(1, at = 10^x_ticks, labels = x_labels, padj = -0.2, cex = 1.25)
    }
  }

  # Add main title
  mtext(main, side = 3, line = 3, outer = TRUE)
  par(mfrow = c(1,1))
}


#' @title Plot the AVLR with the Allan Deviation
#'
#' @description
#' Displays a plot of the Allan deviation (AD) with the CI values and the AD implied by the estimated parameters.
#' @method plot avlr
#' @param x                An \code{avlr} object.
#' @param decomp           A \code{boolean} that determines whether the contributions of each individual model are plotted.
#' @param units            A \code{string} that specifies the units of time plotted on the x axis.
#' @param xlab             A \code{string} that gives a title for the x axis.
#' @param ylab             A \code{string} that gives a title for the y axis.
#' @param main             A \code{string} that gives an overall title for the plot.
#' @param col_ad           A \code{string} that specifies the color of the line allan deviation line.
#' @param col_ci           A \code{string} that specifies the color of the shaded area covered by the confidence intervals.
#' @param nb_ticks_x       An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y       An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param legend_position  A \code{string} that specifies the position of the legend (use \code{legend_position = NA} to remove legend).
#' @param ci_ad            A \code{boolean} that determines whether to plot the confidence interval shaded area.
#' @param point_cex        A \code{double} that specifies the size of each symbol to be plotted.
#' @param point_pch        A \code{double} that specifies the symbol type to be plotted.
#' @param show_scales      A \code{boolean} that specifies if the scales used for each process should be plotted.
#' @param ...              Additional arguments affecting the plot.
#' @return Plot of Allan deviation and relative confidence intervals for each scale.
#' @author Stephane Guerrier and Justin Lee
#' @export
#' @examples
#' \donttest{
#' set.seed(999)
#'
#' N = 100000
#' Xt = rnorm(N) + cumsum(rnorm(N, 0, 3e-3))
#' av = avlr(Xt, wn = 1:7, rw = 12:15)
#'
#' plot.avlr(av)
#' plot.avlr(av, decomp = TRUE, main = "Simulated white noise", xlab = "Scales")
#' plot.avlr(av, units = "sec", legend_position = "topright")
#' plot.avlr(av, col_ad = "darkred", col_ci = "pink")
#' plot(fit, decomp = TRUE, show_scales = TRUE)
#' }
plot.avlr = function(x, decomp = FALSE,
                     units = NULL, xlab = NULL, ylab = NULL, main = NULL,
                     col_ad = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                     legend_position = NULL, ci_ad = NULL, point_cex = NULL,
                     point_pch = NULL, show_scales = FALSE, ...){


  # Labels
  if (is.null(xlab)){
    if (is.null(units)){
      xlab = expression(paste("Scale ", tau, sep =""))
    }else{
      xlab = bquote(paste("Averaging time ", tau, " [", .(units), "]", sep = " "))
    }
  }

  if (is.null(ylab)){
    ylab = expression(paste("Allan Deviation ", phi, sep = ""))
  }else{
    ylab = ylab
  }

  # Main Title
  if (is.null(main)){
    main = "Allan Deviation Representation"
  }

  # Line and CI colors
  if (is.null(col_ad)){
    col_ad = "darkblue"
  }

  if (is.null(col_ci)){
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }

  # Range
  x_range = range(x$av$levels)
  if(length(x$av$levels) >= 10){
    x_low = floor(log10(x_range[1]))
    x_high = ceiling(log10(x_range[2]))
  }else{
    x_low = floor(log2(x_range[1]))
    x_high = ceiling(log2(x_range[2]))
  }

  #compute y range
  y_range = range(cbind(x$av$adev - x$av$adev*x$av$errors, x$av$adev + x$av$adev*x$av$errors))
  y_low = floor(log10(y_range[1]))
  y_high = ceiling(log10(y_range[2]))

  # Axes
  if (is.null(nb_ticks_x)){
    nb_ticks_x = 6
  }

  if (is.null(nb_ticks_y)){
    nb_ticks_y = 5
  }

  x_ticks = seq(x_low, x_high, by = 1)
  if (length(x_ticks) > nb_ticks_x){
    x_ticks = x_low + ceiling((x_high - x_low)/(nb_ticks_x + 1))*(0:nb_ticks_x)
  }

  if(length(x$av$clusters) >= 10){
    x_labels = sapply(x_ticks, function(i) as.expression(bquote(10^ .(i))))
  }else{
    x_labels = sapply(x_ticks, function(i) as.expression(bquote(2^ .(i))))
  }

  y_ticks <- seq(y_low, y_high, by = 1)
  if (length(y_ticks) > nb_ticks_y){
    y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
  }
  y_labels <- sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))

  # Legend Position
  if (is.null(legend_position)){
    #if (which.min(abs(c(y_low, y_high) - log2(x$variance[1]))) == 1){
    #  legend_position = "topleft"
    #}else{
    legend_position = "bottomleft"
    #}
  }

  # Main Plot
  plot(NA, xlim = x_range, ylim = y_range, xlab = xlab, ylab = ylab,
       log = "xy", xaxt = 'n', yaxt = 'n', bty = "n", ann = FALSE)
  win_dim = par("usr")

  #par(new = TRUE)
  #plot(NA, xlim = x_range, ylim = 10^c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])), xlab = xlab, ylab = ylab, log = "xy", xaxt = 'n', yaxt = 'n', bty = "n")
  #win_dim = par("usr")

  # Add Grid
  if(length(x$av$levels) >=10){
    abline(v = 10^x_ticks, lty = 1, col = "grey95")
  }else{
    abline(v = 2^x_ticks, lty = 1, col = "grey95")
  }
  abline(h = 10^y_ticks, lty = 1, col = "grey95")

  #y_ticks = y_ticks[(2^y_ticks) < 10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))]
  y_labels = y_labels[1:length(y_ticks)]

  #x axis
  if(length(x$av$levels) >=10){
    axis(1, at = 10^x_ticks, labels = x_labels, padj = 0.3)
  }else{
    axis(1, at = 2^x_ticks, labels = x_labels, padj = 0.3)
  }

  #yaxis
  axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2)

  # CI for AD
  if (ci_ad == TRUE || is.null(ci_ad)){
    polygon(c(x$av$levels, rev(x$av$levels)), c(x$av$adev - x$av$errors*x$av$adev, rev(x$av$adev + x$av$errors*x$av$adev)),
            border = NA, col = col_ci)
  }

  U = dim(x$implied_ad_decomp)[2]
  col_decomp = hcl(h = seq(100, 375, length = U + 1), l = 65, c = 200, alpha = 1)[1:U]

  # Legend Position
  if (is.null(legend_position)){
    #if (which.min(abs(c(y_low, y_high) - log2(x$variance[1]))) == 1){
    #  legend_position = "topleft"
    #}else{
    legend_position = "bottomleft"
    #}
  }

  if(decomp == TRUE){
    # Plot lines of decomp theo
    for (i in 1:U){
      lines(x$av$levels, x$implied_ad_decomp[,i], col = col_decomp[i])
    }
  }
  # Plot implied AD
  lines(t(x$av$levels),x$implied_ad, type = "b", lwd = 2, col = "#F47F24", pch = 1, cex = 1.5)
  #lines(t(x$av$levels),x$implied_ad, type = "p", lwd = 2, col = "#F47F24", pch = 1, cex = 1.5)

  # Add AD
  lines(x$av$levels, x$av$adev, type = "l", col = col_ad, pch = 16)

  if (is.null(point_pch)){
    point_pch = 16
  }

  if (is.null(point_cex)){
    point_cex = 1.25
  }
  lines(x$av$levels, x$av$adev, type = "p", col = col_ad, pch = point_pch, cex = point_cex)

  if (show_scales){
    process_cols = col_decomp
    counter = 0
    for (i in 1:4){
      if (!is.na(x$scales_used[i,1])){
        counter = counter + 1
        lines(x$av$levels[(x$scales_used[i,1]):(x$scales_used[i,2])],
              x$av$adev[(x$scales_used[i,1]):(x$scales_used[i,2])], type = "p",
              col = process_cols[counter], pch = point_pch, cex = point_cex)
      }
    }
  }

  # Add Title
  x_vec = 10^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
  y_vec = 10^c(win_dim[4], win_dim[4],
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]),
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
  polygon(x_vec, y_vec, col = "grey95", border = NA)
  text(x = 10^mean(c(win_dim[1], win_dim[2])), y = 10^(win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), main)

  # Add line below title
  lines(x_vec[1:2], rep(10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = 1)

  #add box around plot
  box()

  # Add legend
  CI_conf = .95
  ad_title_part1 = "Empirical AD "

  if(decomp == TRUE){
    if (show_scales){
      legend_names = c(as.expression(bquote(paste(.(ad_title_part1), hat(phi)))),
                       as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")"))),"Implied AD",
                       x$process_desc, paste("Scales used for", x$process_desc))
      col_legend = c(col_ad, col_ci,"#F47F24",col_decomp, col_decomp)
      p_cex_legend = c(1.25, 3, 1.5,rep(NA,U), rep(1.35,U))
      lty_legend = c(1, NA, 1, rep(1,U), rep(NA,U))
      pch_legend = c(16,15, 1, rep(NA,U), rep(16,U))
    }else{
      legend_names = c(as.expression(bquote(paste(.(ad_title_part1), hat(phi)))),
                       as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")"))),"Implied AD",
                       x$process_desc)
      col_legend = c(col_ad, col_ci,"#F47F24",col_decomp)
      p_cex_legend = c(1.25, 3, 1.5,rep(NA,U))
      lty_legend = c(1, NA, rep(1,U))
      pch_legend = c(16,15,1,rep(NA,U))
    }

  }else{
    if (show_scales){
      legend_names = c(as.expression(bquote(paste(.(ad_title_part1), hat(phi)))),
                       as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")"))),"Implied AD",
                       paste("Scales used for", x$process_desc))
      col_legend = c(col_ad, col_ci,"#F47F24", col_decomp)
      p_cex_legend = c(1.25, 3, 1.5, rep(1.35,U))
      lty_legend = c(1, NA, 1, rep(NA,U))
      pch_legend = c(16,15, 1, rep(16,U))
    }else{
      legend_names = c(as.expression(bquote(paste(.(ad_title_part1), hat(phi)))),
                       as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")"))),"Implied AV")
      col_legend = c(col_ad, col_ci,"#F47F24")
      p_cex_legend = c(1.25, 3, 1.5)
      lty_legend = c(1, NA)
      pch_legend = c(16,15,1)
    }
  }
  if (!is.na(legend_position)){
    if (legend_position == "topleft"){
      legend_position = 10^c(1.1*win_dim[1], 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
      legend(x = legend_position[1], y = legend_position[2],
             legend = legend_names, pch = pch_legend, lty = lty_legend,
             col = col_legend, cex = 1, pt.cex = p_cex_legend, bty = "n")
    }else{
      if (legend_position == "topright"){
        legend_position = 10^c(0.7*win_dim[2], 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
        legend(x = legend_position[1], y = legend_position[2],
               legend =legend_names, pch = pch_legend, lty = lty_legend,
               col = col_legend, cex = 1, pt.cex = p_cex_legend, bty = "n")
      }else{
        legend(legend_position,
               legend = legend_names, pch = pch_legend, lty = lty_legend,
               col = col_legend, cex = 1, pt.cex = p_cex_legend, bty = "n")
      }
    }
  }
}

#' Compute bootstrap confidence intervals for the AVLR estimator
#'
#' @keywords internal
#' @importFrom stats quantile sd
#' @param model A \code{ts.model} object that was estimated with the avlr function.
#' @param B     A \code{double} for the number of bootsrap replicates to compute the confidence intervals.
#' @param n     A \code{double} with the sample size.
#' @param qn    A \code{vec} specifying on which scales the parameters of a Quantization Noise (QN) was computed.
#' @param wn    A \code{vec} specifying on which scales the parameters of a White Noise (WN) was computed.
#' @param rw    A \code{vec} specifying on which scales the parameters of a Random Wakk (RW) was computed.
#' @param dr    A \code{vec} specifying on which scales the parameters of a Drift (DR) was computed.
#' @param alpha A \code{double} defining the level of the confidence interval (1 - `alpha`).
#' @return   A \code{list} that contains:
#' \itemize{
#'  \item "ci": The 1-\code{alpha} confidence intervals.
#'  \item "sd": The standard deviation of the estimated parameters.
#' }
boostrap_ci_avlr = function(model, B, n, qn, wn, rw, dr, alpha){
  results = matrix(NA, B, model$plength)
  print("Starting bootstrap:")

  for (i in 1:B){
    x_star = gen_gts(n = n, model = model)
    results[i, ] = as.numeric(avlr(x_star, qn = qn, wn = wn, rw = rw,
                                   dr = dr, ci = FALSE)$estimates)
  }

  mean_parameters = rep(NA, model$plength)
  ci_parameters = matrix(NA, model$plength, 2)
  sd_parameters = rep(NA, model$plength)

  for (i in 1:model$plength){
    mean_parameters[i] = mean(results[,i])
    ci_parameters[i, ] = as.numeric(quantile(results[,i], probs = c(alpha/2, 1 - alpha/2)))
    sd_parameters[i] = sd(results[,i])
  }
  list(mu = mean_parameters, ci = ci_parameters, sd = sd_parameters)
}

