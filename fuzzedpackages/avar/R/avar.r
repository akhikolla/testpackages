#' @title Compute the Empirical Allan Variance
#'
#' @description
#' This function estimates the Allan variance.
#' @param x     A \code{vec} of time series observations or an \code{imu} object.
#' @param type  A \code{string} containing either \code{"mo"} for Maximal Overlap or \code{"to"} for Tau Overlap.
#' @param ...   Further arguments passed to other methods.
#' @return  If the input \code{x} is a \code{vec}, then the function returns a \code{list} that contains:
#' \itemize{
#'  \item "levels": The averaging time at each level.
#'  \item "allan": The estimated Allan variance.
#'  \item "type": Type of estimator (\code{mo} or \code{to}).
#' }
#' If the input \code{x} is an \code{imu} object, then the function returns a \code{list} that contains:
#' \itemize{
#'  \item "sensor": Name of the sensor.
#'  \item "freq": The frequency at which the error signal is measured.
#'  \item "n": Sample size of the data.
#'  \item "type": The types of sensors considered in the data.
#'  \item "axis": The axes of sensors considered in the data.
#'  \item "avar": A list containing the computed Allan variance based on the data.
#' }
#' @details
#' The decomposition and the amount of time it takes to perform this function depends on whether you are using
#' the Maximal Overlap or the Tau Overlap.
#'
#' @section Maximal Overlap Allan Variance:
#' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
#' we define \eqn{n} as an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
#' Therefore, \eqn{n} can be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n | n < floor(log2(N))}}
#' Based on the latter, we have \eqn{M = N - 2n} levels of decomposition.
#' The Maximal-overlap estimator is given by:
#' \deqn{\frac{1}{{2\left( {N - 2k + 1} \right)}}\sum\limits_{t = 2k}^N {{{\left[ {{{\bar Y}_t}\left( k \right) - {{\bar Y}_{t - k}}\left( k \right)} \right]}^2}} }{See PDF Manual}
#'
#' where \deqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }{See PDF Manual}.
#'
#' @section Tau-Overlap Allan Variance:
#' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
#' we define \eqn{n} as an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
#' Therefore, \eqn{n} can be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n | n < floor(log2(N))}}
#' Based on the latter, we have \eqn{m = \left\lfloor {\frac{{N - 1}}{n}} \right\rfloor  - 1} levels of decomposition.
#' The tau-overlap estimator is given by:
#'
#' where \eqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }{See PDF Manual}.
#'
#' @references Long-Memory Processes, the Allan Variance and Wavelets, D. B. Percival and P. Guttorp
#' @rdname avar
#' @export
#' @examples
#' set.seed(999)
#' Xt = rnorm(10000)
#' av_mat_mo = avar(Xt, type = "mo", freq = 100)
#' av_mat_tau = avar(Xt, type = "to")
#'
avar = function(x, type = "mo", ...){
  UseMethod("avar")
}

#' @rdname avar
#' @export
avar.default = function(x, type = "mo", freq = 1, ...) {

  if(is.null(x) | length(x) <=1 | dim(as.matrix(x))[2] >1){
    stop("Provide a vector or an 'imu' object")
  }

  if(sum(class(x) == "imu") == 1){
    freq = attributes(x)$freq
    x = as.vector(x)
  }

  if(type == "mo"){
    av = avar_mo_cpp(x)
  }else{
    av = avar_to_cpp(x)
  }

  av = list(levels = av[,1]/freq, allan=av[,2], errors = av[,3])
  av$adev = sqrt(av$allan)
  av$lci = av$adev - 2*av$errors*av$adev
  av$uci = av$adev + 2*av$errors*av$adev
  av$type = type
  av$freq = freq
  av$n = length(x)
  class(av) = c("avar", "list")
  av
}

#' @rdname avar
#' @export
avar.imu = function(x, type = "mo", ...){
  # Retrive sensor name
  if (!is.null(attr(x, "stype"))){
    sensor_name = attr(x, "stype")
  }else{
    warning("Unknown sensor name. IMU object is missing some information.")
    sensor_name = NULL
  }

  # Retrive freq
  if (!is.null(attr(x, "freq"))){
    freq = attr(x, "freq")
  }else{
    warning("Unknown frequency. IMU object is missing some information. Freq is set to 1 by default.")
    freq = 1
  }

  # Retrive sample size
  if (!is.null(attr(x, "dim"))){
    n = attr(x, "dim")[1]
  }else{
    warning("Unknown sample size. IMU object is missing some information.")
    n = NULL
  }

  # Retrive col names
  if (!is.null(attr(x, "dimnames")[[2]])){
    col_names = attr(x, "dimnames")[[2]]
  }else{
    stop("Unknown colunms names. IMU object is missing some information.")
    col_names = NULL
  }

  # Retrive sensor
  if (!is.null(attr(x, "sensor"))){
    sensor = attr(x, "sensor")
  }else{
    warning("Unknown sensor. IMU object is missing some information.")
    sensor = NULL
  }

  # Retrive axis
  if (!is.null(attr(x, "axis"))){
    ax = attr(x, "axis")
  }else{
    warning("Unknown axes. IMU object is missing some information.")
    ax = NULL
  }

  # Compute avar
  m = length(col_names)
  av = list()
  for (i in 1:m){
    av[[i]] = avar.default(x[,i], type = type, freq = freq)
  }
  names(av) = col_names
  out = list(sensor = sensor_name, freq = freq, n = n, type = sensor, axis = ax, avar = av)
  class(out) = "imu_avar"
  invisible(out)
}



#' @title Prints Allan Variance
#'
#' @description Displays the information on the output of the `avar()` function
#' @method print avar
#' @export
#' @param x   A \code{avar} object.
#' @param ... Arguments to be passed to methods
#' @return console output
#' @examples
#' set.seed(999)
#' Xt = rnorm(10000)
#' out = avar(Xt)
#' print(out)
#'
print.avar = function(x, ...) {
  cat("\n Levels: \n")
  print(x$levels, digits=5)
  cat("\n Allan Variances: \n")
  print(x$allan, digits=5)
  cat("\n Type: \n")
  print(x$type, digits=5)
}

#' @title Summary Allan Variance
#'
#' @description Displays the summary table of the output of the `avar()` function
#' @method summary avar
#' @export
#' @param object A \code{avar} object.
#' @param ...    Additional arguments affecting the summary produced.
#' A \code{table} that contains:
#' \itemize{
#'  \item "Time": The averaging time at each level.
#'  \item "AVar": The estimated Allan variance.
#'  \item "ADev": The estimated Allan deviation.
#'  \item "Lower CI": The lower bound of the confidence interval for the Allan deviation (ADev).
#'  \item "Upper CI": The upper bound of the confidence interval for the Allan deviation (ADev).
#' }
#' @examples
#' set.seed(999)
#' Xt = rnorm(10000)
#' out = avar(Xt)
#' summary(out)
#'
summary.avar = function(object, ...) {
  out_matrix = matrix(0, nrow = length(object$levels), ncol = 5)
  colnames(out_matrix) = c("Time", "AVar", "ADev", "Lower CI", "Upper CI")
  out_matrix[,"Time"] = object$levels
  out_matrix[,"AVar"] = object$allan
  out_matrix[,"ADev"] = object$adev
  out_matrix[,"Lower CI"] = object$lci
  out_matrix[,"Upper CI"] = object$uci

  class(out_matrix) = c("summary.avar","matrix")
  out_matrix
}

#' @title Plot Allan Deviation
#'
#' @description
#' Displays a plot of Allan deviation with its corresponding pointwise confidence intervals.
#' @method plot avar
#' @param x                An \code{avar} object.
#' @param units            A \code{string} that specifies the units of time plotted on the x axis.
#' @param xlab             A \code{string} that gives a title for the x axis.
#' @param ylab             A \code{string} that gives a title for the y axis.
#' @param main             A \code{string} that gives an overall title for the plot.
#' @param col_ad           A \code{string} that specifies the color of the line allan deviation line.
#' @param col_ci           A \code{string} that specifies the color of the shaded area covered by the confidence intervals.
#' @param ci_ad            A \code{boolean} that determines whether to plot the confidence interval shaded area.
#' @param nb_ticks_x       An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y       An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param legend_position  A \code{string} that specifies the position of the legend (use \code{legend_position = NA} to remove legend).
#' @param point_pch        A \code{double} that specifies the symbol type to be plotted.
#' @param point_cex        A \code{double} that specifies the size of each symbol to be plotted.
#' @param ...              Additional arguments affecting the plot.
#' @return A plot of the Allan deviation and relative confidence interval for each scale.
#' @author Stephane Guerrier, Nathanael Claussen and Justin Lee
#' @export
#' @examples
#' \donttest{
#' set.seed(999)
#' Xt = rnorm(10000)
#' av = avar(Xt)
#'
#' plot(av)
#' plot(av, main = "Simulated white noise", xlab = "Scales")
#' plot(av, units = "sec", legend_position = "topright")
#' plot(av, col_ad = "darkred", col_ci = "pink")
#'}
plot.avar = function(x, units = NULL, xlab = NULL, ylab = NULL, main = NULL,
                     col_ad = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                     legend_position = NULL, ci_ad = NULL, point_cex = NULL,
                     point_pch = NULL, ...){

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
  x_range = range(x$levels)
  if(length(x$levels) >= 10){
    x_low = floor(log10(x_range[1]))
    x_high = ceiling(log10(x_range[2]))
  }else{
    x_low = floor(log2(x_range[1]))
    x_high = ceiling(log2(x_range[2]))
  }

  y_range = range(cbind(x$adev - x$adev*x$errors, x$adev + x$adev*x$errors))
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

  if(length(x$levels) >= 10){
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

  par(new = TRUE)
  plot(NA, xlim = x_range, ylim = 10^c(win_dim[3], win_dim[4] + 0.45*(win_dim[4] - win_dim[3])),
       xlab = xlab, ylab = ylab, log = "xy", xaxt = 'n', yaxt = 'n', bty = "n")
  win_dim = par("usr")

  # Add Grid
  if(length(x$levels) >=10){
    abline(v = 10^x_ticks, lty = 1, col = "grey95")
  }else{
    abline(v = 2^x_ticks, lty = 1, col = "grey95")
  }

  abline(h = 10^y_ticks, lty = 1, col = "grey95")

  # Add Title
  x_vec = 10^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
  y_vec = 10^c(win_dim[4], win_dim[4],
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]),
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
  polygon(x_vec, y_vec, col = "grey95", border = NA)
  text(x = 10^mean(c(win_dim[1], win_dim[2])), y = 10^(win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), main)

  # Add Axes and Box
  lines(x_vec[1:2], rep(10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = 1)
  #y_ticks = y_ticks[(2^y_ticks) < 10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))]
  y_labels = y_labels[1:length(y_ticks)]
  box()
  if(length(x$levels) >=10){
    axis(1, at = 10^x_ticks, labels = x_labels, padj = 0.3)
  }else{
    axis(1, at = 2^x_ticks, labels = x_labels, padj = 0.3)
  }
  axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2)

  # CI for AD
  if (ci_ad == TRUE || is.null(ci_ad)){
    polygon(c(x$levels, rev(x$levels)), c(x$adev - x$errors*x$adev, rev(x$adev + x$errors*x$adev)),
            border = NA, col = col_ci)
  }

  # Add legend
  CI_conf = .95

  ad_title_part1 = "Empirical AD "

  if (!is.na(legend_position)){
    if (legend_position == "topleft"){
      legend_position = 10^c(1.1*win_dim[1], 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
      legend(x = legend_position[1], y = legend_position[2],
             legend = c(as.expression(bquote(paste(.(ad_title_part1), hat(phi)))),
                        as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")")))),
             pch = c(16, 15), lty = c(1, NA), col = c(col_ad, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
    }else{
      if (legend_position == "topright"){
        legend_position = 10^c(0.7*win_dim[2], 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
        legend(x = legend_position[1], y = legend_position[2],
               legend = c(as.expression(bquote(paste(.(ad_title_part1), hat(phi)))),
                          as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")")))),
               pch = c(16, 15), lty = c(1, NA), col = c(col_ad, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
      }else{
        legend(legend_position,
               legend = c(as.expression(bquote(paste(.(ad_title_part1), hat(phi)))),
                          as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")")))),
               pch = c(16, 15), lty = c(1, NA), col = c(col_ad, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
      }
    }
  }

  # Add AD
  if (is.null(point_pch)){
    point_pch = 16
  }

  if (is.null(point_cex)){
    point_cex = 1.25
  }
  lines(x$levels, x$adev, type = "l", col = col_ad, pch = 16)
  lines(x$levels, x$adev, type = "p", col = col_ad, pch = point_pch, cex = point_cex)
}

#' @title Plot Allan Deviation based on IMU Data
#'
#' @description
#' Displays a plot of Allan deviation based on IMU data with its corresponding pointwise confidence intervals.
#' @method plot imu_avar
#' @param x                An \code{avar} object.
#' @param xlab             A \code{string} that gives a title for the x axis.
#' @param ylab             A \code{string} that gives a title for the y axis.
#' @param main             A \code{string} that gives an overall title for the plot.
#' @param col_ad           A \code{string} that specifies the color of the line allan deviation line.
#' @param col_ci           A \code{string} that specifies the color of the shaded area covered by the confidence intervals.
#' @param ci_ad            A \code{boolean} that determines whether to plot the confidence interval shaded area.
#' @param nb_ticks_x       An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y       An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param point_pch        A \code{double} that specifies the symbol type to be plotted.
#' @param point_cex        A \code{double} that specifies the size of each symbol to be plotted.
#' @param ...              Additional arguments affecting the plot.
#' @return A plot of the Allan deviation and relative confidence interval for each scale.
#' @author Stephane Guerrier and Yuming Zhang
#' @export
#' @examples
#' data("navchip_av")
#' plot(navchip_av)
plot.imu_avar = function(x, xlab = NULL, ylab = NULL, main = NULL,
                         col_ad = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                         ci_ad = NULL, point_pch = NULL, point_cex = NULL, ...){
  type = unique(x$type)

  if ("Gyroscope" %in% type){
    gyro_index = which(x$type == "Gyroscope")
  }else{
    gyro_index = NULL
  }

  if ("Accelerometer" %in% type){
    accel_index = which(x$type == "Accelerometer")
  }else{
    accel_index = NULL
  }

  ncol = length(unique(x$axis))
  nrow = length(type)

  m = length(x$avar)
  J = length(x$avar[[1]]$allan)

  # remove negative CI values
  index_to_remove = c()
  for (i in 1:m) {
    if(length(which(x$avar[[i]]$lci<0)) > 0){
      index_to_remove = c(index_to_remove, which(x$avar[[i]]$lci<0))
    }
  }
  if (!is.null(index_to_remove)){
    index_to_remove = unique(index_to_remove)
    index_to_keep = which(seq(1:J) != index_to_remove)
  }else{
    index_to_keep = 1:J
  }

  J = length(index_to_keep)
  scales = x$avar[[1]]$levels[index_to_keep]

  ci_up = ci_lw = av = matrix(NA, J, m)

  for (i in 1:m){
    ci_up[,i] = x$avar[[i]]$uci[index_to_keep]
    ci_lw[,i] = x$avar[[i]]$lci[index_to_keep]
    av[,i] = x$avar[[i]]$allan[index_to_keep]
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
    main = paste("Allan Variance Representation - ", x$sensor, " @ ", x$freq, " Hz", sep="")
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
    y_range = c(min(ci_lw[,gyro_index]), max(ci_up[,gyro_index]))
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

      mtext(paste("Axis - ", x$axis[gyro_index][i], sep = ""), 3, line = 0.5)

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
    }
  }

  # Accel
  if (!is.null(accel_index)){
    y_range = c(min(ci_lw[,accel_index]), max(ci_up[,accel_index]))
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
        mtext(paste("Axis - ", x$axis[gyro_index][i], sep = ""), 3, line = 0.5)
      }

      abline(h = 10^y_ticks, col = "grey85")
      abline(v = 10^x_ticks, col = "grey85")

      # CI for AD
      if(ci_ad == TRUE || is.null(ci_ad)){
        polygon(c(scales, rev(scales)), c(ci_lw[,accel_index[i]], rev(ci_up[,accel_index[i]])),
                border = NA, col = col_ci)
      }

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

