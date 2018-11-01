# Copyright (C) 2017 - 2018 Stéphane Guerrier and Roberto Molinari
#
# This file is part of the `avar` R Methods Package
#
# The `avar` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#' Compute the Empirical Allan Variance
#'
#' This function estimates the Allan variance.
#' @param x     A \code{vec} of time series observations or an \code{imu} object.
#' @param type  A \code{string} containing either \code{"mo"} for Maximal Overlap or \code{"to"} for Tau Overlap.
#' @param freq  An \code{integer} with the frequency at which the error signal is measured.
#' @return  A \code{list} that contains:
#' \itemize{
#'  \item{"levels": }{The averaging time at each level.}
#'  \item{"allan": }{The estimated Allan variance.}
#'  \item{"type": }{Type of estimator (\code{mo} or \code{to}).}
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
#' @examples
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Simulate a Gaussian white noise
#' Xt = rnorm(10000)
#'
#' # Maximal overlap
#' av_mat_mo = avar(Xt, type = "mo", freq = 100)
#'
#' # Tau overlap
#' av_mat_tau = avar(Xt, type = "to")
#'
avar = function(x, type = "mo", freq = 1) {

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
  av$n = length(x)
  class(av) = c("avar", "list")
  av
}


#' Prints Allan Variance
#'
#' Displays the information on the output of the `avar()` function
#' @method print avar
#' @export
#' @param x   A \code{avar} object.
#' @param ... Arguments to be passed to methods
#' @return console output
#' @examples
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Simulate a Gaussian white noise
#' Xt = rnorm(10000)
#'
#' # Compute Allan
#' out = avar(Xt)
#'
#' # Print results
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

#' Summary Allan Variance
#'
#' Displays the summary table of the output of the `avar()` function
#' @method summary avar
#' @export
#' @param object A \code{avar} object.
#' @param ...    Additional arguments affecting the summary produced.
#' A \code{table} that contains:
#' \itemize{
#'  \item{"Time": }{The averaging time at each level.}
#'  \item{"AVar": }{The estimated Allan variance.}
#'  \item{"ADev": }{The estimated Allan deviation.}
#'  \item{"Lower CI": }{The lower bound of the confidence interval for the Allan deviation (ADev).}
#'  \item{"Upper CI": }{The upper bound of the confidence interval for the Allan deviation (ADev).}
#' }
#' @examples
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Simulate a Gaussian white noise
#' Xt = rnorm(10000)
#'
#' # Compute Allan
#' out = avar(Xt)
#'
#' # Summary
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
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Simulate a Gaussian white noise
#' Xt = rnorm(10000)
#'
#' # Compute Allan
#' av = avar(Xt)
#'
#' # Plot example
#' plot(av)
#' plot(av, main = "Simulated white noise", xlab = "Scales")
#' plot(av, units = "sec", legend_position = "topright")
#' plot(av, col_ad = "darkred", col_ci = "pink")
#'
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
  lines(x$levels, x$adev, type = "l", col = col_ad, pch = 16)

  if (is.null(point_pch)){
    point_pch = 16
  }

  if (is.null(point_cex)){
    point_cex = 1.25
  }
  lines(x$levels, x$adev, type = "p", col = col_ad, pch = point_pch, cex = point_cex)
}

#' @title Computes the Allan Variance Linear Regression estimator
#'
#' @description
#' Estimate the parameters of time series models based on the Allan Variance Linear Regression (AVLR) approach
#' @param x     A \code{vec} of time series observations or an \code{imu} object.
#' @param qn    A \code{vec} specifying on which scales the parameters of a Quantization Noise (QN) should be computed.
#' @param wn    A \code{vec} specifying on which scales the parameters of a White Noise (WN) should be computed.
#' @param rw    A \code{vec} specifying on which scales the parameters of a Random Wakk (RW) should be computed.
#' @param dr    A \code{vec} specifying on which scales the parameters of a Drift (DR) should be computed.
#' @param ci    A \code{boolean} to compute parameter confidence intervals.
#' @param B     A \code{double} for the number of bootstrap replicates to compute the parameter confidence intervals.
#' @param alpha A \code{double} defining the level of the confidence interval (1 - `alpha`).
#' @return avlr   A \code{list} that contains:
#' \itemize{
#'  \item{"estimates"}{The estimated value of the parameters.}
#'  \item{"implied_ad"}{The Allan deviation implied by the estimated parameters.}
#'  \item{"implied_ad_decomp"}{The Allan deviation implied by the estimated parameters for each individual model (if more than one is specified).}
#'  \item{"av"}{The \code{avar} object computed from the provided data.}
#' }
#' @importFrom stats dnorm
#' @examples
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Simulate a Gaussian WN and RW
#' N = 100000
#' Xt = rnorm(N) + cumsum(rnorm(N, 0, 3e-3))
#'
#' # Compute Allan variance
#' av = avar(Xt)
#' plot(av)
#'
#' # Parameter estimation
#' fit = avlr(Xt, wn = 1:8, rw = 11:15)
#' plot(fit, decomp = TRUE)
#'
#' # Point estimates
#' fit
#'
#' # Compute confidence intervals (this step is time-demanding)
#' # fit = avlr(x, wn = 1:8, rw = 10:15, ci = TRUE, B = 30)
#'
#' # Estimated confidence intervals and standard deviations
#' # fit$ci
avlr = function(x, qn = NULL, wn = NULL, rw = NULL, dr = NULL,
                ci = FALSE, B = 100, alpha = 0.05){

  if (ci == TRUE){
    stop("Method of confidence is not currently supported as it depends on an external package.")
  }

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

  n_processes = 4 - sum(sapply(list(qn,wn,rw,dr), is.null))

  process = rep(NA,n_processes)
  param = rep(NA,n_processes)
  implied = matrix(NA,length(x$levels),n_processes)

  counter = 0

  if(!is.null(wn)){
    if(length(wn) < 1 || !is.whole(wn) || min(wn) < 1 || max(wn) > length(x$allan)){
      stop("wn incorrectly formatted.")
    }
    counter = counter + 1
    process[counter] = "WN"
    param[counter] = exp(mean(log(x$adev[wn]) + log(x$levels[wn])/2))
    implied[,counter] = param[counter]/sqrt(x$levels)

    if (counter == 1){
      model_estimated = WN(sigma2 = (param[counter])^2)
    }else{
      model_estimated = model_estimated + WN(sigma2 = (param[counter])^2)
    }
  }

  if(!is.null(qn)){
    if(length(qn) < 1 || !is.whole(qn) || min(qn) < 1 || max(qn) > length(x$allan)){
      stop("qn incorrectely formatted.")
    }
    counter = counter + 1
    process[counter] = "QN"
    param[counter] = (1/sqrt(3))*exp(mean(log(x$adev[qn]) + log(x$levels[qn])))
    implied[,counter] = sqrt(3)*param[counter]/(x$levels)

    if (counter == 1){
      model_estimated = QN(q2 = (param[counter])^2)
    }else{
      model_estimated = model_estimated + QN(q2 = (param[counter])^2)
    }
  }

  if(!is.null(rw)){
    if(length(rw) < 1 || !is.whole(rw) || min(rw) < 1 || max(rw) > length(x$allan)){
      stop("rw incorrectely formatted.")
    }
    counter = counter + 1
    process[counter] = "RW"
    param[counter] = sqrt(3)*exp(mean(log(x$adev[rw]) - log(x$levels[rw])/2))
    implied[,counter] = param[counter]*sqrt(x$levels/3)

    if (counter == 1){
      model_estimated = RW(gamma2 = (param[counter])^2)
    }else{
      model_estimated = model_estimated + RW(gamma2 = (param[counter])^2)
    }
  }

  if(!is.null(dr)){
    if(length(dr) < 1 || !is.whole(dr) || min(dr) < 1 || max(dr) > length(x$allan)){
      stop("dr incorrectely formatted.")
    }
    counter = counter + 1
    process[counter] = "DR"
    param[counter] = sqrt(2)*exp(mean(log(x$adev[dr]) - log(x$levels[dr])))
    implied[,counter] = param[counter]*x$levels/2

    if (counter == 1){
      model_estimated = DR(omega = param[counter])
    }else{
      model_estimated = model_estimated + DR(omega = param[counter])
    }
  }

  implied_ad = apply(implied, 1, sum)

  estimates = t(t(param))
  rownames(estimates) = process
  colnames(estimates) = "Value"

  # Bootstrap parameters
  if (ci == TRUE){
    out_boot = boostrap_ci_avlr(model = model_estimated,
                                B = B, n = x$n, qn = qn,
                                wn = wn, rw = rw, dr = dr,
                                alpha = alpha)
    param = 2*param - out_boot$mu
    out_boot$ci = cbind(param - dnorm(1-alpha/2)*out_boot$sd, param + dnorm(1-alpha/2)*out_boot$sd)
    print("Parameter estimates corrected for bias via bootstrap")
  }else{
    out_boot = NULL
  }

  x = structure(list(estimates = param,
                     process_desc = process,
                     implied_ad = implied_ad,
                     implied_ad_decomp = implied,
                     av = x,
                     model = model_estimated,
                     ci = out_boot), class = "avlr")
  invisible(x)
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
#' \dontrun{
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Simulate time series
#' N = 100000
#' Xt = rnorm(N) + cumsum(rnorm(N, 0, 3e-3))
#'
#' # Maximal overlap
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
#' @param ci_ad            A \code{boolean} that determines whether to plot the confidence interval shaded area.
#' @param nb_ticks_x       An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y       An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param legend_position  A \code{string} that specifies the position of the legend (use \code{legend_position = NA} to remove legend).
#' @param point_pch        A \code{double} that specifies the symbol type to be plotted.
#' @param point_cex        A \code{double} that specifies the size of each symbol to be plotted.
#' @param ...              Additional arguments affecting the plot.
#' @return Plot of Allan deviation and relative confidence intervals for each scale.
#' @author Stephane Guerrier and Justin Lee
#' @export
#' @examples
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Simulate time series
#' N = 100000
#' Xt = rnorm(N) + cumsum(rnorm(N, 0, 3e-3))
#' av = avlr(Xt, wn = 1:7, rw = 12:15)
#'
#' # Plot example
#' plot.avlr(av)
#' plot.avlr(av, decomp = TRUE, main = "Simulated white noise", xlab = "Scales")
#' plot.avlr(av, units = "sec", legend_position = "topright")
#' plot.avlr(av, col_ad = "darkred", col_ci = "pink")
plot.avlr = function(x, decomp = FALSE,
                     units = NULL, xlab = NULL, ylab = NULL, main = NULL,
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
  x_range = range(x$av$levels)
  if(length(x$av$levels) >= 10){
    x_low = floor(log10(x_range[1]))
    x_high = ceiling(log10(x_range[2]))
  }else{
    x_low = floor(log2(x_range[1]))
    x_high = ceiling(log2(x_range[2]))
  }


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

  par(new = TRUE)
  plot(NA, xlim = x_range, ylim = 10^c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])),
       xlab = xlab, ylab = ylab, log = "xy", xaxt = 'n', yaxt = 'n', bty = "n")
  win_dim = par("usr")

  # Add Grid
  if(length(x$av$levels) >=10){
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
  if(length(x$av$levels) >=10){
    axis(1, at = 10^x_ticks, labels = x_labels, padj = 0.3)
  }else{
    axis(1, at = 2^x_ticks, labels = x_labels, padj = 0.3)
  }
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
  lines(t(x$av$levels),x$implied_ad, type = "l", lwd = 3, col = "#F47F24", pch = 1, cex = 1.5)
  lines(t(x$av$levels),x$implied_ad, type = "p", lwd = 2, col = "#F47F24", pch = 1, cex = 1.5)

  # Add AD
  lines(x$av$levels, x$av$adev, type = "l", col = col_ad, pch = 16)

  if (is.null(point_pch)){
    point_pch = 16
  }

  if (is.null(point_cex)){
    point_cex = 1.25
  }
  lines(x$av$levels, x$av$adev, type = "p", col = col_ad, pch = point_pch, cex = point_cex)

  # Add legend
  CI_conf = .95
  ad_title_part1 = "Empirical AD "

  if(decomp == TRUE){
    legend_names = c(as.expression(bquote(paste(.(ad_title_part1), hat(phi)))),
                     as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")"))),"Implied AD",
                     x$process_desc)
    col_legend = c(col_ad, col_ci,"#F47F24",col_decomp)
    p_cex_legend = c(1.25, 3, 1.5,rep(NA,U))
    lty_legend = c(1, NA, rep(1,U))
    pch_legend = c(16,15,1,rep(NA,U))
  }else{
    legend_names = c(as.expression(bquote(paste(.(ad_title_part1), hat(phi)))),
                     as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")"))),"Implied AV")
    col_legend = c(col_ad, col_ci,"#F47F24")
    p_cex_legend = c(1.25, 3, 1.5)
    lty_legend = c(1, NA)
    pch_legend = c(16,15,1)
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
#'  \item{"ci"}{The 1-\code{alpha} confidence intervals.}
#'  \item{"sd"}{The standard deviation of the estimated parameters.}
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

