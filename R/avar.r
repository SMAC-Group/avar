# Copyright (C) 2018 - 2018 Stephane Guerrier, Gaetan Bakalli, James Balamuta
#
# This file is part of av R Methods Package
#
# The `av` R package is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# The `av` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Calculate the Allan Variance
#'
#' Computes the Allan Variance
#' @param x     A \code{vec} of time series observations.
#' @param type  A \code{string} containing either \code{"mo"} for Maximal Overlap or \code{"to"} for Tau Overlap
#' @return av   A \code{list} that contains:
#' \itemize{
#'  \item{"clusters"}{The size of the cluster}
#'  \item{"allan"}{The Allan variance}
#'  \item{"errors"}{The error associated with the variance estimation.}
#' }
#' @details
#' The decomposition and the amount of time it takes to perform it depends on whether you are using
#' the Tau Overlap or the Maximal Overlap.
#'
#' @section Maximal Overlap Allan Variance:
#' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
#' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
#' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n | n < floor(log2(N))}}
#' Then, \eqn{M = N - 2n} samples exist.
#' The Maximal-overlap estimator is given by:
#' \deqn{\frac{1}{{2\left( {N - 2k + 1} \right)}}\sum\limits_{t = 2k}^N {{{\left[ {{{\bar Y}_t}\left( k \right) - {{\bar Y}_{t - k}}\left( k \right)} \right]}^2}} }{See PDF Manual}
#'
#' where \deqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }{See PDF Manual}.
#'
#' @section Tau-Overlap Allan Variance:
#' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
#' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
#' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n | n < floor(log2(N))}}
#' Then, a sampling of \eqn{m = \left\lfloor {\frac{{N - 1}}{n}} \right\rfloor  - 1} samples exist.
#' The tau-overlap estimator is given by:
#'
#' where \eqn{ {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} }{See PDF Manual}.
#'
#' @references Long-Memory Processes, the Allan Variance and Wavelets, D. B. Percival and P. Guttorp
#' @examples
#' # Load simts package
#' library(simts)
#'
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Simulate time series
#' N = 100000
#' ts = gen_gts(N, WN(sigma2 = 2) + RW(gamma2 = 1))
#'
#' # Maximal overlap
#' av_mat_mo = avar(ts, type = "mo")
#'
#' # Tau overlap
#' av_mat_tau = avar(ts, type = "to")
avar = function(x, type = "mo") {
  x = as.vector(x)

  if(type == "mo"){
    av = avar_mo_cpp(x)
  }else{
    av = avar_to_cpp(x)
  }

  av = list(clusters = av[,1], allan=av[,2], errors=av[,3])
  av$adev = sqrt(av$allan)
  av$lci = av$adev - 2*av$errors*av$adev
  av$uci = av$adev + 2*av$errors*av$adev
  av$type = type
  class(av) = c("avar", "list")
  av
}


#' Prints Allan Variance
#'
#' Displays the allan variance information
#' @method print avar
#' @export
#' @param x   A \code{avar} object.
#' @param ... Arguments to be passed to methods
#' @return console output
#' @examples
#' # Load simts package
#' library(simts)
#'
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Generate time series
#' x = gen_gts(100, WN(sigma2 = 1))
#'
#' # Compute Allan
#' out = avar(x)
#'
#' # Print results
#' print( out )
print.avar = function(x, ...) {
  cat("\n Clusters: \n")
  print(x$clusters, digits=5)
  cat("\n Allan Variances: \n")
  print(x$allan, digits=5)
  cat("\n Errors: \n")
  print(x$errors, digits=5)
}

#' Summary Allan Variance
#'
#' Displays the summary table of allan variance
#' @method summary avar
#' @export
#' @param object A \code{avar} object.
#' @param ...    Additional arguments affecting the summary produced.
#' @return Summary table
#' @examples
#' # Load simts package
#' library(simts)
#'
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Generate time series
#' x = gen_gts(100, WN(sigma2 = 1))
#'
#' # Compute Allan
#' out = avar(x)
#'
#' # Summary
#' summary( out )
summary.avar = function(object, ...) {
  out_matrix = matrix(0, nrow = length(object$clusters), ncol = 6)
  colnames(out_matrix) = c("Time", "AVAR", "ADEV", "Lower CI", "Upper CI", "Error")
  out_matrix[,"Time"] = object$clusters
  out_matrix[,"AVAR"] = object$allan
  out_matrix[,"ADEV"] = object$adev
  out_matrix[,"Lower CI"] = object$lci
  out_matrix[,"Upper CI"] = object$uci
  out_matrix[,"Error"] = object$errors

  class(out_matrix) = c("summary.avar","matrix")
  out_matrix
}

#' @title Plot Allan Deviation
#'
#' @description
#' Displays a plot of allan deviation accounting for CI values.
#' @method plot avar
#' @param x                A \code{avar} object.
#' @param units            A \code{string} that specifies the units of time plotted on the x axis.
#' @param xlab             A \code{string} that gives a title for the x axis.
#' @param ylab             A \code{string} that gives a title for the y axis.
#' @param main             A \code{string} that gives an overall title for the plot.
#' @param col_wv           A \code{string} that specifies the color of the wavelet variance line.
#' @param col_ci           A \code{string} that specifies the color of the confidence interval polygon.
#' @param ci_wv            A \code{boolean} that determines whether a confidence interval polygon will be drawn.
#' @param nb_ticks_x       An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y       An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param legend_position  A \code{string} that specifies the position of the legend (use \code{legend_position = NA} to remove legend).
#' @param point_pch        A \code{double} that specifies the symbol type to be plotted.
#' @param point_cex        A \code{double} that specifies the size of each symbol to be plotted.
#' @param ...              Additional arguments affecting the plot.
#' @return Plot of allan deviation and confidence interval for each scale.
#' @author Stephane Guerrier, Nathanael Claussen, and Justin Lee
#' @export
#' @examples
#' # Load simts package
#' library(simts)
#'
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Generate time series
#' x = gen_gts(100, WN(sigma2 = 1))
#'
#' # Compute Allan
#' av = avar(x)
#'
#' # Plot example
#' plot(av)
#' plot(av, main = "Simulated white noise", xlab = "Scales")
#' plot(av, units = "sec", legend_position = "topright")
#' plot(av, col_wv = "darkred", col_ci = "pink")
plot.avar = function(x, units = NULL, xlab = NULL, ylab = NULL, main = NULL,
                     col_wv = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                     legend_position = NULL, ci_wv = NULL, point_cex = NULL,
                     point_pch = NULL, ...){

  # Labels
  if (is.null(xlab)){
    if (is.null(units)){
      xlab = expression(paste("Scale ", tau, sep =""))
    }else{
      xlab = bquote(paste("Clustering time ", tau, " [", .(units), "]", sep = " "))
    }
  }

  if (is.null(ylab)){
    ylab = expression(paste("Allan Deviation ", phi, sep = ""))
  }else{
    ylab = ylab
  }

  # Main Title
  if (is.null(main)){
    main = "Allan Variance Representation"
  }

  # Line and CI colors
  if (is.null(col_wv)){
    col_wv = "darkblue"
  }

  if (is.null(col_ci)){
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }

  # Range
  x_range = range(x$clusters)
  x_low = floor(log2(x_range[1]))
  x_high = ceiling(log2(x_range[2]))

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
  x_labels = sapply(x_ticks, function(i) as.expression(bquote(2^ .(i))))

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
  abline(v = 2^x_ticks, lty = 1, col = "grey95")
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
  axis(1, at = 2^x_ticks, labels = x_labels, padj = 0.3)
  axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2)

  # CI for WV
  if (ci_wv == TRUE || is.null(ci_wv)){
    polygon(c(x$cluster, rev(x$cluster)), c(x$adev - x$errors*x$adev, rev(x$adev + x$errors*x$adev)),
            border = NA, col = col_ci)
  }

  # Add legend
  CI_conf = .95

  wv_title_part1 = "Empirical AD "

  if(add_legend == TRUE){
    if (!is.na(legend_position)){
      if (legend_position == "topleft"){
        legend_position = 10^c(1.1*win_dim[1], 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
        legend(x = legend_position[1], y = legend_position[2],
               legend = c(as.expression(bquote(paste(.(wv_title_part1), hat(phi)))),
                          as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")")))),
               pch = c(16, 15), lty = c(1, NA), col = c(col_wv, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
      }else{
        if (legend_position == "topright"){
          legend_position = 10^c(0.7*win_dim[2], 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
          legend(x = legend_position[1], y = legend_position[2],
                 legend = c(as.expression(bquote(paste(.(wv_title_part1), hat(phi)))),
                            as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")")))),
                 pch = c(16, 15), lty = c(1, NA), col = c(col_wv, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
        }else{
          legend(legend_position,
                 legend = c(as.expression(bquote(paste(.(wv_title_part1), hat(phi)))),
                            as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")")))),
                 pch = c(16, 15), lty = c(1, NA), col = c(col_wv, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
        }
      }
    }
  }

  # Add WV
  lines(x$clusters, x$adev, type = "l", col = col_wv, pch = 16)

  if (is.null(point_pch)){
    point_pch = 16
  }

  if (is.null(point_cex)){
    point_cex = 1.25
  }
  lines(x$clusters, x$adev, type = "p", col = col_wv, pch = point_pch, cex = point_cex)
}

#' @title Computes the Allan Variance Linear Regression estimator
#'
#' @description
#' Compute the latent processes parameters estimator based on the Allan Variance
#' @param x     A \code{vec} of time series observations or an \code{avar} object.
#' @param qn    A \code{vec} specifying on which scales the parameters of a Quantization Noise (QN) should be computed.
#' @param wn    A \code{vec} specifying on which scales the parameters of a White Noise (WN) should be computed.
#' @param rw    A \code{vec} specifying on which scales the parameters of a Random Wakk (RW) should be computed.
#' @param dr    A \code{vec} specifying on which scales the parameters of a Drift (DR) should be computed.
#' @param type  A \code{string} containing either \code{"mo"} for Maximal Overlap or \code{"to"} for Tau Overlap
#' @return avlr   A \code{list} that contains:
#' \itemize{
#'  \item{"estimates"}
#'  \item{"implied_ad"}{The Allan Deviation implied by the estimated parameter.}
#'  \item{"implied_ad_decomp"}{The Allan Deviation implied by the estimated parameter for the sub-processes.}
#'  \item{"av"}{The \code{avar} object provided or corresponding to the data provided.}
#' }

#' @examples
#' # Load simts package
#' library(simts)
#'
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Simulate time series
#' N = 100000
#' x = gen_gts(N, WN(sigma2 = 1) + RW(gamma2 = 1e-7))
#'
#' # Maximal overlap
#' fit1 = avlr(x, wn = 1:12, rw = 12:15)
#'
avlr = function(x, qn = NULL, wn = NULL, rw = NULL, dr = NULL, type = NULL){

  if(is.null(x)){
    stop("Please provide a time series vector or a 'avar' objet")
  }else if(class(x)[1] != "avar"){
    if(is.null(type)){
      x = avar(x, type = "mo")
    }else{
      x = avar(x)
    }
  }

  if(sum(sapply(list(qn,wn,rw,dr), is.null)) == 4){
    stop("Please specify a least one process")
  }

  n_processes = sum(sapply(list(qn,wn,rw,dr), is.null))

  process = rep(NA,n_processes)
  param = rep(NA,n_processes)
  implied = matrix(NA,length(x$clusters),n_processes)

  counter = 1

    if(!is.null(wn)){
      if(length(wn) < 1 || !is.whole(wn) || min(wn) < 1 || max(wn) > length(x$allan)){
        stop("wn incorrectely formatted.")
      }
      process[counter] = "WN"
      param[counter] = exp(mean(log(x$adev[wn]) + log(x$clusters[wn])/2))
      implied[,counter] = param[counter]/sqrt(x$clusters)
      counter = counter + 1
    }

    if(!is.null(qn)){
      if(length(qn) < 1 || !is.whole(qn) || min(qn) < 1 || max(qn) > length(x$allan)){
        stop("qn incorrectely formatted.")
      }
      process[counter] = "QN"
      param[counter] = (1/sqrt(3))*exp(mean(log(x$adev[qn]) + log(x$clusters[qn])))
      implied[,counter] = sqrt(3)*param[counter]/(x$clusters)
      counter = counter + 1
    }

    if(!is.null(rw)){
      if(length(rw) < 1 || !is.whole(rw) || min(rw) < 1 || max(rw) > length(x$allan)){
        stop("rw incorrectely formatted.")
      }
      process[counter] = "RW"
      param[counter] = sqrt(3)*exp(mean(log(x$adev[rw]) - log(x$clusters[rw])/2))
      implied[,counter] = param[counter]*sqrt(x$clusters/3)
      counter = counter + 1
    }

    if(!is.null(dr)){
      if(length(dr) < 1 || !is.whole(dr) || min(dr) < 1 || max(dr) > length(x$allan)){
        stop("dr incorrectely formatted.")
      }
      process[counter] = "DR"
      param[counter] = sqrt(2)*exp(mean(log(x$adev[dr]) - log(x$clusters[dr])))
      implied[,counter] = param[counter]*x$clusters/2
      counter = counter + 1
    }

  implied_ad = apply(implied, 1, sum)

  estimates = t(t(param))
  rownames(estimates) = process
  colnames(estimates) = "Value"

  x = structure(list(estimates = estimates,
                     process_desc = process,
                     implied_ad = implied_ad,
                     implied_ad_decomp = implied,
                     av = x), class = "avlr")
  invisible(x)
}

#' Print gmwm object
#'
#' Displays information about GMWM object
#' @method print avlr
#' @export
#' @keywords internal
#' @param x   A \code{avlr} object
#' @param ... Other arguments passed to specific methods
#' @return Text output via print
#' @examples
#' \dontrun{
#' # Load simts package
#' library(simts)
#'
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Simulate time series
#' N = 100000
#' x = gen_gts(N, WN(sigma2 = 1) + RW(gamma2 = 1e-7))
#'
#' # Maximal overlap
#' fit1 = avlr(x, wn = 1:12, rw = 12:15)
#' print(fit1)
#' }
print.avlr = function(x, ...) {
  cat("\n Estimates: \n")
  print(x$estimates)
}

#' @title Plot Allan Variance Linear Regression Fit
#'
#' @description
#' Displays a plot of allan deviation accounting for CI values with the AD implied by the estimated parameters
#' @method plot avar
#' @param x                A \code{avlr} object.
#' @param decomp           A \code{boolean} that determines whether the latent proceses individual contributions are plotted.
#' @param units            A \code{string} that specifies the units of time plotted on the x axis.
#' @param xlab             A \code{string} that gives a title for the x axis.
#' @param ylab             A \code{string} that gives a title for the y axis.
#' @param main             A \code{string} that gives an overall title for the plot.
#' @param col_wv           A \code{string} that specifies the color of the wavelet variance line.
#' @param col_ci           A \code{string} that specifies the color of the confidence interval polygon.
#' @param ci_wv            A \code{boolean} that determines whether a confidence interval polygon will be drawn.
#' @param nb_ticks_x       An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y       An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param legend_position  A \code{string} that specifies the position of the legend (use \code{legend_position = NA} to remove legend).
#' @param point_pch        A \code{double} that specifies the symbol type to be plotted.
#' @param point_cex        A \code{double} that specifies the size of each symbol to be plotted.
#' @param ...              Additional arguments affecting the plot.
#' @return Plot of allan deviation and confidence interval for each scale.
#' @author Stephane Guerrier, Nathanael Claussen, and Justin Lee
#' @export
#' @examples
#' # Load simts package
#' library(simts)
#'
#' # Set seed for reproducibility
#' set.seed(999)
#'
#' # Simulate time series
#' N = 100000
#' ts = gen_gts(N, WN(sigma2 = 1) + RW(gamma2 = 1e-7))
#'
#' x = avlr(ts, wn = 1:12, rw = 12:15)
#'
#' # Plot example
#' plot.avlr(x)
#' plot.avlr(x, decomp = TRUE, main = "Simulated white noise", xlab = "Scales")
#' plot.avlr(x, units = "sec", legend_position = "topright")
#' plot.avlr(x, col_wv = "darkred", col_ci = "pink")
plot.avlr = function(x, decomp = FALSE,
                      units = NULL, xlab = NULL, ylab = NULL, main = NULL,
                      col_wv = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                      legend_position = NULL, ci_wv = NULL, point_cex = NULL,
                      point_pch = NULL, ...){


  # Labels
  if (is.null(xlab)){
    if (is.null(units)){
      xlab = expression(paste("Scale ", tau, sep =""))
    }else{
      xlab = bquote(paste("Clustering time ", tau, " [", .(units), "]", sep = " "))
    }
  }

  if (is.null(ylab)){
    ylab = expression(paste("Allan Deviation ", phi, sep = ""))
  }else{
    ylab = ylab
  }

  # Main Title
  if (is.null(main)){
    main = "Allan Variance Representation"
  }

  # Line and CI colors
  if (is.null(col_wv)){
    col_wv = "darkblue"
  }

  if (is.null(col_ci)){
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }

  # Range
  x_range = range(x$av$clusters)
  x_low = floor(log2(x_range[1]))
  x_high = ceiling(log2(x_range[2]))

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
  x_labels = sapply(x_ticks, function(i) as.expression(bquote(2^ .(i))))

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
  abline(v = 2^x_ticks, lty = 1, col = "grey95")
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
  axis(1, at = 2^x_ticks, labels = x_labels, padj = 0.3)
  axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2)

  # CI for WV
  if (ci_wv == TRUE || is.null(ci_wv)){
    polygon(c(x$av$cluster, rev(x$av$cluster)), c(x$av$adev - x$av$errors*x$av$adev, rev(x$av$adev + x$av$errors*x$av$adev)),
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
      lines(x$av$clusters, x$implied_ad_decomp[,i], col = col_decomp[i])
    }
  }
  # Plot implied AD
  lines(t(x$av$clusters),x$implied_ad, type = "l", lwd = 3, col = "#F47F24", pch = 1, cex = 1.5)
  lines(t(x$av$clusters),x$implied_ad, type = "p", lwd = 2, col = "#F47F24", pch = 1, cex = 1.5)

  # Add WV
  lines(x$av$clusters, x$av$adev, type = "l", col = col_wv, pch = 16)

  if (is.null(point_pch)){
    point_pch = 16
  }

  if (is.null(point_cex)){
    point_cex = 1.25
  }
  lines(x$av$clusters, x$av$adev, type = "p", col = col_wv, pch = point_pch, cex = point_cex)

  # Add legend
  CI_conf = .95
  wv_title_part1 = "Empirical AV "

  if(decomp == TRUE){
    legend_names = c(as.expression(bquote(paste(.(wv_title_part1), hat(phi)))),
                     as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")"))),"Implied AV",
                     x$process_desc)
    col_legend = c(col_wv, col_ci,"#F47F24",col_decomp)
    p_cex_legend = c(1.25, 3, 1.5,rep(NA,U))
    lty_legend = c(1, NA, rep(1,U))
    pch_legend = c(16,15,1,rep(NA,U))
  }else{
    legend_names = c(as.expression(bquote(paste(.(wv_title_part1), hat(phi)))),
                     as.expression(bquote(paste("CI(",hat(phi),", ",.(CI_conf),")"))),"Implied AV")
    col_legend = c(col_wv, col_ci,"#F47F24")
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

