#' @title Allan variance (maximum-overlap estimator) of 4 Hours from a navchip sensor
#' @description This data set contains Allan variance of gyroscope and accelerometer data from a navchip-based sensors.
#' @format An imu with 3,105,250 observations and 6 sensors. The sensors are defined as follows:
#' \describe{
#'  \item{\code{sensor}}{Name of the sensor.}
#'  \item{\code{freq}}{The frequency at which the error signal is measured.}
#'  \item{\code{n}}{Sample size of the data.}
#'  \item{\code{type}}{The types of sensors considered in the data.}
#'  \item{\code{axis}}{The axes of sensors considered in the data.}
#'  \item{\code{avar}}{A list containing the computed Allan variance based on the data.}
#' }
#' @export
#' @author James Balamuta, Stephane Guerrier and Yuming Zhang
#' @source The IMU data of the navchip sensor comes from Geodetic Engineering Laboratory (TOPO) and Swiss Federal Institute of Technology Lausanne (EPFL).
"navchip_avar"
