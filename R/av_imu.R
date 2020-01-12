avar_imu = function(x, type = "mo"){
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
    av[[i]] = avar(x[,i], type = type, freq = freq)
  }
  names(av) = col_names
  out = list(sensor = sensor_name, freq = freq, n = n, type = sensor, axis = ax, avar = av)
  class(out) = "imu_avar"
  invisible(out)
}
