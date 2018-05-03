#' @export
avar = function(x, type = "mo") {
  x = as.vector(x)

  if(type == "mo"){
    av = avar_mo_cpp(x)
  }else{
    av = avar_to_cpp(x)
  }

  av = list(clusters = av[,1], allan=av[,2], errors=av[,3])
  av$adev = sqrt(av$allan)
  av$lci = av$adev - av$errors*av$adev
  av$uci = av$adev + av$errors*av$adev
  av$type = type
  class(av) = c("avar", "list")
  av
}


#' @export
print.avar = function(x, ...) {
  cat("\n Clusters: \n")
  print(x$clusters, digits=5)
  cat("\n Allan Variances: \n")
  print(x$allan, digits=5)
  cat("\n Errors: \n")
  print(x$errors, digits=5)
}

#' @export
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
