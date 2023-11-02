#' Title
#'
#' @param hess numerical matrix
#' @param lambda numerical vector
#'
#' @return list with standard error, z-statistics and pvalues associated with each estimate
#' @export
#'
#' @examples
#' ######
Inference_x <- function(hess, lambda) {

  if(matrixcalc::is.singular.matrix(hess, tol = 1e-05) == FALSE){
    hessinv <- solve(-hess)

    var <- diag(hessinv)

    ifelse(any(var<0),
      return(l = list(
        warning = 1,
        se = '.',
        zstat = '.',
        pvalue = '.'
      )),
      return(l = list(
        warning = 0,
        se = sqrt(var),
        zstat = lambda / sqrt(var),
        pvalue =  2 * (1 - stats::pnorm(abs(lambda / sqrt(var))))
      ))
    )

  }else{
    return(l = list(
      warning = 1,
      se = '.',
      zstat = '.',
      pvalue = '.'
    ))
  }

}
