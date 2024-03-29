#Initial Values p. 387 de Berchtold (2001) - (This function is adapted from march::march.mtd.construct())
#' Title
#'
#' @param f a numerical array
#' @param split_result a numerical vector, to split indexes of array
#'
#' @return a numerical vector
#' @export
#'
#' @examples
#' set.seed(1234)
#' s1 <- sample(c(1,2), size=500, replace=TRUE)
#' s2 <- sample(c(1,2), size=500, replace=TRUE)
#' f <- ProbMatrixF(cbind(s1,s2))
#' CalculateInitialValues(f, c(1,2)) ##for first equation
#' CalculateInitialValues(f, c(3,4)) ##for second equation
CalculateInitialValues <- function(f, split_result = split_result) {

  u <- apply(f, 3, FUN = function(cg){
    tc <- sum(cg)
    sr <- rowSums(cg)
    sc <- colSums(cg)

    num <- sum(cg*log2(sc) + t(t(cg)*log2(sr)) - cg*log2(cg) - cg*log2(tc))

    den <- sum(sc*log2(sc) - sc*log2(tc))

    u <- num/den
    u[is.infinite(u)] <- 0

    return(u)
  })

  lambda <- lapply(split_result, FUN = function(x) u[x]/sum(u[x]))

  return(lambda = unlist(lambda))
}
