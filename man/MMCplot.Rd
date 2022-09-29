\name{MMCplot}
\alias{MMCplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting Multivariate Markov Chains' Network
}
\description{
Plots the Multivariate Markov Chains obtained through the estimation of the Generalized Multivariate Markov Chain Model, through \code{mmcx()}, for a specific value of the exogeneous variable considered.
}
\usage{
MMCplot(s, x, value = max(x), eq = 1, result)
}
\arguments{
  \item{s}{
Matrix of categorical data sequences.
}
  \item{x}{
Matrix of covariates (exogeneous variables).
}
  \item{value}{
Value of \code{x} for which the network is plotted.
}
  \item{eq}{
Number of the equation to be considered in the network.
}
  \item{result}{
Output object obtained from \code{mmcx()}.
}
}
\value{
The result of the function is a plot of the igraph object produced from the Markov chain, considering the user-defined value of the exogeneous variable.
}
\references{
Ching, W. K., E. S. Fung, and M. K. Ng (2002). A multivariate Markov chain modelfor categorical data sequences and its applications in demand predictions. IMA Journal of Management Mathematics, 13(3), 187-199. \doi{10.1093/imaman/13.3.187}

Raftery, A. E. (1985). A Model for High-Order Markov Chains. Journal of the Royal Statistical Society. Series B (Methodological), 47(3), 528-539. \url{http://www.jstor.org/stable/2345788}

Vasconcelos, C. and Damasio, B. (2022). GenMarkov: Modeling Generalized Multivariate Markov Chains in R. \url{https://arxiv.org/abs/2202.00333v1}
}
\author{
Carolina Vasconcelos and Bruno Damasio
}
\examples{
set.seed(1234)
s1 <- sample(c(1,2), 500, replace=TRUE)
s2 <- sample(c(1,2), 500, replace=TRUE)
x <- rnorm(500)
res <- mmcx(y = cbind(s1,s2), x = cbind(x), initial=c(1,1))

MMCplot(s=cbind(s1,s2), x=cbind(x), value=2, eq=1, result = res)
}
