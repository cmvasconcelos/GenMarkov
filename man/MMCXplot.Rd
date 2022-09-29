\name{MMCXplot}
\alias{MMCXplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting Non-homogeneous Multivariate Markov Chains
}
\description{
Plots the Non-homogeneous Multivariate Markov Chains considered in the estimation of the Generalized Multivariate Markov Chain Model, through \code{mmcx()}.
}
\usage{
MMCXplot(s, x, mov = TRUE, order = 5, wdth = 7, hght = 7, unts = "in")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
Matrix of categorical data sequences.
}
  \item{x}{
Matrix of covariates (exogeneous variables).
}
  \item{mov}{
Logical value. If \emph{TRUE}, then it will transform the non-homogeneous Markov chains through a moving average.
}
  \item{order}{
Order of the moving average transformation.
}
  \item{wdth}{
Width of the output plots.
}
  \item{hght}{
Height of the output plots.
}
  \item{unts}{
Units considered for the height and width of the plots. For more information see \code{\link{gridExtra}}.
}
}
\value{
The function returns the plots (with png extension) at the working directory defined by the user.
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
#set.seed(1234)
#s1 <- sample(c(1,2), 500, replace=TRUE)
#s2 <- sample(c(1,2), 500, replace=TRUE)
#x <- rnorm(500)
#mmcx(y = cbind(s1,s2), x = cbind(x), initial=c(1,1))

#MMCXplot(s=cbind(s1,s2), x=cbind(x), mov=FALSE)

}

