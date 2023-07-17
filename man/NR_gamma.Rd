\name{NR_gamma}
\alias{NR_gamma}
\title{
Newton-Raphson method for estimation of gamma
}
\description{
This function estimates gamma with Newton-Raphson method
}
\usage{
NR_gamma(Y, X, beta, gamma, alpha, n.iter)
}
\arguments{
  \item{Y}{
Observation matrix
}
  \item{X}{
Design matrix
}
  \item{beta}{
Initial beta vector
}
  \item{gamma}{
Initial gamma vector
}
  \item{alpha}{
Initial overdispertion parameter
}
  \item{n.iter}{
Number of iterations of the algorithm. Default=100
}
}
\value{
\item{gamma}{
Estimated gamma vector
}
}
\references{
M. Gomtsyan "Variable selection in a specific regression time series of counts.", arXiv:2307.00929
}
\author{
Marina Gomtsyan

Maintainer: Marina Gomtsyan <mgomtsian@gmail.com>
}
\examples{
n = 50
p = 30
X = matrix(NA,(p+1),n)
f = 1/0.7
for(t in 1:n){X[,t] = c(1,cos(2*pi*(1:(p/2))*t*f/n),sin(2*pi*(1:(p/2))*t*f/n))}
gamma0 = c(0)
data(Y)
glm_nb = glm.nb(Y~t(X)[,2:(p+1)])
beta0 = as.numeric(glm_nb$coefficients)
alpha0 = glm_nb$theta
gamma_est = NR_gamma(Y, X, beta0, gamma0, alpha0, n.iter=100)
}
