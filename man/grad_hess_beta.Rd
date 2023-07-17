\name{grad_hess_beta}
\alias{grad_hess_beta}
\title{
Gradient and Hessian of the log-likelihood with respect to beta
}
\description{
This function calculates the gradient and Hessian of the log-likelihood with
respect to beta.
}
\usage{
grad_hess_beta(Y, X, beta, gamma, alpha)
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
}
\value{
\item{grad_L_beta}{
Vector of the gradient of L with respect to beta
}
\item{hess_L_beta}{
Matrix of the Hessian of L with respect to beta
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
result = grad_hess_beta(Y, X, beta0, gamma0, alpha0)
grad = result$grad_L_beta
Hessian = result$hess_L_beta
}