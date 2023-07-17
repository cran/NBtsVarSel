\name{variable_selection}
\alias{variable_selection}
\title{
Variable selection
}
\description{
This function performs variable selection, estimates new vectors of beta and gamma and a new alpha 
}
\usage{
variable_selection(Y, X, gamma.init, alpha.init = NULL, k.max = 1, method = "cv", 
tr = 0.3, n.iter = 100, n.rep = 1000)
}
\arguments{
  \item{Y}{
Observation matrix
}
  \item{X}{
Design matrix
}
  \item{gamma.init}{
Initial gamma vector
}
  \item{alpha.init}{
Optional initial alpha value. The default is NULL
}
  \item{k.max}{
Number of iteration to repeat the whole algorithm
}
  \item{method}{
Stability selection method: "min" or "cv". In "min" the smallest lambda is chosen, in "cv" cross-validation lambda is chosen for stability selection. The default is "cv"
}
  \item{tr}{
Threshold for stability selection. The default is 0.3
}
  \item{n.iter}{
Number of iteration for Newton-Raphson algorithm. The default is 100
}
  \item{n.rep}{
Number of replications in stability selection step. The default is 1000
}
}
\value{
\item{estim_active }{Estimated active coefficients}
\item{beta_est}{Vector of estimated beta values}
\item{gamma_est}{Vector of estimated gamma values}
\item{alpha_est}{Estimation of alpha}
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
result = variable_selection(Y, X, gamma.init=gamma0, alpha.init=NULL, k.max=1, method="cv", 
tr=0.3, n.iter=100, n.rep=1000)
beta_est = result$beta_est
Estim_active = result$estim_active
gamma_est = result$gamma_est
alpha_est = result$alpha_est
}
