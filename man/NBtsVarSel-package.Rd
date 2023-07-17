\name{NBtsVarSel-package}
\alias{NBtsVarSel-package}
\alias{NBtsVarSel}
\docType{package}
\title{
\packageTitle{NBtsVarSel}
}
\description{
NBtsVarSel consists of four functions: "variable_selection.R", "grad_hess_beta.R", "grad_hess_gamma.R" and "NR_gamma.R"
For further information on how to use these functions, 
we refer the reader to the vignette of the package.
}
\details{
This package consists of four functions: "variable_selection.R", "grad_hess_beta.R", "grad_hess_gamma.R" and "NR_gamma.R"
For further information on how to use these functions, 
we refer the reader to the vignette of the package.
}
\author{
Marina Gomtsyan

Maintainer: Marina Gomtsyan <mgomtsian@gmail.com>
}
\references{
M. Gomtsyan "Variable selection in a specific regression time series of counts.", arXiv:2307.00929
}
\keyword{ package }
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
