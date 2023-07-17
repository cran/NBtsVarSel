grad_hess_beta <-
function(Y, X, beta, gamma, alpha)
{
  q = length(gamma)
  p = length(X[,1])-1
  n = length(Y)
  mu = numeric(n) 
  E = numeric(n) 
  Z = numeric(n) 
  W = numeric(n) 
  
  grad_W_beta = matrix(NA,nrow=(p+1),ncol=n)
  hess_W_beta_liste = list()
  hess_L_beta = matrix(0,ncol=(p+1),nrow=(p+1))
  
  Z[1] = 0
  W[1] = beta %*% X[,1]
  mu[1] = exp(W[1])
  E[1] = (Y[1] - mu[1])/(mu[1] + mu[1]^2 / alpha)
  hess_W_beta_liste[[1]] = matrix(0,ncol=(p+1),nrow=(p+1))
  grad_W_beta[,1] = X[,1]
  
  for (t in 2:n)
  {
    hess_W_beta_liste[[t]] = matrix(NA,ncol=(p+1),nrow=(p+1))
    jsup = min(q,(t-1))
    W[t] = beta%*%X[,t] + sum(gamma[1:jsup]*E[(t-1):(t-jsup)])
    for (k in 1:(p+1))
    {
      grad_W_beta_1 = E[(t-1):(t-jsup)] + (1 + E[(t-1):(t-jsup)] * mu[(t-1):(t-jsup)] /alpha) / (1 + mu[(t-1):(t-jsup)] /alpha)
      grad_W_beta[k,t] = X[k,t]-sum(gamma[1:jsup] * grad_W_beta_1 * grad_W_beta[k,(t-1):(t-jsup)])
      
      for (j in k:(p+1))
      {
        A_1 = E[(t-1):(t-jsup)] + (1 + E[(t-1):(t-jsup)] * mu[(t-1):(t-jsup)] /alpha) / (1 + mu[(t-1):(t-jsup)] /alpha)
        A_2_1_den = alpha * (1 + mu[(t-1):(t-jsup)]/alpha)^2
        A_2_1 = 2 * (E[(t-1):(t-jsup)] * exp(2*W[(t-1):(t-jsup)])/alpha + Y[(t-1):(t-jsup)])
        A_2_2 = (1 - E[(t-1):(t-jsup)] * mu[(t-1):(t-jsup)] /alpha) / (1 + mu[(t-1):(t-jsup)] /alpha)
        A_2 = E[(t-1):(t-jsup)] + A_2_1 / A_2_1_den + A_2_2
        
        B = A_2*grad_W_beta[k,(t-1):(t-jsup)]*grad_W_beta[j,(t-1):(t-jsup)]
        C = -A_1*as.numeric((lapply(hess_W_beta_liste[((t-1):(t-jsup))],'[',k,j)))
        hess_W_beta_liste[[t]][k,j] = sum(gamma[1:jsup]*(B+C))
      }
    }
    mu[t] = exp(W[t])
    E[t] = (Y[t]-mu[t])/(mu[t] + mu[t]^2 / alpha)
    
    term1_1_hess = (alpha + Y[t]) * exp(W[t]) / (alpha + exp(W[t]))
    terme1_hess = hess_W_beta_liste[[t]] * (Y[t] - term1_1_hess)
    term2_1_hess = term1_1_hess * (1 - exp(W[t]) / (alpha + exp(W[t])))
    terme2_hess = term2_1_hess * (grad_W_beta[,t] %*%t (grad_W_beta[,t]))
    hess_L_beta = hess_L_beta + terme1_hess - terme2_hess
  }
  hess_L_beta[lower.tri(hess_L_beta)] = hess_L_beta[upper.tri(hess_L_beta)]
  term1_grad = (alpha + Y) / (alpha + exp(W))
  grad_L_beta = grad_W_beta %*% (Y-exp(W) * term1_grad)
  return(list(grad_L_beta=grad_L_beta, hess_L_beta=hess_L_beta))
}
