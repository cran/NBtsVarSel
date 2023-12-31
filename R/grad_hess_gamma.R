grad_hess_gamma <-
function(Y, X, beta, gamma, alpha)
{
  q = length(gamma)
  p = length(X[,1])-1
  n = length(Y)
  mu = numeric(n)
  E = numeric(n) 
  Z = numeric(n) 
  W = numeric(n) 
  
  grad_W = matrix(0,nrow=q,ncol=n)
  hess_W_liste = list()
  hess_L_gamma = matrix(0,ncol=q,nrow=q)
  
  Z[1] = 0
  W[1] = beta %*% X[,1]
  mu[1] = exp(W[1])
  E[1] = (Y[1] - mu[1])/(mu[1] + mu[1]^2 / alpha)
  hess_W_liste[[1]] = matrix(0,ncol=q,nrow=q)
  
  for (t in 2:n)
  {
    hess_W_liste[[t]] = matrix(NA,ncol=q,nrow=q)
    jsup = min(q,(t-1))
    W[t] = beta %*% X[,t] + sum(gamma[1:jsup]*E[(t-1):(t-jsup)])
    
    for (l in 1:q)
    {
      if (l<=jsup) 
      {
        grad_W_1 = E[(t-1):(t-jsup)] + (1 + E[(t-1):(t-jsup)] * mu[(t-1):(t-jsup)] /alpha) / (1 + mu[(t-1):(t-jsup)] /alpha) 
        grad_W[l,t] = E[(t-l)] - sum(gamma[1:jsup]*grad_W_1*grad_W[l,(t-1):(t-jsup)])
      }
      for (m in l:q)
      {
        if (m+l<t) 
        {
          A_1 = - (E[(t-l)] + (1 + E[(t-l)] * mu[(t-l)] /alpha) / (1 + mu[(t-l)] /alpha)) * grad_W[m,(t-l)]
          A_2 = - (E[(t-m)] + (1 + E[(t-m)] * mu[(t-m)] /alpha) / (1 + mu[(t-m)] /alpha)) * grad_W[l,(t-m)]
        }
        else 
        {
          A_1=0
          A_2=0
        }
        A = A_1+A_2
        
        C_1 = E[(t-1):(t-jsup)] + (1 + E[(t-1):(t-jsup)] * mu[(t-1):(t-jsup)] /alpha) / (1 + mu[(t-1):(t-jsup)] /alpha)
        B_2_1_den = alpha * (1 + mu[(t-1):(t-jsup)] /alpha)^2
        B_2_1 = 2 * (E[(t-1):(t-jsup)] * exp(2*W[(t-1):(t-jsup)])/alpha + Y[(t-1):(t-jsup)])
        B_2_2 = (1 - E[(t-1):(t-jsup)] * mu[(t-1):(t-jsup)] /alpha) / (1 + mu[(t-1):(t-jsup)] /alpha)
        B_2 = E[(t-1):(t-jsup)] + B_2_1 / B_2_1_den + B_2_2
        
        B = B_2*grad_W[l,(t-1):(t-jsup)]*grad_W[m,(t-1):(t-jsup)]
        C = -C_1*as.numeric((lapply(hess_W_liste[((t-1):(t-jsup))],'[',l,m)))
        hess_W_liste[[t]][l,m] = A + sum(gamma[1:jsup]*(B+C))
      }
    }
    mu[t] = exp(W[t])
    E[t] = (Y[t]-mu[t])/(mu[t] + mu[t]^2 / alpha)
    
    term1_1_hess = (alpha + Y[t]) * exp(W[t]) / (alpha + exp(W[t]))
    terme1_hess = hess_W_liste[[t]] * (Y[t] - term1_1_hess)
    term2_1_hess = term1_1_hess * (1 - exp(W[t]) / (alpha + exp(W[t])))
    terme2_hess = term2_1_hess * (grad_W[,t] %*%t (grad_W[,t]))
    hess_L_gamma = hess_L_gamma + terme1_hess - terme2_hess
  }
  hess_L_gamma[lower.tri(hess_L_gamma)] = hess_L_gamma[upper.tri(hess_L_gamma)]
  term1_grad = (alpha + Y) * exp(W) / (alpha + exp(W))
  grad_L_gamma = grad_W %*% (Y- term1_grad)
  return(list(grad_L_gamma=grad_L_gamma, hess_L_gamma=hess_L_gamma))
}
