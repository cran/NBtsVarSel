variable_selection <-
function(Y, X, gamma.init, alpha.init=NULL, k.max=1, method="cv", tr=0.3, n.iter=100, n.rep=1000)
{
  n = length(Y)
  p = dim(X)[1]-1
  q = length(gamma.init)
  
  if((p+1)<n){
    glm_nb = glm.nb(Y~t(X)[,2:(p+1)])
    beta.init = as.numeric(glm_nb$coefficients)
    theta = glm_nb$theta
  }else{
    new_data = data.frame(cbind(Y, t(X)[,2:(p+1)]))
    lambda_cv = cv.glmregNB(Y~., new_data, alpha=1)$lambda.optim
    glmnet_nb = glmregNB(Y~t(X), alpha = 1, lambda = lambda_cv)
    theta = glmnet_nb$theta
    beta.init = as.numeric(glmnet_nb$beta)
    beta.init[1] = as.numeric(glmnet_nb$b0)
  }
  
  if(is.null(alpha.init)){
    alpha.init = theta
  }
  
  for(k in 1:k.max){
    gamma_est = NR_gamma(Y, X, beta.init, gamma.init, alpha.init, n.iter)
    
    grad_hess_res_est = grad_hess_beta(Y, X, beta.init, gamma.init, alpha.init)
    Gradient = grad_hess_res_est[[1]][1:(p+1)]
    Hessienne = grad_hess_res_est[[2]]
    
    res_svd = svd(-Hessienne) 
    U = res_svd$u
    ind_vp_not_null = which(round(res_svd$d,digits=6)!=0)
    Lambda_rac = diag(sqrt(res_svd$d[ind_vp_not_null])) 
    Lambda_rac_inv = diag(1/sqrt(res_svd$d[ind_vp_not_null]))
    Y_eta = Lambda_rac_inv%*%t(U[,ind_vp_not_null])%*%Gradient+Lambda_rac%*%t(U[,ind_vp_not_null])%*%beta.init
    X_eta = Lambda_rac%*%t(U[,ind_vp_not_null])
    
    
    if(method=='min'){
      all_lambda=glmnet(X_eta,Y_eta,family="gaussian",alpha=1,parallel=TRUE)$lambda
      lambda_cv = c(min(all_lambda))
    }else if(method=='cv'){
      lambda_cv = cv.glmnet(X_eta,Y_eta,family="gaussian",alpha=1,parallel=TRUE)$lambda.min
    }
    
    res.cum = rep(0,(p+1))
    b_sort_matrix = matrix(NA, nrow = n.rep, ncol = floor((p+1)/2))
    for(j in 1:n.rep){
      b_sort_matrix[j,] = sort(sample(1:(p+1),floor((p+1)/2)))
    }
    for(j in 1:n.rep){
      b_sort = b_sort_matrix[j,]
      resultat_glmnet = glmnet(X_eta[b_sort,],Y_eta[b_sort],family="gaussian",alpha=1,lambda=lambda_cv)
      ind_glmnet = which(as.numeric(resultat_glmnet$beta)!=0)
      res.cum = res.cum + tabulate(ind_glmnet,(p+1))
    }
    freq = res.cum/n.rep
    Estim_active = which(freq>=tr)
    
    beta_est=rep(0,(p+1))
    if(length(Estim_active) == 1){
      glm_nb_1 = glm.nb(Y~t(X)[,Estim_active]-1)
      alpha_est = glm_nb_1$theta
      beta_est[Estim_active]<-as.numeric(glm_nb_1$coefficients)
    }else if((length(Estim_active)<n) & (rankMatrix(t(X)[,Estim_active])[1]==length(Estim_active))){
      glm_nb_1 = glm.nb(Y~t(X)[,Estim_active]-1)
      alpha_est = glm_nb_1$theta
      beta_est[Estim_active] = as.numeric(glm_nb_1$coefficients)
    }else{
      new_data_1 = data.frame(cbind(Y, t(X)[,Estim_active]))
      lambda_cv_1 = cv.glmregNB(Y~., alpha=0.01, new_data_1, plot.it=FALSE)$lambda.optim
      glmnet_nb_1 = glmregNB(Y~t(X)[,Estim_active], alpha = 0, lambda = lambda_cv_1)
      alpha_est = glmnet_nb_1$theta
      beta_est[Estim_active]<-as.numeric(glmnet_nb_1$beta)
      if((1%in% Estim_active)==TRUE){
        beta_est[1] = as.numeric(glmnet_nb_1$b0)
      }
    }
    beta.init = beta_est
    gamma.init = gamma_est
    alpha.init = alpha_est
  }
  
  return(list(estim_active=Estim_active,beta_est=beta_est,gamma_est=gamma_est, alpha_est=alpha_est))
}
