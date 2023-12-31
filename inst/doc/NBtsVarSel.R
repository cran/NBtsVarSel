## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE)
library(NBtsVarSel)
library(ggplot2)
library(formatR)
set.seed(555)

## ----Y------------------------------------------------------------------------
data(Y)

## ----dimensions, echo=FALSE, eval=TRUE----------------------------------------
n = length(Y)
p = 30

## ----gamma_alpha, echo=FALSE, eval=TRUE---------------------------------------
gamma = c(0.5)
alpha = 2

## ----beta, echo=FALSE, eval=TRUE----------------------------------------------
active=c(1,3,17,21,23)
beta_t_pos=c(1.73,0.38,0.29,-0.64,-0.13)
beta = rep(0,(p+1))
beta[active] = beta_t_pos

## ----X, echo=FALSE, eval=TRUE-------------------------------------------------
X = matrix(NA,(p+1),n)
f = 1/0.7
for(t in 1:n){X[,t]<-c(1,cos(2*pi*(1:(p/2))*t*f/n),sin(2*pi*(1:(p/2))*t*f/n))}

## ----initialization-----------------------------------------------------------
gamma0 = c(0)

glm_nb = glm.nb(Y~t(X)[,2:(p+1)])
beta0<-as.numeric(glm_nb$coefficients)

alpha0 = glm_nb$theta

## ----gammaEst-----------------------------------------------------------------
gamma_est_nr = NR_gamma(Y, X, beta0, gamma0, alpha0, n.iter=100)
gamma_est_nr

## ----variableSelection, tidy=TRUE, tidy.opts=list(width.cutoff=60)------------
result = variable_selectionresult = variable_selection(Y, X, gamma.init=gamma0, alpha.init=NULL, k.max=1, method="cv", t=0.3, n.iter=100, n.rep=1000)
beta_est = result$beta_est
Estim_active = result$estim_active
gamma_est = result$gamma_est
alpha_est = result$alpha_est

## ----print, echo=FALSE, eval=TRUE---------------------------------------------
cat("Estimated active coefficients: ", Estim_active, "\n")
cat("Estimated gamma: ", gamma_est, "\n")
cat("Estimated alphaa: ", alpha_est, "\n")

## ----plot, fig.width=6, fig.height=4, tidy=TRUE, tidy.opts=list(width.cutoff=54)----
#First, we make a dataset of estimated betas
beta_data = data.frame(beta_est)
colnames(beta_data)[1] <- "beta"
beta_data$Variable = seq(1, (p+1), 1)
beta_data$y = 0
beta_data = beta_data[beta_data$beta!=0,]
#Next, we make a dataset of true betas
beta_t_data = data.frame(beta)
colnames(beta_t_data)[1] <- "beta"
beta_t_data$Variable = seq(1, (p+1), 1)
beta_t_data$y = 0
beta_t_data = beta_t_data[beta_t_data$beta!=0,]
#Finally, we plot the result
plot = ggplot()+
  geom_point(data = beta_data, aes(x=Variable, y=y, color=beta), pch=16, size=5, stroke = 2)+
  geom_point(data= beta_t_data, aes(x=Variable, y=y, color=beta), pch=4, size=6, stroke = 2)+
  scale_color_gradient2(name=expression(hat(beta)), midpoint=0, low="steelblue", mid = "white", high ="red")+
  scale_x_continuous(breaks=c(1, seq(10, (p+1), 10)), limits=c(0, (p+1)))+
  scale_y_continuous(breaks=c(), limits=c(-1, 1))+
  theme(legend.title = element_text(color = "black", size = 12, face="bold"), legend.text = element_text(color = "black", size = 10))+
  theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold"), axis.title.y=element_blank())
plot

