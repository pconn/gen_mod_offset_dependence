
library(mvtnorm)
library(lme4)
set.seed(12345)

mu_lambda_log = log(10)
mu_exposure_log = log(1)
sigma_lambda_log = 0.5
sigma_exposure_log = 0.5
cor_link = .5
cov_link = sigma_lambda_log*sigma_exposure_log*cor_link
Sigma = matrix(c(sigma_lambda_log^2,cov_link,cov_link,sigma_exposure_log^2),2,2)

n_sim = 10000
n_samp = 100
Bias_reg = Bias_ht = Bias_mle = Bias_wtd = Bias_glmm = Bias_glmm_wtd =  Bias_glmm_nobiascorr = Bias_glmm_ht = Cor = rep(0,n_sim)

pois_negloglik <- function(par,offset,data){
  lambda=exp(par)
  -sum(data*(par+log(offset))-lambda*offset)
}


set.seed(12345)
for(isim in 1:n_sim){
  Pars = rmvnorm(n_samp,c(mu_lambda_log,mu_exposure_log),Sigma)
  Lam_i = exp(Pars[,1])
  exposure_i = exp(Pars[,2])
  y_i = rpois(n_samp,Lam_i*exposure_i)
  ht_i = y_i/exposure_i
  Data = data.frame("y_i"=y_i,"exposure_i"=exposure_i,"ht_i"=ht_i)
  
  pois_glm <- glm(y_i ~ offset(log(exposure_i)),data=Data,family="poisson")
  New_Data=Data
  New_Data$exposure_i = 1
  pred_glm = exp(pois_glm$coefficients)
  
  #pois_glm2 <- glm(y_i ~ offset(log(exposure_i)),data=Data,family="poisson",method=my.glm.fit)
  #pred_glm2 = exp(pois_glm2$coefficients)
  
  pois_est_mle <- nlminb(0,pois_negloglik,offset=Data$exposure_i,data=Data$y_i)
  pred_mle = exp(pois_est_mle$par)
  
  True_Lam = mean(Lam_i)
  ht_glm <- glm(ht_i ~ 1, family="poisson")
  pred_ht = exp(ht_glm$coefficients)
  
  Data$obs = c(1:n_samp)
  pois_glmm <- glmer(y_i~offset(log(exposure_i))+(1 | obs),family="poisson",data=Data)
  pred_glmm = exp(fixef(pois_glmm)+VarCorr(pois_glmm)$obs[1]/2)
  
  Cor[isim]=cor(Lam_i,exposure_i)
  Bias_reg[isim] = (pred_glm-True_Lam)/True_Lam
  #Bias_reg2[isim] = (pred_glm2-True_Lam)/True_Lam
  Bias_ht[isim] = (pred_ht - True_Lam)/True_Lam
  Bias_mle[isim] = (pred_mle - True_Lam)/True_Lam
  Bias_glmm[isim]=(pred_glmm - True_Lam)/True_Lam
  
  pois_reg_wtd <- glm(y_i ~ offset(log(exposure_i)),weights=1/exposure_i,data=Data,family="poisson")
  pred_glm_wtd = exp(pois_reg_wtd$coefficients)
  Bias_wtd[isim]= (pred_glm_wtd - True_Lam)/True_Lam

  pois_glmm_wtd <- glmer(y_i~offset(log(exposure_i))+(1 | obs),family="poisson",data=Data,weights=1/exposure_i)
  pred_glmm_wtd = exp(fixef(pois_glmm_wtd)+VarCorr(pois_glmm_wtd)$obs[1]/2)
  Bias_glmm_wtd[isim]=(pred_glmm_wtd - True_Lam)/True_Lam
  
  pred_glmm_nobiascorr = exp(fixef(pois_glmm_wtd))
  Bias_glmm_nobiascorr[isim]=(pred_glmm_nobiascorr - True_Lam)/True_Lam
  
  #skip since glmer seems to have strange behavior when poisson responses are non-integer
  #pois_glmm_ht <- glmer(ht_i~(1 | obs),family="poisson",data=Data)
  #pred_glmm_ht = exp(fixef(pois_glmm_ht)+VarCorr(pois_glmm_ht)$obs[1]/2)
  #Bias_glmm_ht[isim]=(pred_glmm_ht - True_Lam)/True_Lam
  
}

#plot(p_i,cooks.distance(pois_glm))

mean(Bias_reg)
mean(Bias_ht)
mean(Bias_mle)
mean(Bias_wtd)
mean(Bias_glmm)
mean(Bias_glmm_wtd)

sqrt(var(Bias_reg)/10000)
sqrt(var(Bias_wtd)/10000)
sqrt(var(Bias_ht)/10000)
sqrt(var(Bias_glmm)/10000)
sqrt(var(Bias_glmm_wtd)/10000)




library(ggplot2)
#scatter plot
Plot_df = data.frame("lambda_i"=exp(Pars[,1]),"o_i"=exp(Pars[,2]))
scatter_plot = ggplot(Plot_df)+geom_point(aes(x=lambda_i,y=o_i))+
  xlab(expression(paste("Poisson intensity (",lambda[i],")")))+
  ylab(expression(paste("Effort offset (",o[i],")")))+
  theme(text = element_text(size = 14))

png(filename="scatter_pois.png",width=6,height=6,units='in',res=600)
  scatter_plot
dev.off()


Plot_df = data.frame("Bias"=c(Bias_reg,Bias_glmm,Bias_ht,Bias_wtd),
                     "Model"=rep(c("GLM","GLMM","GLM-adj","GLM-wtd"),each=n_sim))
ggplot(Plot_df)+geom_boxplot(aes(y=Bias,color=Model))+
  scale_x_discrete(labels=c("GLM","GLM-adj","GLM-wtd","GLMM"))


#



