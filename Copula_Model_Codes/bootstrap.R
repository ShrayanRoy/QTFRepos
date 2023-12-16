rm(list = ls(all = T))

G_sample <- function(y,x,betas,sigmas){
  
  myCop <- tCopula(df = 5,dim = 5,param = cormat0[[i]][lower.tri(cormat0[[i]])],
                     dispstr = "un")
  myMvd <- mvdc(copula=myCop, margins=c("ALD","ALD","ALD","ALD","ALD"),
                  paramMargins=list(list(mu = sum(c(1,x)*betas[,1]), 
                                         sigma = sigmas[1],p = tau_fixed),
                                    list(mu = sum(c(1,x)*betas[,2]), 
                                         sigma = sigmas[2],p = tau_fixed), 
                                    list(mu = sum(c(1,x)*betas[,3]), 
                                         sigma = sigmas[3],p = tau_fixed),
                                    list(mu = sum(c(1,x)*betas[,4]), 
                                         sigma = sigmas[4],p = tau_fixed),
                                    list(mu = sum(c(1,x)*betas[,5]), 
                                         sigma = sigmas[5],p = tau_fixed)))
  return(rMvdc(y,myMvd))
}



