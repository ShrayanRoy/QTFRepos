rm(list = ls(all = T))

new_gamma_h <- results$gamma
new_m <- results$gamma[1,]
new_hidden_var <- results$state
new_Q <- Reduce(`+`,results$epsilon)
new_Q <- new_Q/rowSums(new_Q)
K <- nrow(new_Q)

new_betas <- vector(mode = 'list',length = K)
new_sigmas <- vector(mode = "list",length = K)

for(i in 1:K){
  m1 = rq(BTC ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data,
          tau = 0.5,weights = new_gamma_h[,i])
  m2 = rq(ETH ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data,
          tau = 0.5,weights = new_gamma_h[,i])
  m3 = rq(USDT ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data,
          tau = 0.5,weights = new_gamma_h[,i])
  m4 = rq(BNB ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data,
          tau = 0.5,weights = new_gamma_h[,i])
  m5 = rq(XRP ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data,
          tau = 0.5,weights = new_gamma_h[,i])
  new_betas[[i]] <- cbind(m1$coefficients,m2$coefficients,m3$coefficients,m4$coefficients,
                        m5$coefficients)
  new_sigmas[[i]] <- c(mean(abs(m1$residuals)),mean(abs(m2$residuals)),mean(abs(m3$residuals)),
                     mean(abs(m4$residuals)),mean(abs(m5$residuals)))
}


#Initial Correlation Matrix
new_cormat <- vector(mode = 'list',length = K)
for(i in 1:K){
  new_cormat[[i]] <- cor(our_data[new_hidden_var == i,2:6])
}

new_G <- function(y,x,betas,sigmas){
  dens_val <- NULL;k <- length(sigmas)
  for(i in 1:k){
    #myCop <- tCopula(param = 5,dim = 5,dispstr = "un", correlation = cormat0[[i]])
    myCop <- tCopula(df = 5,dim = 5,param = cormat0[[i]][lower.tri(new_cormat[[i]])],dispstr = "un")
    myMvd <- mvdc(copula=myCop, margins=c("ALD","ALD","ALD","ALD","ALD"),
                  paramMargins=list(list(mu = sum(c(1,x[1],x[2],x[3],x[4])*betas[[i]][,1]), sigma = sigmas[[i]][1],p = 0.5),
                                    list(mu = sum(c(1,x)*betas[[i]][,2]), sigma = sigmas[[i]][2],p = 0.5), 
                                    list(mu = sum(c(1,x)*betas[[i]][,3]), sigma = sigmas[[i]][3],p = 0.5),
                                    list(mu = sum(c(1,x)*betas[[i]][,4]), sigma = sigmas[[i]][4],p = 0.5),
                                    list(mu = sum(c(1,x)*betas[[i]][,5]), sigma = sigmas[[i]][5],p = 0.5)))
    dens_val <- c(dens_val,dMvdc(y,myMvd))
  }
  return(dens_val)
}
