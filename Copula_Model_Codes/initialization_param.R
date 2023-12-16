rm(list = ls(all = T))
library(markovchain)
library(quantreg)
library(copula)
library(ald)

y_dummy <- our_data[-1,2:6]
x_dummy <- our_data[-nrow(our_data),7:10]
our_data <- data.frame(t = 1:nrow(y_dummy),y,x)
  
tau_fixed <- 0.05
T <- nrow(our_data); K <- 2
hidden_chain <- sample(1:K,size = T,replace = T)
transition_prob <- as.matrix(markovchainFit(hidden_chain)$estimate[1:K,])
transition_prob

beta0_1 <- vector(mode = 'list',length = K)
sigma0_1 <- vector(mode = "list",length = K)

#CQHMM
for(i in 1:K){
  m1 = lm(BTC ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data[hidden_chain == i,])
  m2 = lm(ETH ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data[hidden_chain == i,])
  m3 = lm(USDT ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data[hidden_chain == i,])
  m4 = lm(BNB ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data[hidden_chain == i,])
  m5 = lm(XRP ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data[hidden_chain == i,])
  beta0_1[[i]] <- cbind(m1$coefficients,m2$coefficients,m3$coefficients,m4$coefficients,
                       m5$coefficients)
  sigma0_1[[i]] <- c(sum(m1$residuals^2)/m1$df.residual,sum(m2$residuals^2)/m2$df.residual,
                     sum(m3$residuals^2)/m3$df.residual,sum(m4$residuals^2)/m4$df.residual,
                     sum(m5$residuals^2)/m5$df.residual)
  
}

#CEHMM
beta0_2 <- vector(mode = 'list',length = K)
sigma0_2 <- vector(mode = "list",length = K)
for(i in 1:K){
  m1 = rq(BTC ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data[hidden_chain == i,],tau = 0.5)
  m2 = rq(ETH ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data[hidden_chain == i,],tau = 0.5)
  m3 = rq(USDT ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data[hidden_chain == i,],tau = 0.5)
  m4 = rq(BNB ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data[hidden_chain == i,],tau = 0.5)
  m5 = rq(XRP ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data[hidden_chain == i,],tau = 0.5)
  beta0_2[[i]] <- cbind(m1$coefficients,m2$coefficients,m3$coefficients,m4$coefficients,
                        m5$coefficients)
  sigma0_2[[i]] <- c(mean(abs(m1$residuals)),mean(abs(m2$residuals)),mean(abs(m3$residuals)),
                     mean(abs(m4$residuals)),mean(abs(m5$residuals)))
}

#Initial Correlation Matrix
cormat0 <- vector(mode = 'list',length = K)
for(i in 1:K){
  cormat0[[i]] <- cor(our_data[hidden_chain == i,2:6])
}

G <- function(y,x,betas,sigmas){
  dens_val <- NULL;k <- length(sigmas)
  for(i in 1:k){
   #myCop <- tCopula(param = 5,dim = 5,dispstr = "un", correlation = cormat0[[i]])
   myCop <- tCopula(df = 5,dim = 5,param = cormat0[[i]][lower.tri(cormat0[[i]])],dispstr = "un")
   myMvd <- mvdc(copula=myCop, margins=c("ALD","ALD","ALD","ALD","ALD"),
                 paramMargins=list(list(mu = sum(c(1,x)*betas[[i]][,1]), sigma = sigmas[[i]][1],p = tau_fixed),
                                   list(mu = sum(c(1,x)*betas[[i]][,2]), sigma = sigmas[[i]][2],p = tau_fixed), 
                                   list(mu = sum(c(1,x)*betas[[i]][,3]), sigma = sigmas[[i]][3],p = tau_fixed),
                                   list(mu = sum(c(1,x)*betas[[i]][,4]), sigma = sigmas[[i]][4],p = tau_fixed),
                                   list(mu = sum(c(1,x)*betas[[i]][,5]), sigma = sigmas[[i]][5],p = tau_fixed)))
    dens_val <- c(dens_val,dMvdc(y,myMvd))
  }
  return(dens_val)
}

G(as.vector(as.matrix(our_data[4,2:6])),as.vector(as.matrix(our_data[3,7:10])),beta0_1,sigma0_1)



