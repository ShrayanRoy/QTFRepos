rm(list = ls(all = T))  #removes all objects

library(markovchain)
library(quantreg)
library(copula)
library(ald)

tau_fixed <- 0.95
y_dummy <- our_data[-1,2:6]
x_dummy <- our_data[-nrow(our_data),7:10]
our_data <- data.frame(t = 1:nrow(y_dummy),y_dummy,x_dummy)
  
T <- nrow(our_data); K <- 3
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

# Modified forward-backward algorithm with parameter estimation
forward_backward_em_algorithm <- function(y,x,Q,G,m,betas,sigmas){
  
  T <- nrow(y)  # Number of time points
  K <- nrow(Q)    # Number of states
  prox <- rep(0.0001,K)
  
  # Initialization
  alpha <- matrix(0,nrow = T, ncol = K)
  beta <- matrix(0,nrow = T, ncol = K)
  
  # E-step (Forward-Backward)
  alpha[1, ] <- m * G(as.vector(as.matrix(y[1,])),as.vector(as.matrix(x[1,])),betas,sigmas) + 0.0001
  for (t in 2:T) {
    g <- G(as.vector(as.matrix(y[t,])),as.vector(as.matrix(x[t,])),betas,sigmas)
    alpha[t, ] <- apply(Q * crossprod(t(alpha[t - 1, ]),g + prox), 2, sum) + 0.0001
    alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
  }
  
  beta[T, ] <- rep(1, K)
  for (t in (T - 1):1) {
    g <- G(as.vector(as.matrix(y[t+1,])),as.vector(as.matrix(x[t+1,])),betas,sigmas)
    beta[t, ] <- apply(Q * crossprod(t(beta[t + 1, ]), (g + prox)), 1, sum) + 0.0001
    beta[t, ] <- beta[t, ] / sum(beta[t, ])
  }
  
  # Combine forward and backward probabilities to get smoothed probabilities
  gamma_h <- alpha * beta
  gamma_h <- gamma_h / rowSums(gamma_h)
  
  epsilon <- vector(mode='list', length = length(y))
  epsilon[[1]] = t(replicate(K,gamma_h[1,]))
  for(t in 2:T){
    g <- G(as.vector(as.matrix(y[t,])),as.vector(as.matrix(x[t,])),betas,sigmas)
    epsilon[[t]] = t(t(Q)*alpha[t-1,]*(g + prox)*beta[t,]) + 0.0001
  }
  
  epsilon <- lapply(epsilon,FUN = function(mat){mat/sum(mat)})
  hidden_var <- apply(gamma_h,1,FUN = function(arr){which.max(arr)})
  
  return(list(alpha = alpha, beta = beta, gamma = gamma_h,epsilon = epsilon,state = hidden_var))
}

y = our_data[,2:6]
x = our_data[,7:10]
Q = transition_prob
m = rep(1,K)/K
results = forward_backward_em_algorithm(y,x,Q,G,m,beta0_1,sigma0_1) 

for(epoch in 1:20){
  
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
            tau = tau_fixed,weights = new_gamma_h[,i])
    m2 = rq(ETH ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data,
            tau = tau_fixed,weights = new_gamma_h[,i])
    m3 = rq(USDT ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data,
            tau = tau_fixed,weights = new_gamma_h[,i])
    m4 = rq(BNB ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data,
            tau = tau_fixed,weights = new_gamma_h[,i])
    m5 = rq(XRP ~ X.GSPC + DX.Y.NYB + CL.F + GC.F ,data = our_data,
            tau = tau_fixed,weights = new_gamma_h[,i])
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
                    paramMargins=list(list(mu = sum(c(1,x)*betas[[i]][,1]), sigma = sigmas[[i]][1],p = tau_fixed),
                                      list(mu = sum(c(1,x)*betas[[i]][,2]), sigma = sigmas[[i]][2],p = tau_fixed), 
                                      list(mu = sum(c(1,x)*betas[[i]][,3]), sigma = sigmas[[i]][3],p = tau_fixed),
                                      list(mu = sum(c(1,x)*betas[[i]][,4]), sigma = sigmas[[i]][4],p = tau_fixed),
                                      list(mu = sum(c(1,x)*betas[[i]][,5]), sigma = sigmas[[i]][5],p = tau_fixed)))
      dens_val <- c(dens_val,dMvdc(y,myMvd))
    }
    return(dens_val)
  }
  
  results = forward_backward_em_algorithm(y,x,new_Q,new_G,new_m,new_betas,new_sigmas) 
}

new_betas
new_sigmas
new_cormat

table(results$state)


#========================================

graph_d <- bind_cols(Crypt_colwise[-1,],state = results$state)
ggpairs(graph_d,columns = 2:6, mapping = ggplot2::aes(color = as.factor(results$state))) + 
  theme_bw(14) + defined_theme
