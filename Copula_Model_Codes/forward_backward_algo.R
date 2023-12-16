
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
  epsilon[[1]] = rbind(gamma_h[1,],gamma_h[1,],gamma_h[1,])
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
results = forward_backward_em_algorithm(y,x,transition_prob,G,c(1/3,1/3,1/3),beta0_1,sigma0_1)



