#####################################################################################
#####################################################################################
# Define the ODE models
#####################################################################################
# note that in all these models:

# t denotes time
# x is the prey population
# theta is the parameter _vector_
# pred is the predator population.

ode.Hollings.type.II.alt <- function(t, x, theta, pred) {
  f <- - theta[1]*x*pred/(1+theta[1]*theta[2]*x);
  list(f)
}


### ODE for SSS equation model
ode.SSS <- function(t, x, theta, pred) {
  part1 <- 1 + theta[1]*x*(theta[2] + theta[3])
  part2 <- 1 + theta[1]*x*( 2*(theta[2]+theta[3]) + theta[1]*x*(theta[2]-theta[3])^2 )
  f <- - pred*(part1 - sqrt(part2))/(2*theta[1]*theta[2]*theta[3]*x)
  list(f)
}

# Holling's type II -- need to solve the ODE numerically.
solve.ode.Hollings.type.II.alt <- function(theta, N0, n.pred, t) {
  
  Hol.t <- seq(0, t, len = 5) # define the time range in which the ODE will be solved
  res <- ode(N0, Hol.t, ode.Hollings.type.II.alt, theta, pred=n.pred)[,-1] # sole the ODE
  res[length(Hol.t)] # only look at the final end point since this is what the data tells us. 
  
}

# SSS
solve.ode.SSS <- function(theta, N0, n.pred, t) {
  
  BD.t <- seq(0, t, len = 5)
  res <- ode(N0, BD.t, ode.SSS, theta, pred = n.pred)[,-1]
  res[length(BD.t)]
  
}



# This functions transforms the original dataset into a dataset where the computations are 
# minimized since we have repeated measurements essentially. 

transform.data <- function(data) {
  
  #first look for unique values in the first column
  N0.values <- unique(data[,1])
  
  # find the number of repetitions. --- add one for the first column which will have the unique values
  n.col <- length(which(N0.values[1]==data[,1])) + 1
  
  # create a table which will store the new data.
  new.data <- matrix(NA, nrow=length(N0.values), ncol=n.col)
  
  # assign to first column the initial values
  new.data[,1] <- N0.values
  
  # fill in the matrix
  for (i in 1:length(N0.values)) {
    new.data[i,2:ncol(new.data)] <- data[which(data[,1] == N0.values[i]),2]
  }
  new.data
}

new.transform.data <- function(data) {
  
  #first look for unique values in the first column
  N0.values <- unique(data[,1])
  
  # find the number of repetitions. --- add one for the first column which will have the unique values
  #n.col <- length(which(N0.values[1]==data[,1]))
  
  new.data.list <- vector('list', length(N0.values))
  
  # fill in the list
  for (i in 1:length(N0.values)) {
    new.data.list[[i]] <- data[which(data[,1] == N0.values[i]),2]
  }
  
  return(list("NO"=N0.values, "Ne"=new.data.list))
}



####################################################################################################

                    ####### log-likelihoods and posterior densities functions ########

####################################################################################################


# Holling's type II
loglikelihood.Hollings.type.II.alt <- function(theta, data, n.pred, obs.time) {
  
  loglik <- 0
  
  # create a vector of the probabilities which will 
  # correspond to the elements of the first columns
  
  prob.vec <- (data[,1] - sapply(data[,1], function(x) solve.ode.Hollings.type.II.alt(theta, N0=x, n.pred = n.pred, obs.time)))/data[,1]
  
  # check if any of the probabilities are > 1 or < 0 and if they are exit
  if ( sum(prob.vec > 0) == length(prob.vec) && sum(prob.vec <1) == length(prob.vec) && sum(is.na(prob.vec))==0  ) {
    
    for (i in 1:length(prob.vec)) {
      n.eaten.vec <- data[i,-1]
      loglik <- loglik + sum(sapply(n.eaten.vec, function(x) dbinom(x, data[i,1], prob.vec[i], log=TRUE)))
    }
    
  }
  else {
    loglik <- -Inf
  }
  return(loglik)
}

# Holling's type II
loglikelihood.Hollings.type.II.unbalanced <- function(theta, data, n.pred, obs.time) {
  
  loglik <- 0
  
  # create a vector of the probabilities which will 
  # correspond to the elements of the first columns
  
  prob.vec <- (data[[1]] - sapply(data[[1]], function(x) solve.ode.Hollings.type.II.alt(theta, N0=x, n.pred = n.pred, obs.time)))/data[[1]]
  #print(prob.vec)
  # check if any of the probabilities are > 1 or < 0 and if they are exit
  if ( sum(prob.vec > 0) == length(prob.vec) && sum(prob.vec <1) == length(prob.vec) && sum(is.na(prob.vec))==0  ) {
    
    for (i in 1:length(prob.vec)) {
      n.eaten.vec <- data[[2]][i]
      loglik <- loglik + sum(sapply(n.eaten.vec, function(x) dbinom(x, data[[1]][i], prob.vec[i], log=TRUE)))
    }
    
  }
  else {
    loglik <- -Inf
  }
  return(loglik)
}


### SSS
loglikelihood.SSS <- function(theta, data, n.pred, obs.time) {
  
  if (sum(theta > 0)==length(theta)){
    
    
    loglik <- 0
    
    # create a vector of the probabilities which will 
    # correspond to the elements of the first columns
    
    prob.vec <- (data[,1] - sapply(data[,1], function(x) solve.ode.SSS(theta, N0=x, n.pred = n.pred, obs.time)))/data[,1]
    
    # check if any of the probabilities are > 1 or < 0 and if they are exit
    if ( sum(prob.vec > 0) == length(prob.vec) && sum(prob.vec <1) == length(prob.vec) && sum(is.na(prob.vec))==0) {
      
      for (i in 1:length(prob.vec)) {
        n.eaten.vec <- data[i,-1]
        loglik <- loglik + sum(sapply(n.eaten.vec, function(x) dbinom(x, data[i,1], prob.vec[i], log=TRUE)))
      }
      
    }
    else {
      loglik <- -Inf
    }
    
  }
  else {
    loglik <- -Inf
  }
  return(loglik)
}

# log posterior densities
log.post.den.SSS <- function(theta, data, nu, lambda, n.pred, obs.time) {
  
  loglik <- 0
  
  # create a vector of the probabilities which will 
  # correspond to the elements of the first columns
  
  prob.vec <- (data[,1] - sapply(data[,1], function(x) solve.ode.SSS(theta, N0=x, n.pred = n.pred, obs.time)))/data[,1]
  
  # check if any of the probabilities are > 1 or < 0 and if they are exit
  if ( sum(prob.vec > 0) == length(prob.vec) && sum(prob.vec <1) == length(prob.vec) && sum(is.na(prob.vec))==0) {
    
    for (i in 1:length(prob.vec)) {
      n.eaten.vec <- data[i,-1]
      loglik <- loglik + sum(sapply(n.eaten.vec, function(x) dbinom(x, data[i,1], prob.vec[i], log=TRUE)))
    }
    
  }
  else {
    loglik <- -Inf
  }
  
  loglik +  dgamma(theta[1], nu[1], lambda[1], log=TRUE) + dgamma(theta[2], nu[2], lambda[2], log=TRUE) + dgamma(theta[3], nu[3], lambda[3], log=TRUE)
  
}




####################################################################################################

                              ####### MCMC FUNCTIONS ########

####################################################################################################


mcmc.SSS.indep <- function(data, n.pred = 1, nu, lambda, iter, t, init, theta.hat, var.cov.hat) {
  
  # vector to store the results
  res <- matrix(NA, nrow=iter, ncol=4);
  
  # initial values for the parameters.
  theta.cur <- init;
  
  theta1.cur <- theta.cur[1]
  theta2.cur <- theta.cur[2]
  theta3.cur <- theta.cur[3]
  
  # current value of the loglik
  loglik.cur <- loglikelihood.SSS(theta.cur, data, n.pred, t);
  log.prior.cur <- dgamma(theta1.cur, nu[1], lambda[1],log=TRUE) + dgamma(theta2.cur, nu[2], lambda[2],log=TRUE) + dgamma(theta3.cur, nu[3], lambda[3],log=TRUE)
  
  # store output
  res[1,] <- c(theta.cur, loglik.cur);
  
  
  # mcmc iter
  for (i in 2:iter) {
    
    # propose theta from the Normal Approximation to the MLE
    theta.can <- rmvnorm(1, mean=theta.hat, sigma=var.cov.hat)
    
    # check that the proposed values are combatible with the priors.
    if (sum(theta.can > 0) == length(theta.can)) {
      
      loglik.can <- loglikelihood.SSS(theta.can, data, n.pred, t);
      log.prior.can <- dgamma(theta.can[1], nu[1], lambda[1],log=TRUE) + dgamma(theta.can[2], nu[2], lambda[2],log=TRUE) + dgamma(theta.can[3], nu[3], lambda[3],log=TRUE)
      
      log.q.num <- dmvnorm(theta.cur, mean=theta.hat, sigma=var.cov.hat, log=TRUE)
      log.q.denom <- dmvnorm(theta.can, mean=theta.hat, sigma=var.cov.hat, log=TRUE)
      
      log.q.ratio <- log.q.num - log.q.denom
      
      u <- runif(1)
      if (log(u) < loglik.can - loglik.cur + log.prior.can - log.prior.cur +  log.q.ratio) {
        
        loglik.cur <- loglik.can
        theta.cur <- theta.can
        log.prior.cur <- log.prior.can
        
      }
      
    }
    
    # store the sample
    res[i,] <- c(theta.cur, loglik.cur);
  }
  
  return(res)    
}


mcmc.Hollings.type.II.unif.priors.indep <- function(data, n.pred = 1, lambda.upper, iter, t, init, theta.hat, var.cov.hat) {
  
  # vector to store the results
  res <- matrix(NA, nrow=iter, ncol=3);
  
  # initial values for the parameters.
  theta.cur <- init;
  
  theta1.cur <- theta.cur[1]
  theta2.cur <- theta.cur[2]
  
  # current value of the loglik
  loglik.cur <- loglikelihood.Hollings.type.II.unbalanced(theta.cur, data, n.pred, t);
  print(loglik.cur)
  
  # store output
  res[1,] <- c(theta.cur, loglik.cur);
  
  
  # mcmc iter
  for (i in 2:iter) {
    
    # propose theta from the Normal Approximation to the MLE
    theta.can <- rmvnorm(1, mean=theta.hat, sigma=var.cov.hat)
    
    # check that the proposed values are combatible with the priors.
    if (theta.can[1] < lambda.upper[1] && theta.can[1] > 0 && theta.can[2] < lambda.upper[2] && theta.can[2] > 0) {
      
      loglik.can <- loglikelihood.Hollings.type.II.unbalanced(theta.can, data, n.pred, t);
      
      log.q.num <- dmvnorm(theta.cur, mean=theta.hat, sigma=var.cov.hat, log=TRUE)
      log.q.denom <- dmvnorm(theta.can, mean=theta.hat, sigma=var.cov.hat, log=TRUE)
      
      log.q.ratio <- log.q.num - log.q.denom
      
      u <- runif(1)
      if (log(u) < loglik.can - loglik.cur + log.q.ratio) {
        
        loglik.cur <- loglik.can
        theta.cur <- theta.can
        
      }
      
    }
    
    # store the sample
    res[i,] <- c(theta.cur, loglik.cur);
  }
  
  return(res)    
}


mom.gamma <- function(x) {
  mu = mean(x)
  mu.sq = mu^2
  sigma.sq = var(x)
  
  list("shape"= mu.sq/sigma.sq, "rate"=mu/sigma.sq)
  
}
