# clean workspace
rm(list=ls())

source("all-functions.R")

library(deSolve)
library(mvtnorm)
library(MASS)
library(ggplot2)

# read the data
data <- read.table("data/data_acari_3h.txt")
head(data)

# transform them 
new.data <- new.transform.data(data)
head(new.data)

# plot
plot(data, xlab=expression(N[0]), expression(N[e]),main=expression(paste("observed data")))

# set numbers of predators and observation time. 
nP <- 1; # number of predators
T <- 3; # observation time (in hours)

res.mle.type.II.alt <- optim(c(1, 1), loglikelihood.Hollings.type.II.unbalanced, data = new.data, n.pred = nP, obs.time = T, control=list(fnscale=-1,maxit=1000),hessian=TRUE)
# mles
res.mle.type.II.alt$par
# value of the log-lik when maximised
res.mle.type.II.alt$val
# var-covar matrix
obs.Fisher.alt <- solve(-res.mle.type.II.alt$hessian)
# confidence internvals for theta1 (first column) and theta2 (second column)
rbind(res.mle.type.II.alt$par - 1.96*sqrt(diag(solve(-res.mle.type.II.alt$hessian))), res.mle.type.II.alt$par + 1.96*sqrt(diag(solve(-res.mle.type.II.alt$hessian))))
# it is clear that the Normal approximation to the MLEs is not working since the lower bound of the confidence interval
# is negative. However, the confidence interval for the handling time (note that this is theta2 is not too bad.) 


# do mcmc using uniform priors on theta1 [0,4] and [0,2] for theta2 
res.P2.type.II.alt <- mcmc.Hollings.type.II.unif.priors.indep(data = new.data, n.pred=1, lambda.upper = c(4,2), iter = 10000, t = T, 
                                                              init = res.mle.type.II.alt$par,theta.hat = res.mle.type.II.alt$par, var.cov.hat = obs.Fisher.alt)

# trace plots
burn.in <- 0.10*dim(res.P2.type.II.alt)[1]
par(mfrow=c(2,1))
plot(res.P2.type.II.alt[-c(1:burn.in),1],type='l', main="theta_1", ylab="");
abline(h=res.mle.type.II.alt$par[1],col=2)
plot(res.P2.type.II.alt[-c(1:burn.in),2],type='l', main="theta_2", ylab="");
abline(h=res.mle.type.II.alt$par[2],col=2)

# marginal posterior distributions with priors
hist(res.P2.type.II.alt[-c(1:burn.in),1],prob=TRUE, xlab="", main="theta_1"); 
abline(v=res.mle.type.II.alt$par[1],col=2, lwd=2)
hist(res.P2.type.II.alt[-c(1:burn.in),2],prob=TRUE, xlab="", main="theta_2"); 
abline(v=res.mle.type.II.alt$par[2],col=2)

################################################
# ------ 24h experiment -- fit SSS equation
################################################

# read the data
data <- read.table("data/data_acari_24h.txt")
new.data <- transform.data(data)

# plot
par(mfrow=c(1,1)); plot(data, xlab=expression(N[0]), expression(N[e]),
main=expression(paste("observed data")))

# parameters
np <- 1
T <- 24

# maximum likelihood
res.mle.SSS <- optim(c(0.2, res.mle.type.II.alt$par[2],0.1), loglikelihood.SSS, data = new.data, n.pred = nP, obs.time = T, control=list(fnscale=-1, maxit=1000),hessian=FALSE)

# fit a gamma distribution to the posterior distribution of handling time from the Type II model 
prior.theta2.SSS <- mom.gamma(res.P2.type.II.alt[,2])

# have a look at how goog the approximation is
hist(res.P2.type.II.alt[,2] , prob=TRUE, xlab="handling time", main="prior on handling time")
lines(density(rgamma(10^5, prior.theta2.SSS$shape, prior.theta2.SSS$rate)),col=2,lwd=2)

# maximise the posterior density (ie no MCMC)
res.MAP.SSS <- optim(0.2*res.mle.SSS$par, log.post.den.SSS, data = new.data, n.pred = nP, 
                     nu=c(1, prior.theta2.SSS$shape, 1), lambda=c(0.1, prior.theta2.SSS$rate, 0.1), 
                     obs.time = T, control=list(fnscale=-1, maxit=1000), hessian=TRUE);

# MAP estimates
res.MAP.SSS
res.MAP.SSS$par
obs.Fisher.SSS <- solve(-res.MAP.SSS$hessian)

# plot observed versus fitted
x.values <- seq(0, 20, len = 100)
y.values.MAP.SSS <- rep(NA, length(x.values))
for ( i in 1:length(y.values.MAP.SSS)) y.values.MAP.SSS[i] <- x.values[i] - solve.ode.SSS(res.MAP.SSS$par, x.values[i], n.pred=nP, t=T)

# do the plot
plot(data, xlab=expression(N[0]), expression(N[e]),main=expression(paste("observed data"))) 
lines(x.values, y.values.MAP.SSS,col = 6, lwd=3)

# APPROXIMATE 95% credible intervals
res.MAP.SSS$par
rbind(res.MAP.SSS$par - 1.96*sqrt(diag(obs.Fisher.SSS)), res.MAP.SSS$par + 1.96*sqrt(diag(obs.Fisher.SSS)))


# mcmc on SSS model
# DO NOT RUN THIS # res.SSS.mcmc <- mcmc.SSS(data = new.data, n.pred = 1, nu = c(1,prior.theta2.SSS$shape, 1),lambda = c(0.1, prior.theta2.SSS$rate, 0.1), t = T, sigma = c(0.25, 0.25, 0.25),iter = 30000,init = res.MAP.SSS$par)

res.SSS.mcmc <- mcmc.SSS.indep(data = new.data, n.pred = 1, nu = c(1, prior.theta2.SSS$shape, 1), lambda=c(0.1, prior.theta2.SSS$rate, 0.1), iter = 20000, t = T, init = res.MAP.SSS$par, theta.hat = res.MAP.SSS$par, var.cov.hat = obs.Fisher.SSS)

# trace plots
par(mfrow=c(2,2))
plot(res.SSS.mcmc[,1],type='l')
abline(h=res.MAP.SSS$par[1],col=2, lwd=2)
plot(res.SSS.mcmc[,2],type='l')
abline(h=res.MAP.SSS$par[2],col=2, lwd=2)
plot(res.SSS.mcmc[,3],type='l')
abline(h=res.MAP.SSS$par[3],col=2, lwd=2)

# histograms
par(mfrow=c(2,2))
hist(res.SSS.mcmc[,1],prob=TRUE)
hist(res.SSS.mcmc[,2],prob=TRUE); lines(density(rgamma(10^5, prior.theta2.SSS$rate, prior.theta2.SSS$rate)))
hist(res.SSS.mcmc[,3],prob=TRUE)

# 2.5% 50% and 97.5% quantiles
quantile(res.SSS.mcmc[,1], c(0.025, 0.5, 0.975))
quantile(res.SSS.mcmc[,2], c(0.025, 0.5, 0.975))
quantile(res.SSS.mcmc[,3], c(0.025, 0.5, 0.975))
