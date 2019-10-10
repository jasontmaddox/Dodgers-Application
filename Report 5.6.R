library(R2jags); library(mcmcplots); library(coda); library(xtable)
y1 <- rnorm(25, 20, 5)
y2 <- rnorm(25, 21, 5)
y3 <- rnorm(25, 25, 8)

mu3ci <- function(sigmax, taumax){
  cod <- function(){
    for(i in 1:length(y)){
      y[i] ~ dnorm(mu[i], pres.sig)
      mu[i] ~ dnorm(beta[i], pres.tau)
      beta[i] ~ dnorm(mub[i], 1/(sigb[i]*sigb[i]))
    }
    sigma ~ dunif(0, sigmax)
    tau ~ dunif(0, taumax)
    pres.sig <- 25/(sigma*sigma)
    pres.tau <- 1/(tau*tau)
    p3 <- step(mu[3] - max(mu[1], mu[2]))
  }
  
  dat <- list(y = c(mean(y1), mean(y2), mean(y3)), mub = c(19, 21, 15),
              sigb = c(sqrt(5), sqrt(5), 5), sigmax = sigmax, taumax = taumax)
  
  bayesfit <- jags(dat, NULL, c('mu', 'beta', 'tau', 'sigma', 'p3'), cod,
                   n.chains = 4, n.iter = 210000, n.burnin = 10000, n.thin = 20)
  
  return(c(bayesfit$BUGSoutput$summary[7,1], bayesfit$BUGSoutput$summary[7, 3],
           bayesfit$BUGSoutput$summary[7, 7]))
}

mu3cissig <- sapply(1:20*10, mu3ci, taumax = 10)

mu3cissig <- as.data.frame(t(mu3cissig))

ggplot(mu3cissig, aes(x = 1:20*10, y = V1)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = V2, ymin = V3)) +
  scale_x_continuous(expression(B[sigma])) +
  scale_y_continuous(expression(mu[3]))

mu3cistau <- sapply(1:20+9, mu3ci, sigmax = 150)

mu3cistau <- as.data.frame(t(mu3cistau))

ggplot(mu3cistau, aes(x = 1:20+9, y = V1)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = V2, ymin = V3)) +
  scale_x_continuous(expression(B[tau])) +
  scale_y_continuous(expression(mu[3]))


cod <- function(){
  for(i in 1:length(y)){
    y[i] ~ dnorm(mu[i], pres.sig)
    mu[i] ~ dnorm(beta[i], pres.tau)
    beta[i] ~ dnorm(mub[i], 1/(sigb[i]*sigb[i]))
  }
  sigma ~ dunif(0, 150)
  tau ~ dunif(0, 25)
  pres.sig <- 25/(sigma*sigma)
  pres.tau <- 1/(tau*tau)
  p3 <- step(mu[3] - max(mu[1], mu[2]))
}

dat <- list(y = c(mean(y1), mean(y2), mean(y3)), mub = c(19, 21, 15),
            sigb = c(sqrt(5), sqrt(5), 5))

bayesfit <- jags(dat, NULL, c('mu', 'beta', 'tau', 'sigma', 'p3'), cod,
                 n.chains = 4, n.iter = 210000, n.burnin = 10000, n.thin = 20)


par(mfrow = c(5, 2))

#coda::traceplot(as.mcmc(bayesfit))
traceplot(bayesfit, varname = 'sigma', ask = F)
traceplot(bayesfit, varname = 'mu', ask = F)
autocorr.plot(as.mcmc(bayesfit), ask = F)

p3 <- cumsum(Reduce(c, as.list(as.mcmc(bayesfit)[,8])))/1:40000
par(mfrow = c(1, 1))
plot(1:40000*5, p3, type = 'l', ylim = c(0, 1), xlab = 'Iterations')

plot(as.mcmc(bayesfit), trace = F)
print(bayesfit)

