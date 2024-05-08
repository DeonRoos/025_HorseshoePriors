# Script to simulate data for simple logistic regression - w/ HS priors


# Sim data ----------------------------------------------------------------
# 2 real covariates
# 2 red herring covariates

N <- 100
K <- 1

set.seed(666)

x1 <- rnorm(N, 0, 1)
x2 <- rnorm(N, 0, 1)
# Red herring variables
x3 <- rnorm(N, 0, 1)
x4 <- rnorm(N, 0, 1)

b0 <- -1
b1 <- 0.08
b2 <- -0.05

p <- plogis(b0 + b1 * x1 + b2 * x2)
y <- rbinom(N, K, prob = p)

# Nimble data prep --------------------------------------------------------

# Data and constants
my_data <- list(y = y)
my_constants <- list(N = N,
                     x1 = x1,
                     x2 = x2,
                     x3 = x3,
                     x4 = x4)

# Parameters to monitor
parameters.to.save <- c("b0", "b1", "b2", "b3", "b4")

# Sampler settings
n.iter <- 2000
n.burnin <- 200
n.chains <- 3


# Initial values ----------------------------------------------------------

# Uninformed flat priors
initial.values <- function(){list(
  b0 = runif(1, -10, 10),
  b1 = runif(1, -10, 10),
  b2 = runif(1, -10, 10),
  b3 = runif(1, -10, 10),
  b4 = runif(1, -10, 10)
)}


# Standard flat -10,10 ----------------------------------------------------
# Uniform -10, 10 priors
uni_mod <- nimbleCode({
  
  # Priors
  b0 ~ dunif(-10, 10)
  b1 ~ dunif(-10, 10)
  b2 ~ dunif(-10, 10)
  b3 ~ dunif(-10, 10)
  b4 ~ dunif(-10, 10)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

uni_out <- nimbleMCMC(code = uni_mod, 
                      constants = my_constants,
                      data = my_data,              
                      inits = initial.values,
                      monitors = parameters.to.save,
                      niter = n.iter,
                      nburnin = n.burnin, 
                      nchains = n.chains)

MCMCsummary(uni_out, digits = 2)

MCMCtrace(uni_out,
          pdf = FALSE)

MCMCplot(uni_out,
         ci = c(50, 89),
         HPD = TRUE)
# Wide flat -100,100 ----------------------------------------------------

# Uniform -10, 10 priors
wuni_mod <- nimbleCode({
  
  # Priors
  b0 ~ dunif(-100, 100)
  b1 ~ dunif(-100, 100)
  b2 ~ dunif(-100, 100)
  b3 ~ dunif(-100, 100)
  b4 ~ dunif(-100, 100)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

wuni_out <- nimbleMCMC(code = wuni_mod, 
                      constants = my_constants,
                      data = my_data,              
                      inits = initial.values,
                      monitors = parameters.to.save,
                      niter = n.iter,
                      nburnin = n.burnin, 
                      nchains = n.chains)

MCMCsummary(wuni_out, digits = 2)

MCMCtrace(wuni_out,
          pdf = FALSE)

# Wide Normal priors ------------------------------------------------------

# Wide Normal 0, 10 priors
wn_mod <- nimbleCode({
  
  # Priors
  b0 ~ dnorm(0, 10)
  b1 ~ dnorm(0, 10)
  b2 ~ dnorm(0, 10)
  b3 ~ dnorm(0, 10)
  b4 ~ dnorm(0, 10)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

wn_out <- nimbleMCMC(code = wn_mod, 
                     constants = my_constants,
                     data = my_data,              
                     inits = initial.values,
                     monitors = parameters.to.save,
                     niter = n.iter,
                     nburnin = n.burnin, 
                     nchains = n.chains)

MCMCsummary(wn_out, digits = 2)

MCMCtrace(wn_out,
          pdf = FALSE)

MCMCplot(wn_out,
         ci = c(50, 89),
         HPD = TRUE)

# Normal Normal priors ------------------------------------------------------

# Wide Normal 0, 10 priors
nn_mod <- nimbleCode({
  
  # Priors
  b0 ~ dnorm(0, 1)
  b1 ~ dnorm(0, 1)
  b2 ~ dnorm(0, 1)
  b3 ~ dnorm(0, 1)
  b4 ~ dnorm(0, 1)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

nn_out <- nimbleMCMC(code = nn_mod, 
                     constants = my_constants,
                     data = my_data,              
                     inits = initial.values,
                     monitors = parameters.to.save,
                     niter = n.iter,
                     nburnin = n.burnin, 
                     nchains = n.chains)

MCMCsummary(nn_out, digits = 2)

MCMCtrace(nn_out,
          pdf = FALSE)


# All together ------------------------------------------------------------

MCMCplot(nn_out,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wn_out,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(uni_out,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wuni_out,
         ci = c(50, 89),
         HPD = TRUE)


# Initial values ----------------------------------------------------------

# Narrow normal inits
initial.values <- function(){list(
  b0 = rnorm(1, 0, 1),
  b1 = rnorm(1, 0, 1),
  b2 = rnorm(1, 0, 1),
  b3 = rnorm(1, 0, 1),
  b4 = rnorm(1, 0, 1)
)}


# Standard flat -10,10 ----------------------------------------------------
# Uniform -10, 10 priors
uni_mod2 <- nimbleCode({
  
  # Priors
  b0 ~ dunif(-10, 10)
  b1 ~ dunif(-10, 10)
  b2 ~ dunif(-10, 10)
  b3 ~ dunif(-10, 10)
  b4 ~ dunif(-10, 10)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

uni_out2 <- nimbleMCMC(code = uni_mod2, 
                      constants = my_constants,
                      data = my_data,              
                      inits = initial.values,
                      monitors = parameters.to.save,
                      niter = n.iter,
                      nburnin = n.burnin, 
                      nchains = n.chains)

MCMCsummary(uni_out2, digits = 2)

MCMCtrace(uni_out2,
          pdf = FALSE)

# Wide flat -100,100 ----------------------------------------------------

# Uniform -10, 10 priors
wuni_mod2 <- nimbleCode({
  
  # Priors
  b0 ~ dunif(-100, 100)
  b1 ~ dunif(-100, 100)
  b2 ~ dunif(-100, 100)
  b3 ~ dunif(-100, 100)
  b4 ~ dunif(-100, 100)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

wuni_out2 <- nimbleMCMC(code = wuni_mod2, 
                       constants = my_constants,
                       data = my_data,              
                       inits = initial.values,
                       monitors = parameters.to.save,
                       niter = n.iter,
                       nburnin = n.burnin, 
                       nchains = n.chains)

MCMCsummary(wuni_out2, digits = 2)

MCMCtrace(wuni_out2,
          pdf = FALSE)

# Wide Normal priors ------------------------------------------------------

# Wide Normal 0, 10 priors
wn_mod2 <- nimbleCode({
  
  # Priors
  b0 ~ dnorm(0, 10)
  b1 ~ dnorm(0, 10)
  b2 ~ dnorm(0, 10)
  b3 ~ dnorm(0, 10)
  b4 ~ dnorm(0, 10)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

wn_out2 <- nimbleMCMC(code = wn_mod2, 
                     constants = my_constants,
                     data = my_data,              
                     inits = initial.values,
                     monitors = parameters.to.save,
                     niter = n.iter,
                     nburnin = n.burnin, 
                     nchains = n.chains)

MCMCsummary(wn_out2, digits = 2)

MCMCtrace(wn_out2,
          pdf = FALSE)

# Normal Normal priors ------------------------------------------------------

# Wide Normal 0, 10 priors
nn_mod2 <- nimbleCode({
  
  # Priors
  b0 ~ dnorm(0, 1)
  b1 ~ dnorm(0, 1)
  b2 ~ dnorm(0, 1)
  b3 ~ dnorm(0, 1)
  b4 ~ dnorm(0, 1)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

nn_out2 <- nimbleMCMC(code = nn_mod2, 
                     constants = my_constants,
                     data = my_data,              
                     inits = initial.values,
                     monitors = parameters.to.save,
                     niter = n.iter,
                     nburnin = n.burnin, 
                     nchains = n.chains)

MCMCsummary(nn_out2, digits = 2)

MCMCtrace(nn_out2,
          pdf = FALSE)


# All together ------------------------------------------------------------
par(mfrow = c(4,2))
MCMCplot(nn_out,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wn_out,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(uni_out,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wuni_out,
         ci = c(50, 89),
         HPD = TRUE)
# More "Informed" initial values
MCMCplot(nn_out2,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wn_out2,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(uni_out2,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wuni_out2,
         ci = c(50, 89),
         HPD = TRUE)


# Initial values ----------------------------------------------------------

# Narrow normal inits
initial.values <- function(){list(
  b0 = rnorm(1, 0, 10),
  b1 = rnorm(1, 0, 10),
  b2 = rnorm(1, 0, 10),
  b3 = rnorm(1, 0, 10),
  b4 = rnorm(1, 0, 10)
)}


# Standard flat -10,10 ----------------------------------------------------
# Uniform -10, 10 priors
uni_mod3 <- nimbleCode({
  
  # Priors
  b0 ~ dunif(-10, 10)
  b1 ~ dunif(-10, 10)
  b2 ~ dunif(-10, 10)
  b3 ~ dunif(-10, 10)
  b4 ~ dunif(-10, 10)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

uni_out3 <- nimbleMCMC(code = uni_mod3, 
                       constants = my_constants,
                       data = my_data,              
                       inits = initial.values,
                       monitors = parameters.to.save,
                       niter = n.iter,
                       nburnin = n.burnin, 
                       nchains = n.chains)

MCMCsummary(uni_out3, digits = 2)

MCMCtrace(uni_out3,
          pdf = FALSE)

# Wide flat -100,100 ----------------------------------------------------

# Uniform -10, 10 priors
wuni_mod3 <- nimbleCode({
  
  # Priors
  b0 ~ dunif(-100, 100)
  b1 ~ dunif(-100, 100)
  b2 ~ dunif(-100, 100)
  b3 ~ dunif(-100, 100)
  b4 ~ dunif(-100, 100)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

wuni_out3 <- nimbleMCMC(code = wuni_mod3, 
                        constants = my_constants,
                        data = my_data,              
                        inits = initial.values,
                        monitors = parameters.to.save,
                        niter = n.iter,
                        nburnin = n.burnin, 
                        nchains = n.chains)

MCMCsummary(wuni_out3, digits = 2)

MCMCtrace(wuni_out3,
          pdf = FALSE)

# Wide Normal priors ------------------------------------------------------

# Wide Normal 0, 10 priors
wn_mod3 <- nimbleCode({
  
  # Priors
  b0 ~ dnorm(0, 10)
  b1 ~ dnorm(0, 10)
  b2 ~ dnorm(0, 10)
  b3 ~ dnorm(0, 10)
  b4 ~ dnorm(0, 10)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

wn_out3 <- nimbleMCMC(code = wn_mod3, 
                      constants = my_constants,
                      data = my_data,              
                      inits = initial.values,
                      monitors = parameters.to.save,
                      niter = n.iter,
                      nburnin = n.burnin, 
                      nchains = n.chains)

MCMCsummary(wn_out3, digits = 2)

MCMCtrace(wn_out3,
          pdf = FALSE)

# Normal Normal priors ------------------------------------------------------

# Wide Normal 0, 10 priors
nn_mod3 <- nimbleCode({
  
  # Priors
  b0 ~ dnorm(0, 1)
  b1 ~ dnorm(0, 1)
  b2 ~ dnorm(0, 1)
  b3 ~ dnorm(0, 1)
  b4 ~ dnorm(0, 1)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

nn_out3 <- nimbleMCMC(code = nn_mod3, 
                      constants = my_constants,
                      data = my_data,              
                      inits = initial.values,
                      monitors = parameters.to.save,
                      niter = n.iter,
                      nburnin = n.burnin, 
                      nchains = n.chains)

MCMCsummary(nn_out3, digits = 2)

MCMCtrace(nn_out3,
          pdf = FALSE)

# Wide Normal 0, 10 priors
nn_mod3 <- nimbleCode({
  
  # Priors
  b0 ~ dnorm(0, 1)
  b1 ~ dnorm(0, 1)
  b2 ~ dnorm(0, 1)
  b3 ~ dnorm(0, 1)
  b4 ~ dnorm(0, 1)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

nn_out3 <- nimbleMCMC(code = nn_mod3, 
                      constants = my_constants,
                      data = my_data,              
                      inits = initial.values,
                      monitors = parameters.to.save,
                      niter = n.iter,
                      nburnin = n.burnin, 
                      nchains = n.chains)

MCMCsummary(nn_out3, digits = 2)

MCMCtrace(nn_out3,
          pdf = FALSE)

# Very Normal priors ------------------------------------------------------

# Very narrow Normal 0, 0.25 priors
vn_mod3 <- nimbleCode({
  
  # Priors
  b0 ~ dnorm(0, 0.25)
  b1 ~ dnorm(0, 0.25)
  b2 ~ dnorm(0, 0.25)
  b3 ~ dnorm(0, 0.25)
  b4 ~ dnorm(0, 0.25)
  
  # Likelihood
  for(i in 1:N) {
    # True occupancy state
    y[i] ~ dbern(prob = p[i])
    logit(p[i]) <- b0 + b1 * x1[i] + b2 * x2[i] + b3 * x3[i] + b4 * x4[i]
  }
})

vn_out3 <- nimbleMCMC(code = vn_mod3, 
                      constants = my_constants,
                      data = my_data,              
                      inits = initial.values,
                      monitors = parameters.to.save,
                      niter = n.iter,
                      nburnin = n.burnin, 
                      nchains = n.chains)

MCMCsummary(vn_out3, digits = 2)

MCMCtrace(vn_out3,
          pdf = FALSE)

# All together ------------------------------------------------------------
par(mfrow = c(3,4))
MCMCplot(nn_out,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wn_out,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(uni_out,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wuni_out,
         ci = c(50, 89),
         HPD = TRUE)
# More "Informed" initial values
MCMCplot(nn_out2,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wn_out2,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(uni_out2,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wuni_out2,
         ci = c(50, 89),
         HPD = TRUE)
# Less "Informed" initial values
MCMCplot(nn_out3,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wn_out3,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(uni_out3,
         ci = c(50, 89),
         HPD = TRUE)
MCMCplot(wuni_out3,
         ci = c(50, 89),
         HPD = TRUE)