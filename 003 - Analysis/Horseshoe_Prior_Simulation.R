#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Horseshoe Prior Simulation  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Libraries
library(dplyr, warn.conflicts = F)
library(nimble)
# library(jagsUI) # If nimble is collapsing
library(crch) # For truncated T dist
library(LaplacesDemon) # For Half-Cauchy


# Plot range

x <- seq(-10, 10, length.out = 1e3)

## Normal horseshoe prior (Carvalho et al 2009) ----

### Using dhs and rhs from LaplaceDemon

## Lambdas (from halfcauchy)
N <-700

set.seed(1234)

lam <- abs(rcauchy(n = N,location = 0, scale = 1)) # all positive - could use rhalfcauchy(N, 1)

densities <- purrr::map(lam, ~dhs(x, lambda = .x, tau = 0.5))

dens.std <- Reduce(`+`, densities)/N

plot.df <- data.frame(x  = x, density = dens.std)

ggplot(plot.df, aes(x = x, y = density))+
  geom_line()

##