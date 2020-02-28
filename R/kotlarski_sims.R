#-----------------------------------------------------------------------------
# Some simulations to check that everything works
#-----------------------------------------------------------------------------

# load the code
source("~/Dropbox/kotlarski/R/kotlarski.R")

n <- 5000
X <- rnorm(n)
e1 <- rnorm(n)
e2 <- rnorm(n)
X1 <- X + e1
X2 <- X + e2


# run the code to produce the pdf of x
dd <- cf2dens(kotlarski, tgrid, xgrid)

# plot the estimated pdf
plot(xgrid, dd)

# compare to true pdf
curve(dnorm(x), add=TRUE)        
