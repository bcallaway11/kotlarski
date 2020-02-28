#-----------------------------------------------------------------------------
# Tuning parameters
#
# In practice, these may need to be chosen carefully and in a
# theoretically valid way.  Here, I leave to the user to specify.
#-----------------------------------------------------------------------------

# set limits to evaluate integrals
a <- -3
b <- 3

#-----------------------------------------------------------------------------

# set number of points in between limits
nt <- 50

# index of grid points
i <- 0:(nt-1)

# space between grid points
dx <- (b-a)/nt

# grid points for the pdf
xgrid <- a + i * dx

# old
# dt  <- 2*pi/ (nt*dx)
# c <- -nt/2 * dt
# d <- nt/2 * dt
# tgrid <- c + i*dt

#-----------------------------------------------------------------------------
# create grid to evaluate characteristic functions

# Note: be careful here, if you are not careful with points to
# evaluate characteristic functions, can get things that are
# "close" to divide by 0s (this happens when grid is too "wide")

# bounds
Ta <- -3
Tb <- 3

# size of grid
ns  <- 100

# create grid to evaluate characteristic functions
tgrid <- seq(Ta,Tb,length.out=ns)

# space between grid points (relies on equal spacing)
dt <- tgrid[2] - tgrid[1]

# not sure if use this any more
tmax <- Tb

#-----------------------------------------------------------------------------




#-----------------------------------------------------------------------------
# Helper functions

#' @title cf2dens.inner
#' @description a helper function for computing the inverse Fourier transform
#'  of a characteristic function.
#'
#' @param cf a characteristic function (should be evaluated at t)
#' @param t particular value to evaluate characteristic function at
#' @param x a particular value for the pdf
#'
#' @return scalar imaginary number
cf2dens.inner <- function(t, cf, x) {
    exp( -(0+1i)*t*x ) *cf(t)
}

#' @title cf2dens
#' @description compute inverse Fourier transform of a characteristic function
#'  (i.e., turn a characteristic function back into a pdf)
#'
#' @param cf a characteristic function (will be evaluated at all values of
#'  tgrid)
#' @param tgrid vector of all possible values to evaluate characteristic
#'  function
#' @param xgrid vector of values to recover pdf of corresponding characteristic
#'  function cf
#'
#' @return vector of values of pdf for all values of xgrid
cf2dens <- function(cf, tgrid, xgrid) {
    1+1
    innermat <- pbapply::pbsapply(xgrid, function(x) {
        sapply(tgrid, function(t) {
            cf2dens.inner(t, cf=cf, x=x)
        })
    }, cl=4)
    dens.1 <- apply(innermat, 2, sum)
    dens <- (dens.1*dt)/(2*pi)
    return(Re(dens))
    1+1
}

#' @title cfX2
#' @description characteristic function with two variables
#'
#' @param t1 evaluate characteristic function at this value
#' @param t2 evaluate characteristic function at this value
#'
#' @return scalar value of characteristic function at t1 and t2
cfX2 <- function(t1, t2) {
    mean( exp( (0 + 1i)*((t1*X1) + (t2*X2)) ) )
}


#' @title d1cfX2
#' @description derivative of characteristic function with two variables
#'  with respect to first variable
#'
#' @param t1 evaluate derivative at this value
#' @param t2 evaluate derivative at this value
#'
#' @return scalar value of derivative of characteristic function at t1 and t2
d1cfX2 <- function(t1, t2) {
    mean( exp( (0 + 1i)*((t1*X1) + (t2*X2)) ) * (0+1i)*X1 )
}

#' @title kot.inner
#' @description kotlarski theorem result for single value of s;
#'  that is, estimate the characteristic function of the true X
#'  with two measures X_1 = X + e_1, X_2 = X + e_2.  This function is
#'  called by a wrapper that obtains the kotlarski result for a
#'  vector of values of s
#'
#' @param s a particular value to evaluate the charactersitic function at
#'
#' @return scalar value of characteristic function
kot.inner <- function(s) {
    #if (s <= 0) s <- -s
    # this involves calculating an integral from 0 to s
    # so calculation changes depending on if s is positive or negative
    if (s >= 0) {
        thist <- tgrid[ (tgrid >= 0 & tgrid <= s)]
    } else { # handle negative s
        thist <- tgrid[ (tgrid <= 0 & tgrid >= s)]
    }

    # calculate the integral
    inner1 <- sapply(thist, function(t) d1cfX2(0, t))
    inner2  <- sapply(thist, function(t) cfX2(0, t))
    inner  <- (inner1/inner2)*dt
    if (s >= 0) {
        return(exp(sum(inner)))
    } else {
        return(exp(-sum(inner)))
    }  
}

#' @title kotlarski
#' @description Estimate the characteristic function of the true X
#'  with two measures X_1 = X + e_1, X_2 = X + e_2
#'
#' @param tvec a vector of values to estimate the characteristic function at
#'
#' @return vector of the values of the characteristic function of X
#'  for each value of tvec
kotlarski <- function(tvec) {
    sapply(tvec, function(t) kot.inner(t))
}


## case where distribution of measurement error is known
cfX1 <- function(t) {
    mean( exp( (0+1i)*t*X1) )
}

#' @title cfXX
#' @description distribution of X when e_1 is known to be ~N(0,1)
#'
#' @param t value to evaluate the characteristic function at (can be vector)
#'
#' @return vector of the values of the characterstic function of X
#'  for each value of t
cfXX <- function(t) {
    cfX1(t) / cfN(t)
}



#' @title cfN
#' @description characteristic function of N(\mu,\sigma^2)
#'
#' @param s value to evaluate the characteristic function at (can be vector)
#' @param mu mean
#' @param sig2 variance
#'
#' @return vector of the values of the characteristic function for each
#'  value of s
cfN <- function(s, mu=0, sig2=1) {
    exp( (0+1i)*s*mu - 0.5*sig2*s^2 )
}

#' @title cfE
#' @description characteristic function of exponential(\lambda)
#'
#' @param s value to evaluate the characteristic function at (can be vector)
#' @param lam rate parameter
#'
#' @return vector of the values of the characteristic function for each
#'  value of s
cfE <- function(s, lam=1) {
    ( 1 - (0+1i)*s*lam^(-1) )^(-1)
}

#' @title cfX
#' @description Estimate the characteristic function of some random variable X
#'
#' @param tvec vector of values to evaluate the characteristic function at
#'
#' @return vector of the values of the characteristic function of X for
#'  each value of tvec
cfX <- function(tvec) {
    sapply(tvec, function(t) {
        mean( exp( (0+1i)*t*X) )
    })
}





