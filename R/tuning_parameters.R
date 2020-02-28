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
