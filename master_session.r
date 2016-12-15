rm(list = ls())

source_list = list(
  'dhyperb.r',
  'hyperb.r',
  'rsleasqr.r'
)

for (src in source_list) {
  source(file.path(getwd(), src))
}


# set up some dummy values
vmax <- 1;
km   <- 1;
p    <- t(c(vmax, km))
x    <- seq(0, 2, by=0.1)
y    <- matrix(0, 21, 1)

# generate some fake data
out <- hyperb(p, x, y)
ss  <- out$ss
r   <- out$r
ff  <- out$ff

# request for error analysis
niter <- 50
stol  <- 0.001
erra  <- 1;

# make some noise
y      <- ff+0.01 * rnorm(21)
dpr    <- matrix(1, 2, 1)
dpa    <- dpr
ponoff <- NULL

# run it
out <- rsleasqr('hyperb', 'dhyperb', p, dpr, dpa, ponoff, stol, niter, erra, x, y)

# display results
print('Parameters  Std err  Dependencies')
print(matrix(c(out$pbest, out$se, out$dep),2))
plot(x,y)
matplot(x,out$fbest,'l',add=T)
