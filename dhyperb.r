dhyperb <- function(p, ...){    # dhyperb computes
  vargin <- list(...)           # ...
  x    <- vargin[[1]]  
  vmax <- p[1]                  # ...  
  km   <- p[2]                  # dr(i)/dp(j), i.e.
  dr1  <- -x / (km+x)           # dr/dvmax ...
  dr2  <- vmax * x / ((km+x)^2) # and dr/dkm
  dr   <- matrix(c(dr1, dr2), length(dr2))
  return(dr)
}