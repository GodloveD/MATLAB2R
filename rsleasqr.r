rsleasqr <- function(resid, drdp, p, dpr, dpa, ponoff, stol, niter, erra, ...){
# rsleasqr <- function(resid, drdp, p, dpr, dpa, ponoff, stol, niter, erra){
  np    <- length(p)
  p     <- matrix(p, np)
  pbest <- p
  dpr   <- matrix(dpr, np)
  dpa   <- matrix(dpa, np)
  zcol  <- matrix(0, np, 1)
  zrow  <- matrix(0,  1, np)
  
  if (length(ponoff) == 0){
    ponoff <- matrix(1, np, 1)
  }
  
  fixp <- (dpr==0 & dpa==0) | ponoff==0 
  sc <- matrix(0, np, np)
  id <- diag(1, np)
  epstab <- c(1e-10, 0.001, 0.1, 1, 10, 100)
  
  vargin <- list(...)
  out    <-do.call(resid, c(list(p), vargin))
  sbest  <- out$ss
  rbest  <- out$r
  fbest  <- out$ff
  
  if(sbest < 0) return(0)
  nr <- length(rbest)
  
  for (iter in 1:max(c(1, niter))) {
    pprev <- pbest
    prt   <- do.call(drdp, c(list(p), vargin))
    qq    <- t(prt) %*% prt
    sprev <- sbest
    sgoal <- (1 - stol)*sprev
    r     <- rbest
    g     <- t(prt) %*% matrix(r, length(r)) 
    
    for (jj in 1:np) {
      
      if (fixp[jj] | qq[jj,jj]==0){
        qq[,jj]   <- 0 
        qq[jj,]   <- 0
        qq[jj,jj] <- 1
        g[jj]     <- 0
      } else {
        sc[jj,jj] <- 1/sqrt(qq[jj,jj])
      }
    }
    
    qq <- sc %*% qq %*% sc
    g  <- sc %*% g
    
    if (niter==0) break
    
    for (epsl in epstab) {
      p   <- -sweep(sc %*% (solve(qq+epsl*id,g)), 2, pprev)
      out <-do.call(resid, c(list(p), vargin))
      ss  <- out$ss
      rr  <- out$r
      ff  <- out$f
      
      if (ss >= 0) {
        if (ss < sbest) {
          sbest <- ss
          pbest <- p
          rbest <- rr
          fbest <- ff
        }
        
        if (ss <= sgoal) break
      }
    }
    
    if (sbest > sgoal | sbest==0){
      cnvrg <- TRUE
      break
    }
  }
  
  if (!erra) return(0)
  qi  <- solve((qq + 1e-14) %*% id)
  sef <- sqrt(sbest / max(c(1, nr-np + sum(fixp))))
  ze  <- matrix(0, length(fixp))
  ze[!fixp] <- 1
  if (length(qq)==1){
    se      <- sef * sc * sqrt(qi) * ze
    dep     <- ze
  } else {
    se  <- sef * diag(sc) * sqrt(diag(qi) * ze)
    dep <- sqrt(diag(qi) * ze * diag(qq)) 
  }
  qi < qi * sef^2
  for (jj in 1:np) {
    if (fixp[jj]) {
      qi[jj,jj] <- 0    
    }
  }
  
  output <- list(pbest=pbest,
                 se=se,
                 dep=dep,
                 rbest=rbest,
                 fbest=fbest,
                 sbest=sbest,
                 iter=iter,
                 cnvrg=cnvrg,
                 qi=qi) 
  return(output)              
}