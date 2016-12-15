rsexn <- function(tau, tt, Y, base){
  # rsexn: Fitting function for sums of exponentials.
  # Given exponential rate constants p, fit Y vs. t to sum of exponentials
  #    plus optional base level or base line.
  # Coefficients (C) of the exponentials and base function (level or line)
  #    are computed by separation of linear parameters.
  # -------------- Input Parameters ---------------------------------------
  # tau  = column of exponential time constants in increasing order.
  # tt   = column vector of time.
  # Y    = matrix of observed curve(s) y(t) in rows to be fitted by:
  #        sum_on_j [ c(q,j)*exp(-p(j)*t) ] + optional base level or
  #        function, where q is the index of the row of Y.
  # base = scalar base line flag,
  #      = 0 if no base level or line,
  #      = 1 for base level b1,
  #      = 2 for base function, currently b1 + b2/sqrt(t+1).
  # SSBEST_ & TAUBEST_ are global.  They must be initialized 
  #        by the user to initial sum of squares and 
  #        parameters before the fit begins.   They are also 
  #        output.   See below.
  # ----------------Output Parameters -------------------------------------
  # ss = sum of squares.
  # R  = matrix of residuals Y-F.
  # F  = matrix of computed curve(s) to be compared to Y.
  # C  = matrix of coefficients of exps. in the above fit,
  #      with coefficients of each curve in columns.
  # X  = matrix of exponentials.
  # SSBEST_ & TAUBEST_ are best sum of squares and parameters
  #      so far, updated whenever the sum of squares improves.
  #      This is a safeguard.   In case of an interrupted run,
  #      these results will not be destroyed.
  #---------------- Begin Program --------------------------------------- #
  tt <- t(tt)                                           # tt = row vector.       
  k  <- 1 / tau                                         # k = column.       
  
  # multiline switch case start                                                       
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~             
  X <- switch(base+1,                                   # Check base.            
              t(exp(tt %*% -k)),             # case 1   # If base = 0, matrix of exps.          
              #         
              rbind(t(exp(tt %*% -k)),       # case 2   # If base=1, exps. +             
                    matrix(0, 1, length(tt))),          #   base level.             
              #             
              rbind(t(exp(tt %*% -k)),       # case 3   # If base = 2, exps. +            
                    matrix(0, 1, length(tt)),           #   base level +             
                    t(1 / sqrt(tt + 1))))               #   base function.        
  #          
  if (is.null(X)) stop('Bad value of base') # otherwise # Error for other values of base.       
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                      
  # multiline switch case end                                                                              
  
  CC <- Y  %*% solve(X)                                 # Outputs.               
  FF <- CC %*% X                                        # Outputs.    
  R  <- Y-FF                                            # Outputs.
  
  ss <- sum(R^2)                                        # Sum of squares.    
  
  if (ss < SSBEST_){                                    # If ss has improved,    
    assign("SSBEST_",  ss,  envir = .GlobalEnv)         # save best ss     
    assign("TAUBEST_", tau, envir = .GlobalEnv)         # and tau.         
  } 
  
  output <- list(CC=CC, FF=FF, R=R)                     # key value pairs
  return(output)                                        # actually return them
}