## Functions for applying test for multiple
## change-points within the exponential change-point
## hazard function


#### overall functions ####
nll <- function(par, dta) {
  par * sum(dta$time) - sum(dta$censor) * log(par)
}

# X: number of events that occur before time t
X <- function(t) {
  sum(dta$time < t & dta$censor == 1)
}

##### functions for one change-point #####

# MLE log-likelihood
nll.1cp <- function(tau) {
  t1.mle <- X(tau)/(sum(dta$time * (dta$time < tau) + 
                          tau * (dta$time >= tau)))
  t2.mle <- (sum(dta$censor) - X(tau))/(sum((dta$time - 
                                               tau) * (dta$time >= tau)))
  ll <- X(tau) * log(t1.mle) + (sum(dta$censor) - 
                                  X(tau)) * log(t2.mle) - sum(dta$censor)
  return(-ll)
}

# LRT
LRT.1cp <- function(tau) {
  t1.mle <- X(tau)/(sum(dta$time * (dta$time < tau) + 
                          tau * (dta$time >= tau)))
  t2.mle <- (sum(dta$censor) - X(tau))/(sum((dta$time - 
                                               tau) * (dta$time >= tau)))
  lrt <- X(tau) * log(t1.mle) + (sum(dta$censor) - 
                                   X(tau)) * log(t2.mle) - sum(dta$censor) * log(sum(dta$censor)/sum(dta$time))
  return(-lrt)  #use minimum for optim fn
}


# BIC
BIC.1cp <- function(tau) {
  3 * log(n) + 2 * nll.1cp(tau)
}


##### functions for two change-points #####

# MLE log-likelihood
nll.2cp <- function(tau) {
  t1.mle <- X(tau[1])/(sum(dta$time * (dta$time < 
                                         tau[1]) + tau[1] * (dta$time >= tau[1])))
  t2.mle <- (X(tau[2]) - X(tau[1]))/(sum((dta$time - 
                                            tau[1]) * (dta$time >= tau[1]) * (dta$time < 
                                                                                tau[2]) + (tau[2] - tau[1]) * (dta$time >= 
                                                                                                                 tau[2])))
  t3.mle <- (sum(dta$censor) - X(tau[2]))/(sum((dta$time - 
                                                  tau[2]) * (dta$time >= tau[2])))
  ll <- X(tau[1]) * log(t1.mle) + (X(tau[2]) - X(tau[1])) * 
    log(t2.mle) + (sum(dta$censor) - X(tau[2])) * 
    log(t3.mle) - sum(dta$censor)
  return(-ll)
}

# LRT
LRT.2cp <- function(tau) {
  t1.mle <- X(tau[1])/(sum(dta$time * (dta$time < 
                                         tau[1]) + tau[1] * (dta$time >= tau[1])))
  t2.mle <- (X(tau[2]) - X(tau[1]))/(sum((dta$time - 
                                            tau[1]) * (dta$time >= tau[1]) * (dta$time < 
                                                                                tau[2]) + (tau[2] - tau[1]) * (dta$time >= 
                                                                                                                 tau[2])))
  t3.mle <- (sum(dta$censor) - X(tau[2]))/(sum((dta$time - 
                                                  tau[2]) * (dta$time >= tau[2])))
  lrt <- X(tau[1]) * log(t1.mle) + (X(tau[2]) - X(tau[1])) * 
    log(t2.mle) + (sum(dta$censor) - X(tau[2])) * 
    log(t3.mle) - sum(dta$censor) * log(sum(dta$censor)/sum(dta$time))
  return(-lrt)
}

# BIC
BIC.2cp <- function(tau) {
  5 * log(n) + 2 * nll.2cp(tau)
}

##### functions for three change-points #####

# MLE log-likelihood
nll.3cp <- function(tau) {
  t1.mle <- X(tau[1])/(sum(dta$time * (dta$time < 
                                         tau[1]) + tau[1] * (dta$time >= tau[1])))
  t2.mle <- (X(tau[2]) - X(tau[1]))/(sum((dta$time - 
                                            tau[1]) * (dta$time >= tau[1]) * (dta$time < 
                                                                                tau[2]) + (tau[2] - tau[1]) * (dta$time >= 
                                                                                                                 tau[2])))
  t3.mle <- (X(tau[3]) - X(tau[2]))/(sum((dta$time - 
                                            tau[2]) * (dta$time >= tau[2]) * (dta$time < 
                                                                                tau[3]) + (tau[3] - tau[2]) * (dta$time >= 
                                                                                                                 tau[3])))
  t4.mle <- (sum(dta$censor) - X(tau[3]))/sum((dta$time - 
                                                 tau[3]) * (dta$time >= tau[3]))
  ll <- X(tau[1]) * log(t1.mle) + (X(tau[2]) - X(tau[1])) * 
    log(t2.mle) + (X(tau[3]) - X(tau[2])) * log(t3.mle) + 
    (sum(dta$censor) - X(tau[3])) * log(t4.mle) - 
    sum(dta$censor)
  return(-ll)
}

# LRT
LRT.3cp <- function(tau) {
  t1.mle <- X(tau[1])/(sum(dta$time * (dta$time < 
                                         tau[1]) + tau[1] * (dta$time >= tau[1])))
  t2.mle <- (X(tau[2]) - X(tau[1]))/(sum((dta$time - 
                                            tau[1]) * (dta$time >= tau[1]) * (dta$time < 
                                                                                tau[2]) + (tau[2] - tau[1]) * (dta$time >= 
                                                                                                                 tau[2])))
  t3.mle <- (X(tau[3]) - X(tau[2]))/(sum((dta$time - 
                                            tau[2]) * (dta$time >= tau[2]) * (dta$time < 
                                                                                tau[3]) + (tau[3] - tau[2]) * (dta$time >= 
                                                                                                                 tau[3])))
  t4.mle <- (sum(dta$censor) - X(tau[3]))/sum((dta$time - 
                                                 tau[3]) * (dta$time >= tau[3]))
  lrt <- X(tau[1]) * log(t1.mle) + (X(tau[2]) - X(tau[1])) * 
    log(t2.mle) + (X(tau[3]) - X(tau[2])) * log(t3.mle) + 
    (sum(dta$censor) - X(tau[3])) * log(t4.mle) - 
    sum(dta$censor) * log(sum(dta$censor)/sum(dta$time))
  return(-lrt)
}

# BIC
BIC.3cp <- function(tau) {
  7 * log(n) + 2 * nll.3cp(tau)
}

##### NLL functions for four change-points #####

# MLE log-likelihood
nll.4cp <- function(tau) {
  t1.mle <- X(tau[1])/(sum(dta$time * (dta$time < 
                                         tau[1]) + tau[1] * (dta$time >= tau[1])))
  t2.mle <- (X(tau[2]) - X(tau[1]))/(sum((dta$time - 
                                            tau[1]) * (dta$time >= tau[1]) * (dta$time < 
                                                                                tau[2]) + (tau[2] - tau[1]) * (dta$time >= 
                                                                                                                 tau[2])))
  t3.mle <- (X(tau[3]) - X(tau[2]))/(sum((dta$time - 
                                            tau[2]) * (dta$time >= tau[2]) * (dta$time < 
                                                                                tau[3]) + (tau[3] - tau[2]) * (dta$time >= 
                                                                                                                 tau[3])))
  t4.mle <- (X(tau[4]) - X(tau[3]))/(sum((dta$time - 
                                            tau[3]) * (dta$time >= tau[3]) * (dta$time < 
                                                                                tau[4]) + (tau[4] - tau[3]) * (dta$time >= 
                                                                                                                 tau[4])))
  t5.mle <- (sum(dta$censor) - X(tau[4]))/sum((dta$time - 
                                                 tau[4]) * (dta$time >= tau[4]))
  ll <- X(tau[1]) * log(t1.mle) + (X(tau[2]) - X(tau[1])) * 
    log(t2.mle) + (X(tau[3]) - X(tau[2])) * log(t3.mle) + 
    (X(tau[4]) - X(tau[3])) * log(t4.mle) + (sum(dta$censor) - 
                                               X(tau[4])) * log(t5.mle) - sum(dta$censor)
  return(-ll)
}

#LRT
LRT.4cp <- function(tau) {
  t1.mle <- X(tau[1])/(sum(dta$time * (dta$time < 
                                         tau[1]) + tau[1] * (dta$time >= tau[1])))
  t2.mle <- (X(tau[2]) - X(tau[1]))/(sum((dta$time - 
                                            tau[1]) * (dta$time >= tau[1]) * (dta$time < 
                                                                                tau[2]) + (tau[2] - tau[1]) * (dta$time >= 
                                                                                                                 tau[2])))
  t3.mle <- (X(tau[3]) - X(tau[2]))/(sum((dta$time - 
                                            tau[2]) * (dta$time >= tau[2]) * (dta$time < 
                                                                                tau[3]) + (tau[3] - tau[2]) * (dta$time >= 
                                                                                                                 tau[3])))
  t4.mle <- (X(tau[4]) - X(tau[3]))/(sum((dta$time - 
                                            tau[3]) * (dta$time >= tau[3]) * (dta$time < 
                                                                                tau[4]) + (tau[4] - tau[3]) * (dta$time >= 
                                                                                                                 tau[4])))
  t5.mle <- (sum(dta$censor) - X(tau[4]))/sum((dta$time - 
                                                 tau[4]) * (dta$time >= tau[4]))
  lrt <- X(tau[1]) * log(t1.mle) + (X(tau[2]) - X(tau[1])) * 
    log(t2.mle) + (X(tau[3]) - X(tau[2])) * log(t3.mle) + 
    (X(tau[4]) - X(tau[3])) * log(t4.mle) + (sum(dta$censor) - 
                                               X(tau[4])) * log(t5.mle) - sum(dta$censor) * 
    log(sum(dta$censor)/sum(dta$time))
  return(-lrt)
}

# BIC
BIC.4cp <- function(tau) {
  9 * log(n) + 2 * nll.4cp(tau)
}