# this version adds in discrete option to simulation function

#### Weibull functions ####
weib.1cp.nll <- function(tau, gamma = 2, dta) {
  G <- sum(dta$censor)
  T <- (1/gamma) * sum(dta$time^gamma)
  G1 <- sum(dta$time < tau)
  G2 <- G - G1
  T2 <- (1/gamma) * sum((dta$time^gamma - tau^gamma) * 
                          (dta$time >= tau))
  T1 <- T - T2
  # using MLEs for thetas
  ll <- G1 * log(G1/T1) + G2 * log(G2/T2) - G1 - 
    G2 + (gamma - 1) * sum(dta$censor * log(dta$time))
  return(-ll)
}

weib.1cp.LRT<-function(tau, gamma = 2, dta) {
  G <- sum(dta$censor)
  T <- (1/gamma) * sum(dta$time ^ gamma)
  G1 <- sum(dta$time < tau)
  G2 <- G - G1
  T2 <- (1/gamma) * sum((dta$time^gamma - tau^gamma) * (dta$time >= tau))
  T1 <- T - T2
  lrt <- G1 * log((G1 / T1) * (T / G)) + G2 * log((G2 / T2) * (T / G))
  return(-lrt) #need to use minimum so that we can compare to optim fn
}

weib.1cp.BIC<-function(tau , gamma = 2, dta) {
  3 * log(dim(dta)[1]) + 2 * weib.1cp.nll(tau, gamma, dta)
}

weib.2cp.nll <- function(tau, gamma = 2, dta) {
  G <- sum(dta$censor)
  T <- (1/gamma) * sum(dta$time^gamma)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum(dta$censor * (dta$time >= tau[2]))
  T1 <- (1/gamma) * sum((dta$time^gamma) * (dta$time < tau[1]) + 
                          (tau[1]^gamma) * (tau[1] <= dta$time))
  T2 <- (1/gamma) * sum((dta$time^gamma - tau[1]^gamma) * 
                          (tau[1] <= dta$time) * (dta$time < tau[2]) + 
                          (tau[2]^gamma - tau[1]^gamma) * (dta$time >= 
                                                             tau[2]))
  T3 <- (1/gamma) * sum((dta$time^gamma - tau[2]^gamma) * 
                          (dta$time >= tau[2]))
  ll <- G1 * log(G1/T1) + G2 * log(G2/T2) + G3 * 
    log(G3/T3) - G1 - G2 - G3 + (gamma - 1) * sum(dta$censor * 
                                                    log(dta$time))
  return(-ll)
}

weib.2cp.LRT <- function(tau, gamma = 2, dta) {
  G <- sum(dta$censor)
  T <- (1/gamma) * sum(dta$time^gamma)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum(dta$censor * (dta$time >= tau[2]))
  T1 <- (1/gamma) * sum((dta$time^gamma) * (dta$time < tau[1]) +
                    (tau[1]^gamma) * (tau[1] <= dta$time))
  T2 <- (1/gamma) * sum((dta$time^gamma - tau[1]^gamma) * (tau[1] <= dta$time) * 
                        (dta$time < tau[2]) +
                      (tau[2]^gamma - tau[1]^gamma) * (dta$time >= tau[2]))
  T3 <- (1/gamma) * sum((dta$time^gamma - tau[2]^gamma) * (dta$time >= tau[2]))
  lrt <- G1 * log((G1 / T1) * (T / G)) + G2 * log((G2 / T2) * (T / G)) +
          G3 * log((G3 / T3) * (T / G))
  return(-lrt)
}

weib.2cp.BIC <- function(tau, gamma = 2, dta) {
  5 * log(dim(dta)[1]) + 2 * weib.2cp.nll(tau, gamma, dta)
}

weib.3cp.nll <- function(tau, gamma = 2, dta) {
  G <- sum(dta$censor)
  T <- (1/gamma) * sum(dta$time^gamma)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum((tau[2] <= dta$time) * (dta$time < tau[3]))
  G4 <- sum((dta$time >= tau[3]) * dta$censor)
  T1 <- (1/gamma) * sum((dta$time^gamma) * (dta$time < 
                                              tau[1]) + (tau[1]^gamma) * (tau[1] <= dta$time))
  T2 <- (1/gamma) * sum((dta$time^gamma - tau[1]^gamma) * 
                          (tau[1] <= dta$time) * (dta$time < tau[2]) + 
                          (tau[2]^gamma - tau[1]^gamma) * (dta$time >= 
                                                             tau[2]))
  T3 <- (1/gamma) * sum((dta$time^gamma - tau[2]^gamma) * 
                          (tau[2] <= dta$time) * (dta$time < tau[3]) + 
                          (tau[3]^gamma - tau[2]^gamma) * (dta$time >= 
                                                             tau[3]))
  T4 <- (1/gamma) * sum((dta$time^gamma - tau[3]^gamma) * 
                          (dta$time >= tau[3]))
  # using MLEs for thetas
  ll <- G1 * log(G1/T1) + G2 * log(G2/T2) + G3 * 
    log(G3/T3) + G4 * log(G4/T4) - G1 - G2 - G3 - 
    G4 + (gamma - 1) * sum(dta$censor * log(dta$time))
  return(-ll)
}

weib.3cp.LRT<-function(tau, gamma = 2, dta) {
  G <- sum(dta$censor)
  T <- (1/gamma) * sum(dta$time^gamma)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum((tau[2] <= dta$time) * (dta$time < tau[3]))
  G4 <- sum((dta$time >= tau[3]) * dta$censor)
  T1 <- (1/gamma) * sum((dta$time^gamma) * (dta$time < tau[1]) + 
                          (tau[1]^gamma) * (tau[1] <= dta$time))
  T2 <- (1/gamma) * sum((dta$time^gamma - tau[1]^gamma) * 
                          (tau[1] <= dta$time) * (dta$time < tau[2]) + 
                          (tau[2]^gamma - tau[1]^gamma) * (dta$time >= tau[2]))
  T3 <- (1/gamma) * sum((dta$time^gamma - tau[2]^gamma) * (tau[2] <= dta$time) * (dta$time < tau[3]) + 
                          (tau[3]^gamma - tau[2]^gamma) * (dta$time >= tau[3]))
  T4 <- (1/gamma) * sum((dta$time^gamma - tau[3]^gamma) * (dta$time >= tau[3]))
  lrt <- G1 * log((G1/T1) * (T/G)) + G2 * log((G2/T2) * (T/G)) + 
          G3 * log((G3/T3) * (T/G)) + G4 * log((G4/T4) * (T/G))
  return(-lrt)
} 

weib.3cp.BIC <- function(tau, gamma = 2, dta) {
  7 * log(dim(dta)[1]) + 2 * weib.3cp.nll(tau, gamma, dta)
}

weib.4cp.nll <- function(tau, gamma = 2, dta) {
  G <- sum(dta$censor)
  T <- (1/gamma) * sum(dta$time^gamma)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum((tau[2] <= dta$time) * (dta$time < tau[3]))
  G4 <- sum((tau[3] <= dta$time) * (dta$time < tau[4]))
  G5 <- sum((dta$time >= tau[4]) * dta$censor)
  T1 <- (1/gamma) * sum((dta$time^gamma) * (dta$time < tau[1]) + 
                          (tau[1]^gamma) * (tau[1] <= dta$time))
  T2 <- (1/gamma) * sum((dta$time^gamma - tau[1]^gamma) * 
                          (tau[1] <= dta$time) * (dta$time < tau[2]) + 
                          (tau[2]^gamma - tau[1]^gamma) * (dta$time >= 
                                                             tau[2]))
  T3 <- (1/gamma) * sum((dta$time^gamma - tau[2]^gamma) * 
                          (tau[2] <= dta$time) * (dta$time < tau[3]) + 
                          (tau[3]^gamma - tau[2]^gamma) * (dta$time >= 
                                                             tau[3]))
  T4 <- (1/gamma) * sum((dta$time^gamma - tau[3]^gamma) * 
                          (tau[3] <= dta$time) * (dta$time < tau[4]) + 
                          (tau[4]^gamma - tau[3]^gamma) * (dta$time >= 
                                                             tau[4]))
  T5 <- (1/gamma) * sum((dta$time^gamma - tau[4]^gamma) * 
                          (dta$time >= tau[4]))
  ll <- G1 * log(G1/T1) + G2 * log(G2/T2) + G3 *log(G3/T3) + 
    G4 * log(G4/T4) + G5 * log(G5/T5) - 
    G1 - G2 - G3 - G4 - G5 + (gamma - 1) * sum(dta$censor * 
                                                 log(dta$time))
  return(-ll)
}

weib.4cp.LRT<-function(tau, gamma = 2) {
  G <- sum(dta$censor)
  T <- (1/gamma) * sum(dta$time^gamma)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum((tau[2] <= dta$time) * (dta$time < tau[3]))
  G4 <- sum((tau[3] <= dta$time) * (dta$time < tau[4]))
  G5 <- sum((dta$time >= tau[4]) * dta$censor)
  T1 <- (1/gamma) * sum((dta$time^gamma) * (dta$time < tau[1]) +
                          (tau[1]^gamma) * (tau[1] <= dta$time))
  T2 <- (1/gamma) * sum((dta$time^gamma - tau[1]^gamma) * (tau[1] <= dta$time) * (dta$time < tau[2]) +
                          (tau[2]^gamma - tau[1]^gamma) * (dta$time >= tau[2]))
  T3 <- (1/gamma) * sum((dta$time^gamma - tau[2]^gamma) * (tau[2] <= dta$time) * (dta$time < tau[3]) +
                          (tau[3]^gamma - tau[2]^gamma) * (dta$time >= tau[3]))
  T4 <- (1/gamma) * sum((dta$time^gamma - tau[3]^gamma) * (tau[3] <= dta$time) * (dta$time < tau[4]) +
                          (tau[4]^gamma - tau[3]^gamma) * (dta$time >= tau[4]))
  T5 <- (1/gamma) * sum((dta$time^gamma - tau[4]^gamma) * (dta$time >= tau[4]))
  lrt<- G1 * log((G1/T1) * (T/G)) + G2 * log((G2/T2) * (T/G)) +
        G3 * log((G3/T3) * (T/G)) + G4 * log((G4/T4) * (T/G)) +
        G5 * log((G5/T5) * (T/G))
  return(-lrt)
} 

weib.4cp.BIC <- function(tau, gamma = 2, dta) {
  9 * log(dim(dta)[1]) + 2 * weib.4cp.nll(tau, gamma, dta)
}

#### Exponential functions ####
exp.1cp.nll <- function(tau, dta) {
  t1.mle <- sum(dta$time < tau & dta$censor == 1) / 
    (sum(dta$time * (dta$time < tau) + tau * (dta$time >= tau)))
  t2.mle <- (sum(dta$censor) - sum(dta$time < tau & dta$censor == 1)) / 
    (sum((dta$time - tau) * (dta$time >= tau)))
  ll <- sum(dta$time < tau & dta$censor == 1) * log(t1.mle) +
    (sum(dta$censor) - sum(dta$time < tau & dta$censor == 1)) * log(t2.mle) - 
    sum(dta$censor)
  return(-ll)
}

exp.1cp.LRT <- function(tau, dta) {
  t1.mle <- sum(dta$time < tau & dta$censor == 1) / 
    (sum(dta$time * (dta$time < tau) + tau * (dta$time >= tau)))
  t2.mle <- (sum(dta$censor) - sum(dta$time < tau & dta$censor == 1)) / 
    (sum((dta$time - tau) * (dta$time >= tau)))
  lrt <- sum(dta$time < tau & dta$censor == 1) * log(t1.mle) + 
    (sum(dta$censor) - sum(dta$time < tau & dta$censor == 1)) * log(t2.mle) -
    sum(dta$censor) * log(sum(dta$censor) / sum(dta$time))
  return(-lrt) #use minimum for optim fn
}

exp.1cp.BIC <- function(tau, dta){
  3 * log(dim(dta)[1]) + 2 * exp.1cp.nll(tau, dta)
}

exp.2cp.nll <- function(tau, dta){
  t1.mle <- sum(dta$time < tau[1] & dta$censor == 1) / 
    (sum(dta$time * (dta$time < tau[1]) + tau[1] * (dta$time >= tau[1])))
  t2.mle <- (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) / 
    (sum((dta$time - tau[1]) * (dta$time >= tau[1]) * (dta$time < tau[2]) +
           (tau[2] - tau[1]) * (dta$time >= tau[2])))
  t3.mle <- (sum(dta$censor) - sum(dta$time < tau[2] & dta$censor == 1)) / 
    (sum((dta$time - tau[2]) * (dta$time >= tau[2])))
  ll <- sum(dta$time < tau[1] & dta$censor == 1) * log(t1.mle) + 
    (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) * 
    log(t2.mle) + (sum(dta$censor) - sum(dta$time < tau[2] & dta$censor == 1)) * 
    log(t3.mle) - sum(dta$censor)
  return(-ll)
}

exp.2cp.LRT <- function(tau, dta) {
  t1.mle <- sum(dta$time < tau[1] & dta$censor == 1) / 
    (sum(dta$time * (dta$time < tau[1]) + tau[1] * (dta$time >= tau[1])))
  t2.mle <- (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) / 
    (sum((dta$time - tau[1]) * (dta$time >= tau[1]) * (dta$time < tau[2]) +
           (tau[2] - tau[1]) * (dta$time >= tau[2])))
  t3.mle <- (sum(dta$censor) - sum(dta$time < tau[2] & dta$censor == 1)) / 
    (sum((dta$time - tau[2]) * (dta$time >= tau[2])))
  lrt <- sum(dta$time < tau[1] & dta$censor == 1) * log(t1.mle) + 
    (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) *
    log(t2.mle)+ (sum(dta$censor) - sum(dta$time < tau[2] & dta$censor == 1)) * log(t3.mle) -
    sum(dta$censor) * log(sum(dta$censor) / sum(dta$time))
  return(-lrt)
}

exp.2cp.BIC <- function(tau, dta) {
  5 * log(dim(dta)[1]) + 2 * exp.2cp.nll(tau, dta)
}

exp.3cp.nll <- function(tau, dta) {
  t1.mle <- sum(dta$time < tau[1] & dta$censor == 1) / 
    sum(dta$time * (dta$time < tau[1]) + tau[1] * (dta$time >= tau[1]))
  t2.mle <- (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) /
    sum((dta$time - tau[1]) * (dta$time >= tau[1]) * (dta$time < tau[2]) +
           (tau[2] - tau[1]) * (dta$time >= tau[2]))
  t3.mle <- (sum(dta$time < tau[3] & dta$censor == 1) - sum(dta$time < tau[2] & dta$censor == 1)) /
    sum((dta$time - tau[2]) * (dta$time >= tau[2]) * (dta$time < tau[3]) +
           (tau[3] - tau[2]) * (dta$time >= tau[3]))
  t4.mle <- (sum(dta$censor) - sum(dta$time < tau[3] & dta$censor == 1)) /
    sum((dta$time - tau[3]) * (dta$time>=tau[3]))
  ll <- sum(dta$time < tau[1] & dta$censor == 1) * log(t1.mle) + 
    (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) * log(t2.mle) +
    (sum(dta$time < tau[3] & dta$censor == 1) - sum(dta$time < tau[2] & dta$censor == 1)) * log(t3.mle) +
    (sum(dta$censor) - sum(dta$time < tau[3] & dta$censor == 1)) * log(t4.mle) - sum(dta$censor)
  return(-ll)
}

exp.3cp.LRT <- function(tau, dta) {
  t1.mle <- sum(dta$time < tau[1] & dta$censor == 1) / 
    sum(dta$time * (dta$time < tau[1]) + tau[1] * (dta$time >= tau[1]))
  t2.mle <- (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) /
    sum((dta$time - tau[1]) * (dta$time >= tau[1]) * (dta$time < tau[2]) +
          (tau[2] - tau[1]) * (dta$time >= tau[2]))
  t3.mle <- (sum(dta$time < tau[3] & dta$censor == 1) - sum(dta$time < tau[2] & dta$censor == 1)) /
    sum((dta$time - tau[2]) * (dta$time >= tau[2]) * (dta$time < tau[3]) +
          (tau[3] - tau[2]) * (dta$time >= tau[3]))
  t4.mle <- (sum(dta$censor) - sum(dta$time < tau[3] & dta$censor == 1)) /
    sum((dta$time - tau[3]) * (dta$time>=tau[3]))
  lrt <- sum(dta$time < tau[1] & dta$censor == 1) * log(t1.mle) +
    (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) * log(t2.mle) +
    (sum(dta$time < tau[3] & dta$censor == 1) - sum(dta$time < tau[2] & dta$censor == 1)) * log(t3.mle) +
    (sum(dta$censor) - sum(dta$time < tau[3] & dta$censor == 1)) * log(t4.mle) -
    sum(dta$censor) * log(sum(dta$censor) / sum(dta$time))
  return(-lrt)
}

exp.3cp.BIC <- function(tau, dta) {
  7 * log(dim(dta)[1]) + 2 * exp.3cp.nll(tau, dta)
}

exp.4cp.nll <- function(tau, dta) {
  t1.mle <- sum(dta$time < tau[1] & dta$censor == 1) / 
    sum(dta$time * (dta$time < tau[1]) + tau[1] * (dta$time >= tau[1]))
  t2.mle <- (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) /
    sum((dta$time - tau[1]) * (dta$time >= tau[1]) * (dta$time < tau[2]) +
          (tau[2] - tau[1]) * (dta$time >= tau[2]))
  t3.mle <- (sum(dta$time < tau[3] & dta$censor == 1) - sum(dta$time < tau[2] & dta$censor == 1)) /
    sum((dta$time - tau[2]) * (dta$time >= tau[2]) * (dta$time < tau[3]) +
          (tau[3] - tau[2]) * (dta$time >= tau[3]))
  t4.mle <- (sum(dta$time < tau[4] & dta$censor == 1) - sum(dta$time < tau[3] & dta$censor == 1)) /
    sum((dta$time - tau[3]) * (dta$time >= tau[3]) * (dta$time < tau[4]) +
          (tau[4] - tau[3]) * (dta$time >= tau[4]))
  t5.mle <- (sum(dta$censor) - sum(dta$time < tau[4] & dta$censor == 1)) /
    sum((dta$time - tau[4]) * (dta$time >= tau[4]))
  ll <- sum(dta$time < tau[1] & dta$censor == 1) * log(t1.mle) +
    (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) * log(t2.mle) +
    (sum(dta$time < tau[3] & dta$censor == 1) - sum(dta$time < tau[2] & dta$censor == 1)) * log(t3.mle) +
    (sum(dta$time < tau[4] & dta$censor == 1) - sum(dta$time < tau[3] & dta$censor == 1)) * log(t4.mle) +
    (sum(dta$censor) - sum(dta$time < tau[4] & dta$censor == 1)) * log(t5.mle) - sum(dta$censor)
  return(-ll)
}

exp.4cp.LRT <- function(tau, dta){
  t1.mle <- sum(dta$time < tau[1] & dta$censor == 1) / 
    sum(dta$time * (dta$time < tau[1]) + tau[1] * (dta$time >= tau[1]))
  t2.mle <- (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) /
    sum((dta$time - tau[1]) * (dta$time >= tau[1]) * (dta$time < tau[2]) +
          (tau[2] - tau[1]) * (dta$time >= tau[2]))
  t3.mle <- (sum(dta$time < tau[3] & dta$censor == 1) - sum(dta$time < tau[2] & dta$censor == 1)) /
    sum((dta$time - tau[2]) * (dta$time >= tau[2]) * (dta$time < tau[3]) +
          (tau[3] - tau[2]) * (dta$time >= tau[3]))
  t4.mle <- (sum(dta$time < tau[4] & dta$censor == 1) - sum(dta$time < tau[3] & dta$censor == 1)) /
    sum((dta$time - tau[3]) * (dta$time >= tau[3]) * (dta$time < tau[4]) +
          (tau[4] - tau[3]) * (dta$time >= tau[4]))
  t5.mle <- (sum(dta$censor) - sum(dta$time < tau[4] & dta$censor == 1)) /
    sum((dta$time - tau[4]) * (dta$time >= tau[4]))
  lrt <- sum(dta$time < tau[1] & dta$censor == 1) * log(t1.mle) +
    (sum(dta$time < tau[2] & dta$censor == 1) - sum(dta$time < tau[1] & dta$censor == 1)) * log(t2.mle) +
    (sum(dta$time < tau[3] & dta$censor == 1) - sum(dta$time < tau[2] & dta$censor == 1)) * log(t3.mle) +
    (sum(dta$time < tau[4] & dta$censor == 1) - sum(dta$time < tau[3] & dta$censor == 1)) * log(t4.mle) +
    (sum(dta$censor) - sum(dta$time < tau[4] & dta$censor == 1)) * log(t5.mle) -
    sum(dta$censor) * log(sum(dta$censor) / sum(dta$time))
  return(-lrt)
}

#BIC
exp.4cp.BIC <- function(tau, dta) {
  9 * log(dim(dta)[1]) + 2 * exp.4cp.nll(tau, dta)
}

#### Functions for simulations ####

### New weibull function simulating from inv cdf that fixes error with no cps
weib_invcdf <- function(n, endtime, gamma, theta, tau = NA) {
  n <- as.numeric(n)
  x <- stats::rexp(n)
  if(is.na(tau[1])){
    cdfcp0 <- function(v){
      ((gamma / theta) * v) ^ (1 / gamma)
    }
    t <- cdfcp0(x)
  }
  
  if(is.na(tau[1]) == FALSE & length(tau) == 1){
    first <- (theta / gamma) * tau ^ gamma
    cdfcp1 <- function(v) {
      ifelse(v < first, ((gamma / theta[1]) * v) ^ (1 / gamma),
             ((gamma / theta[2]) * (v - first) + tau ^ gamma) ^ (1 / gamma))
    }
    t <- cdfcp1(x)
  }
  
  if(length(tau) == 2) {
    first <- (theta[1] / gamma) * tau[1] ^ gamma #set interval 1
    second <- first + (theta[2] / gamma) * (tau[2] ^ gamma - tau[1] ^ gamma) #set interval 2
    cdfcp2 <- function(v) {
      ifelse(v < first, ((gamma / theta[1]) * v) ^ (1 / gamma),
             ifelse(v < second,
                    ((gamma / theta[2]) * (v - first) + tau[1] ^ gamma) ^ (1 / gamma),
                    ((gamma / theta[3]) * (v - second) + tau[2] ^ gamma) ^ (1 / gamma)))
    }
    t <- cdfcp2(x)
  }
  
  if(length(tau) == 3) {
    first <- (theta[1] / gamma) * tau[1] ^ gamma
    second <- first + (theta[2] / gamma) * (tau[2] ^ gamma - tau[1] ^ gamma)
    third <- second + (theta[3] / gamma) * (tau[3] ^ gamma - tau[2] ^ gamma)
    cdfcp3 <- function(v) {
      ifelse(v < first, ((gamma / theta[1]) * v) ^ (1 / gamma),
             ifelse(v < second,
                    ((gamma / theta[2]) * (v - first) + tau[1] ^ gamma) ^ (1 / gamma),
                    ifelse(v < third,
                           ((gamma / theta[3]) * (v - second) + tau[2] ^ gamma) ^ (1 / gamma),
                           ((gamma / theta[4]) * (v - third) + tau[3] ^ gamma) ^ (1 / gamma))))
    }
    t <- cdfcp3(x)
  }
  
  if(length(tau) == 4) {
    first <- (theta[1] / gamma) * tau[1] ^ gamma
    second <- first + (theta[2] / gamma) * (tau[2] ^ gamma - tau[1] ^ gamma)
    third <- second + (theta[3] / gamma) * (tau[3] ^ gamma - tau[2] ^ gamma)
    fourth <- third + (theta[4] / gamma) * (tau[4] ^ gamma - tau[3] ^ gamma)
    cdfcp4 <- function(v) {
      ifelse(v < first, ((gamma / theta[1]) * v) ^ (1 / gamma),
             ifelse(v < second,
                    ((gamma / theta[2]) * (v - first) + tau[1] ^ gamma) ^ (1 / gamma),
                    ifelse(v < third,
                           ((gamma / theta[3]) * (v - second) + tau[2] ^ gamma) ^ (1 / gamma),
                           ifelse(v < fourth,
                                  ((gamma / theta[4]) * (v-third) + tau[3] ^ gamma) ^ (1 / gamma),
                                  ((gamma / theta[5]) * (v-fourth) + tau[4] ^ gamma) ^ (1 / gamma)))))
    }
    t <- cdfcp4(x)
  }
  endtime <- as.numeric(endtime)
  C <- rep(endtime, length(x)) #all censored at endtime
  time <- pmin(t, C)  #observed time is min of censored and true
  censor <- as.numeric(time != endtime) #if not endtime then dropout
  dta <- data.frame(time = time, censor = censor)
  return(dta)
}

exp_invcdf <- function(n, endtime, theta, tau = NA) {
  n <- as.integer(n)
  x <- stats::rexp(n)
  if (is.na(tau[1]) == TRUE) {
    t <- x/theta
  }
  if (is.na(tau[1]) == FALSE & length(tau) == 1) { # this is the one change to this function
    first <- theta[1] * tau
    cdfcp1 <- function(v) {
      ifelse(v < first, v/theta[1], ((v - first)/theta[2]) + 
               tau)
    }
    t <- cdfcp1(x)
  }
  if (length(tau) == 2) {
    first <- theta[1] * tau[1]
    second <- first + theta[2] * (tau[2] - tau[1])
    cdfcp2 <- function(v) {
      ifelse(v < first, v/theta[1], ifelse(v < second, 
                                           ((v - first)/theta[2]) + tau[1], ((v - second)/theta[3]) + 
                                             tau[2]))
    }
    t <- cdfcp2(x)
  }
  if (length(tau) == 3) {
    first <- theta[1] * tau[1]
    second <- first + theta[2] * (tau[2] - tau[1])
    third <- second + theta[3] * (tau[3] - tau[2])
    cdfcp3 <- function(v) {
      ifelse(v < first, v/theta[1], ifelse(v < second, 
                                           ((v - first)/theta[2]) + tau[1], ifelse(v < third, 
                                                                                   ((v - second)/theta[3]) + tau[2], ((v - third)/theta[4]) + 
                                                                                     tau[3])))
    }
    t <- cdfcp3(x)
  }
  if (length(tau) == 4) {
    first <- theta[1] * tau[1]
    second <- first + theta[2] * (tau[2] - tau[1])
    third <- second + theta[3] * (tau[3] - tau[2])
    fourth <- third + theta[4] * (tau[4] - tau[3])
    cdfcp4 <- function(v) {
      ifelse(v < first, v/theta[1], ifelse(v < second, 
                                           ((v - first)/theta[2]) + tau[1], ifelse(v < third, ((v - second)/theta[3]) + tau[2], 
                                                                                   ifelse(v < fourth, ((v - third)/theta[4]) + tau[3],
                                                                                          ((v - fourth)/theta[5]) + tau[4]))))
    }
    t <- cdfcp4(x)
  }
  C <- rep(endtime, length(x))
  time <- pmin(t, C)
  censor <- as.numeric(time != endtime)
  dta <- data.frame(time = time, censor = censor)
  return(dta)
}