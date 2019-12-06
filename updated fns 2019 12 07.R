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


## function for catching errors
is.error <- function(x) inherits(x, "try-error")



test_cp <- function(simset, cp, dist,
                    nsim = 1000,
                    start1 = c(1, simset[2]),
                    start2 = c(4, 16), 
                    start3 = c(4, 10, 16),
                    start4 = c(4, 8, 12, 16),
                    discrete = FALSE) {
  simset <- as.numeric(simset)
  if(dist == "exp") {
      out <- matrix(NA, nsim, 18)
      for(k in 1:nsim) {
        ## simulate data 
        if(cp == 0) {
          dta <- exp_invcdf(n = simset[1], endtime = simset[2], 
                             theta = simset[3])
        }
        if(cp == 1) {
          dta <- exp_invcdf(n = simset[1], endtime = simset[2], 
                             theta = simset[3:4], tau = simset[5])
        }
        if(cp == 2) {
          dta <- exp_invcdf(n = simset[1], endtime = simset[2], 
                             theta = simset[3:5], tau = simset[6:7])
        }
        if(cp == 3) {
          dta <- exp_invcdf(n = simset[1], endtime = simset[2], 
                             theta = simset[3:6], tau = simset[7:9])
        }
        
        if(discrete == TRUE) {
          dta$time <- ifelse(dta$time<1, 1, trunc(dta$time))
        }
        
        #LRT for one change-point
        lrt1cp <- optimize(exp.1cp.LRT, start1, dta = dta)
        out[k, 1:2] <- c(lrt1cp$objective, lrt1cp$minimum)
        out[k, 3] <- exp.1cp.BIC(tau = out[k, 2], dta = dta)
        
        #LRT for two change-points
        startpar <- start2
        lrt2cp <- try(optim(par = startpar, fn = exp.2cp.LRT, dta = dta), silent = TRUE)
        failed <- sum(is.error(lrt2cp))
        while(failed == 1) {
          startpar <- startpar + 0.001
          if(startpar[2] > simset[2]) {
            lrt2cp <- NA
            break
          }
          lrt2cp <- try(optim(par = startpar, fn = exp.2cp.LRT, dta = dta), silent = TRUE)
          failed <- sum(is.error(lrt2cp))
        }
        if(is.na(lrt2cp[1]) == TRUE) {
          out[k, 4:7] <- rep(NA, 4)
        }
        if(is.na(lrt2cp[1]) == FALSE) {
          out[k, 4:6] <- c(lrt2cp$value, lrt2cp$par)
          out[k, 7] <- exp.2cp.BIC(tau = out[k, 5:6], dta = dta)
        }
        
        #LRT for three change-points
        startpar2 <- start3
        lrt3cp <- try(optim(par = startpar2, fn = exp.3cp.LRT, dta = dta), silent = TRUE)
        failed<-sum(is.error(lrt3cp))
        while(failed == 1){
          startpar2 <- startpar2 + 0.001
          if(startpar2[3] > simset[2]) {
            lrt3cp <- NA
            break
          }
          lrt3cp <- try(optim(par = startpar2, fn = exp.3cp.LRT, dta = dta), silent = TRUE)
          failed<-sum(is.error(lrt3cp))
        }
        if(is.na(lrt3cp[1]) == TRUE) {
          out[k, 8:12] <- rep(NA, 5)
        }
        if(is.na(lrt3cp[1]) == FALSE) {
          out[k, 8:11] <- c(lrt3cp$value, lrt3cp$par)
          out[k, 12] <- exp.3cp.BIC(tau = out[k, 9:11], dta = dta)
        }
        
        #LRT for four change-points
        startpar3 <- start4
        lrt4cp <- try(optim(par = startpar3, fn = exp.4cp.LRT, dta = dta), silent = TRUE)
        failed <- sum(is.error(lrt4cp))
        while(failed == 1) {
          startpar3 <- startpar3 + c(-0.0001, 0.0002, -0.0002, 0.0001)
          if(startpar3[4] > simset[2]) {
            lrt4cp <- NA
            break
          }
          lrt4cp <- try(optim(par = startpar3, fn = exp.4cp.LRT, dta = dta), silent = TRUE)
          failed <- sum(is.error(lrt4cp))
        }
        if(is.na(lrt4cp[1]) == TRUE) {
          out[k, 13:18] <- rep(NA, 6)
        }
        if(is.na(lrt4cp[1]) == FALSE) {
          out[k, 13:17] <- c(lrt4cp$value, lrt4cp$par)
          out[k, 18] <- exp.4cp.BIC(tau = out[k, 14:17], dta = dta)
        }
      }
  }
  
  if(dist=="weib") {
    out <- matrix(NA, nsim, 18)
    for(k in 1:nsim) {
      ## simulate data 
      if(cp == 0) {
        dta <- weib_invcdf(n = simset[1], endtime = simset[2], gamma = simset[3], 
                           theta = simset[4])
      }
      if(cp == 1) {
        dta <- weib_invcdf(n = simset[1], endtime = simset[2], gamma = simset[3], 
                           theta = simset[4:5], tau = simset[6])
      }
      if(cp == 2) {
        dta <- weib_invcdf(n = simset[1], endtime = simset[2], gamma = simset[3], 
                           theta = simset[4:6], tau = simset[7:8])
      }
      if(cp == 3) {
        dta <- weib_invcdf(n = simset[1], endtime = simset[2], gamma = simset[3], 
                           theta = simset[4:7], tau = simset[8:10])
      }
      
      if(discrete == TRUE) {
        dta$time <- ifelse(dta$time<1, 1, trunc(dta$time))
      }
      
      #LRT for one change-point
      lrt1cp <- optimize(weib.1cp.LRT, start1, gamma = simset[3], dta = dta)
      out[k, 1:2] <- c(lrt1cp$objective, lrt1cp$minimum)
      out[k, 3] <- weib.1cp.BIC(tau = out[k, 2], gamma = simset[3], dta = dta)
      
      #LRT for two change-points
      startpar <- start2
      lrt2cp <- try(optim(par = startpar, fn = weib.2cp.LRT, gamma = simset[3], dta = dta), silent = TRUE)
      failed <- sum(is.error(lrt2cp))
      while(failed == 1) {
        startpar <- startpar + 0.001
        if(startpar[2] > simset[2]) {
          lrt2cp <- NA
          break
        }
        lrt2cp <- try(optim(par = startpar, fn = weib.2cp.LRT, gamma = simset[3], dta = dta), silent = TRUE)
        failed <- sum(is.error(lrt2cp))
      }
      if(is.na(lrt2cp[1]) == TRUE) {
        out[k, 4:7] <- rep(NA, 4)
      }
      if(is.na(lrt2cp[1]) == FALSE) {
        out[k, 4:6] <- c(lrt2cp$value, lrt2cp$par)
        out[k, 7] <- weib.2cp.BIC(tau = out[k, 5:6], gamma = simset[3], dta = dta)
      }
      
      #LRT for three change-points
      startpar2 <- start3
      lrt3cp <- try(optim(par = startpar2, fn = weib.3cp.LRT, gamma = simset[3], dta = dta), silent = TRUE)
      failed<-sum(is.error(lrt3cp))
      while(failed == 1){
        startpar2 <- startpar2 + 0.001
        if(startpar2[3] > simset[2]) {
          lrt3cp <- NA
          break
        }
        lrt3cp <- try(optim(par = startpar2, fn = weib.3cp.LRT, gamma = simset[3], dta = dta), silent = TRUE)
        failed<-sum(is.error(lrt3cp))
      }
      if(is.na(lrt3cp[1]) == TRUE) {
        out[k, 8:12] <- rep(NA, 5)
      }
      if(is.na(lrt3cp[1]) == FALSE) {
        out[k, 8:11] <- c(lrt3cp$value, lrt3cp$par)
        out[k, 12] <- weib.3cp.BIC(tau = out[k, 9:11], gamma = simset[3], dta = dta)
      }
      
      #LRT for four change-points
      startpar3 <- start4
      lrt4cp <- try(optim(par = startpar3, fn = weib.4cp.LRT, gamma = simset[3], dta = dta), silent = TRUE)
      failed <- sum(is.error(lrt4cp))
      while(failed == 1) {
        startpar3 <- startpar3 + c(-0.0001, 0.0002, -0.0002, 0.0001)
        if(startpar3[4] > simset[2]) {
          lrt4cp <- NA
          break
        }
        lrt4cp <- try(optim(par = startpar3, fn = weib.4cp.LRT, gamma = simset[3], dta = dta), silent = TRUE)
        failed <- sum(is.error(lrt4cp))
      }
      if(is.na(lrt4cp[1]) == TRUE) {
        out[k, 13:18] <- rep(NA, 6)
      }
      if(is.na(lrt4cp[1]) == FALSE) {
        out[k, 13:17] <- c(lrt4cp$value, lrt4cp$par)
        out[k, 18] <- weib.4cp.BIC(tau = out[k, 14:17], gamma = simset[3], dta = dta)
      }

    }
  }
  outdta <- data.frame(out)
  colnames(outdta) <- c("LRT1cp", "tau_1cp", "BIC1cp", 
                        "LRT2cp", "tau1_2cp", "tau2_2cp", "BIC2cp",
                        "LRT3cp", "tau1_3cp", "tau2_3cp", "tau3_3cp", "BIC3cp",
                        "LRT4cp", "tau1_4cp", "tau2_4cp", "tau3_4cp", "tau4_4cp", "BIC4cp")
  return(outdta)
}

power.error <- function(simset, cp, dist, nulldist,
                        nsim = 1000, alpha = 0.05,
                        start1 = c(1, simset[2]),
                        start2 = c(4, 16), 
                        start3 = c(4, 10, 16),
                        start4 = c(4, 8, 12, 16)) {
  output <- test_cp(simset, cp, dist, nsim, start1, start2, start3, start4)
  ### results...
  alpha2<-alpha/2
  alpha3<-alpha/4
  alpha4<-alpha/8
  if(simset[1] == 500) {
    probs1cp.alpha<-quantile(nulldist[[1]]$LRT1cp, probs=alpha)
    probs2cp.alpha<-quantile(nulldist[[1]]$LRT2cp, probs=alpha)
    probs3cp.alpha<-quantile(nulldist[[1]]$LRT3cp, probs=alpha)
    probs4cp.alpha<-quantile(nulldist[[1]]$LRT4cp, probs=alpha)
    
    probs1cp.spending<-quantile(nulldist[[1]]$LRT1cp, probs=alpha)
    probs2cp.spending<-quantile(nulldist[[1]]$LRT2cp, probs=alpha2)
    probs3cp.spending<-quantile(nulldist[[1]]$LRT3cp, probs=alpha3)
    probs4cp.spending<-quantile(nulldist[[1]]$LRT4cp, probs=alpha4)
    
    outstat <- results(output, cp, nsim)
    
    
  }
  if(simset[1] == 1000) {
    probs1cp.alpha<-quantile(nulldist[[2]]$LRT1cp, probs=alpha)
    probs2cp.alpha<-quantile(nulldist[[2]]$LRT2cp, probs=alpha)
    probs3cp.alpha<-quantile(nulldist[[2]]$LRT3cp, probs=alpha)
    probs4cp.alpha<-quantile(nulldist[[2]]$LRT4cp, probs=alpha)
    
    probs1cp.spending<-quantile(nulldist[[2]]$LRT1cp, probs=alpha)
    probs2cp.spending<-quantile(nulldist[[2]]$LRT2cp, probs=alpha2)
    probs3cp.spending<-quantile(nulldist[[2]]$LRT3cp, probs=alpha3)
    probs4cp.spending<-quantile(nulldist[[2]]$LRT4cp, probs=alpha4)
    
    outstat <- results(output, cp, nsim)
  }
  if(simset[1] != 500 & simset[1] != 1000) {
    return("No test statistic for this sample size")
  }
}

#### Functions for saving results ####

results <- function(resdta, cp, nsim=10000) {
  if(cp == 0) {
    nocp.BICres<-data.frame(nocp.BIC.1cp = resdta$BIC1cp,
                            nocp.BIC.2cp = resdta$BIC2cp, 
                            nocp.BIC.3cp = resdta$BIC3cp, 
                            nocp.BIC.4cp = resdta$BIC4cp)
    nocp.LRT.1cp.value <- resdta$LRT1cp
    nocp.LRT.2cp.value <- resdta$LRT2cp
    nocp.LRT.3cp.value <- resdta$LRT3cp
    nocp.LRT.4cp.value <- resdta$LRT4cp
    nocp.LRT.1cp.tau <- resdta$tau_1cp
    nocp.LRT.2cp.tau <- cbind(resdta$tau1_2cp, resdta$tau2_2cp)
    nocp.LRT.3cp.tau <- cbind(resdta$tau1_3cp, resdta$tau2_3cp, resdta$tau3_3cp)
    nocp.LRT.4cp.tau <- cbind(resdta$tau1_4cp, resdta$tau2_4cp, resdta$tau3_4cp, resdta$tau4_4cp)
    
    #### 0cp ####
    ## separate into testing datasets
    nocp.LRT.test1cp<-nocp.LRT.1cp.value[which((nocp.BICres$nocp.BIC.1cp<nocp.BICres$nocp.BIC.2cp & 
                                                  nocp.BICres$nocp.BIC.1cp<nocp.BICres$nocp.BIC.3cp & 
                                                  nocp.BICres$nocp.BIC.1cp<nocp.BICres$nocp.BIC.4cp) |
                                                 #2cp that should be 1 cp
                                                 (nocp.BICres$nocp.BIC.2cp<nocp.BICres$nocp.BIC.1cp & 
                                                    nocp.BICres$nocp.BIC.2cp<nocp.BICres$nocp.BIC.3cp & 
                                                    nocp.BICres$nocp.BIC.2cp<nocp.BICres$nocp.BIC.4cp &
                                                    (nocp.LRT.2cp.tau[,2]-nocp.LRT.2cp.tau[,1]<1 |
                                                       nocp.LRT.2cp.tau[,1]<2 | nocp.LRT.2cp.tau[,2]==20)) |
                                                 #3cp that should be 1cp
                                                 (nocp.BICres$nocp.BIC.3cp<nocp.BICres$nocp.BIC.1cp & 
                                                    nocp.BICres$nocp.BIC.3cp<nocp.BICres$nocp.BIC.2cp & 
                                                    nocp.BICres$nocp.BIC.3cp<nocp.BICres$nocp.BIC.4cp &
                                                    (nocp.LRT.3cp.tau[,2]-nocp.LRT.3cp.tau[,1]<1 |
                                                       nocp.LRT.3cp.tau[,3]-nocp.LRT.3cp.tau[,2]<1 | 
                                                       nocp.LRT.3cp.tau[,1]<2 | nocp.LRT.3cp.tau[,3]==20) &
                                                    (nocp.LRT.2cp.tau[,2]-nocp.LRT.2cp.tau[,1]<1 |
                                                       nocp.LRT.2cp.tau[,1]<2 | nocp.LRT.2cp.tau[,2]==20)) |
                                                 #4cp that should be 1 cp
                                                 (nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.1cp & 
                                                    nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.2cp & 
                                                    nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.3cp &
                                                    (nocp.LRT.4cp.tau[,4]-nocp.LRT.4cp.tau[,3]<1 | 
                                                       nocp.LRT.4cp.tau[,3]-nocp.LRT.4cp.tau[,2]<1 | 
                                                       nocp.LRT.4cp.tau[,2]-nocp.LRT.4cp.tau[,1]<1 | 
                                                       nocp.LRT.4cp.tau[,1]<2 | nocp.LRT.4cp.tau[,4]==20) &
                                                    (nocp.LRT.3cp.tau[,3]-nocp.LRT.3cp.tau[,2]<1 |
                                                       nocp.LRT.3cp.tau[,2]-nocp.LRT.3cp.tau[,1]<1 | 
                                                       nocp.LRT.3cp.tau[,1]<2 | nocp.LRT.3cp.tau[,3]==20) &
                                                    (nocp.LRT.2cp.tau[,2]-nocp.LRT.2cp.tau[,1]<1 |
                                                       nocp.LRT.2cp.tau[,1]<2 | nocp.LRT.2cp.tau[,2]==20)))]
    
    
    
    nocp.LRT.test2cp<-nocp.LRT.2cp.value[which((nocp.BICres$nocp.BIC.2cp<nocp.BICres$nocp.BIC.1cp & 
                                                  nocp.BICres$nocp.BIC.2cp<nocp.BICres$nocp.BIC.3cp & 
                                                  nocp.BICres$nocp.BIC.2cp<nocp.BICres$nocp.BIC.4cp &
                                                  nocp.LRT.2cp.tau[,2]-nocp.LRT.2cp.tau[,1]>1 & 
                                                  nocp.LRT.2cp.tau[,1]>=2 & nocp.LRT.2cp.tau[,2]<20)|
                                                 #3cp that should be 2cp 
                                                 (nocp.BICres$nocp.BIC.3cp<nocp.BICres$nocp.BIC.1cp & 
                                                    nocp.BICres$nocp.BIC.3cp<nocp.BICres$nocp.BIC.2cp & 
                                                    nocp.BICres$nocp.BIC.3cp<nocp.BICres$nocp.BIC.4cp &
                                                    (nocp.LRT.3cp.tau[,2]-nocp.LRT.3cp.tau[,1]<1 |
                                                       nocp.LRT.3cp.tau[,3]-nocp.LRT.3cp.tau[,2]<1 | 
                                                       nocp.LRT.3cp.tau[,1]<2 | nocp.LRT.3cp.tau[,3]==20) &
                                                    nocp.LRT.2cp.tau[,2]-nocp.LRT.2cp.tau[,1]>1 & 
                                                    nocp.LRT.2cp.tau[,1]>=2 & nocp.LRT.2cp.tau[,2]<20)|
                                                 #4cp that should be 2cp
                                                 (nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.1cp & 
                                                    nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.2cp & 
                                                    nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.3cp &
                                                    (nocp.LRT.4cp.tau[,4]-nocp.LRT.4cp.tau[,3]<1 | 
                                                       nocp.LRT.4cp.tau[,3]-nocp.LRT.4cp.tau[,2]<1 | 
                                                       nocp.LRT.4cp.tau[,2]-nocp.LRT.4cp.tau[,1]<1 | 
                                                       nocp.LRT.4cp.tau[,1]<2 | nocp.LRT.4cp.tau[,4]==20) &
                                                    (nocp.LRT.3cp.tau[,3]-nocp.LRT.3cp.tau[,2]<1 |
                                                       nocp.LRT.3cp.tau[,2]-nocp.LRT.3cp.tau[,1]<1 | 
                                                       nocp.LRT.3cp.tau[,1]<2 | nocp.LRT.3cp.tau[,3]==20) &
                                                    nocp.LRT.2cp.tau[,2]-nocp.LRT.2cp.tau[,1]>1 & 
                                                    nocp.LRT.2cp.tau[,1]>=2 & nocp.LRT.2cp.tau[,2]<20))]
    
    
    nocp.LRT.test3cp<-nocp.LRT.3cp.value[which((nocp.BICres$nocp.BIC.3cp<nocp.BICres$nocp.BIC.1cp & 
                                                  nocp.BICres$nocp.BIC.3cp<nocp.BICres$nocp.BIC.2cp & 
                                                  nocp.BICres$nocp.BIC.3cp<nocp.BICres$nocp.BIC.4cp &
                                                  nocp.LRT.3cp.tau[,2]-nocp.LRT.3cp.tau[,1]>1 &
                                                  nocp.LRT.3cp.tau[,3]-nocp.LRT.3cp.tau[,2]>1 & 
                                                  nocp.LRT.3cp.tau[,1] >=2 & nocp.LRT.3cp.tau[,3]<20)|
                                                 #4cp should be 3cp
                                                 (nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.1cp & 
                                                    nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.2cp & 
                                                    nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.3cp &
                                                    #not actually 4 change-points
                                                    (nocp.LRT.4cp.tau[,4]-nocp.LRT.4cp.tau[,3]<1 | 
                                                       nocp.LRT.4cp.tau[,3]-nocp.LRT.4cp.tau[,2]<1 | 
                                                       nocp.LRT.4cp.tau[,2]-nocp.LRT.4cp.tau[,1]<1 | 
                                                       nocp.LRT.4cp.tau[,1] <2 | nocp.LRT.4cp.tau[,4]==20) &
                                                    #but is 3 distinct change-points
                                                    nocp.LRT.3cp.tau[,3]-nocp.LRT.3cp.tau[,2]>1 &
                                                    nocp.LRT.3cp.tau[,2]-nocp.LRT.3cp.tau[,1]>1 & 
                                                    nocp.LRT.3cp.tau[,1] >=2 & nocp.LRT.3cp.tau[,3]<20))]
    
    
    nocp.LRT.test4cp<-nocp.LRT.4cp.value[which(nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.1cp & 
                                                 nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.2cp & 
                                                 nocp.BICres$nocp.BIC.4cp<nocp.BICres$nocp.BIC.3cp &
                                                 nocp.LRT.4cp.tau[,4]-nocp.LRT.4cp.tau[,3]>1 & 
                                                 nocp.LRT.4cp.tau[,3]-nocp.LRT.4cp.tau[,2]>1 & 
                                                 nocp.LRT.4cp.tau[,2]-nocp.LRT.4cp.tau[,1]>1 & 
                                                 nocp.LRT.4cp.tau[,1] >=2 & nocp.LRT.4cp.tau[,4]<20)]
    
    ## Type I error (run for nocp results)
    nocp.alpha<-sum(nocp.LRT.test1cp<probs1cp.alpha, nocp.LRT.test2cp<probs2cp.alpha, 
                    nocp.LRT.test3cp<probs3cp.alpha, nocp.LRT.test4cp<probs4cp.alpha)/nsim
    nocp.spending<-sum(nocp.LRT.test1cp<probs1cp.spending, nocp.LRT.test2cp<probs2cp.spending, 
                       nocp.LRT.test3cp<probs3cp.spending, nocp.LRT.test4cp<probs4cp.spending)/nsim
    
    typeI <- c(nocp.alpha, nocp.spending)
    names(typeI) <- c("Type I error", "Type I error (Sp)")
    
    return(typeI)
  }
  
  if(cp == 1) {
    #### 1cp ####
    cp1.BICres<-data.frame(cp1.BIC.1cp = resdta$BIC1cp,
                           cp1.BIC.2cp = resdta$BIC2cp, 
                           cp1.BIC.3cp = resdta$BIC3cp, 
                           cp1.BIC.4cp = resdta$BIC4cp)
    
    cp1.LRT.1cp.value <- resdta$LRT1cp
    cp1.LRT.2cp.value <- resdta$LRT2cp
    cp1.LRT.3cp.value <- resdta$LRT3cp
    cp1.LRT.4cp.value <- resdta$LRT4cp
    cp1.LRT.1cp.tau <- resdta$tau_1cp
    cp1.LRT.2cp.tau <- cbind(resdta$tau1_2cp, resdta$tau2_2cp)
    cp1.LRT.3cp.tau <- cbind(resdta$tau1_3cp, resdta$tau2_3cp, resdta$tau3_3cp)
    cp1.LRT.4cp.tau <- cbind(resdta$tau1_4cp, resdta$tau2_4cp, resdta$tau3_4cp, resdta$tau4_4cp)
    
    
    ## separate into testing datasets
    cp1.LRT.test1cp<-cp1.LRT.1cp.value[which((cp1.BICres$cp1.BIC.1cp<cp1.BICres$cp1.BIC.2cp & 
                                                cp1.BICres$cp1.BIC.1cp<cp1.BICres$cp1.BIC.3cp & 
                                                cp1.BICres$cp1.BIC.1cp<cp1.BICres$cp1.BIC.4cp) |
                                               #2cp that should be 1 cp
                                               (cp1.BICres$cp1.BIC.2cp<cp1.BICres$cp1.BIC.1cp & 
                                                  cp1.BICres$cp1.BIC.2cp<cp1.BICres$cp1.BIC.3cp & 
                                                  cp1.BICres$cp1.BIC.2cp<cp1.BICres$cp1.BIC.4cp &
                                                  (cp1.LRT.2cp.tau[,2]-cp1.LRT.2cp.tau[,1]<1 |
                                                     cp1.LRT.2cp.tau[,1]<2 | cp1.LRT.2cp.tau[,2]==20)) |
                                               #3cp that should be 1cp
                                               (cp1.BICres$cp1.BIC.3cp<cp1.BICres$cp1.BIC.1cp & 
                                                  cp1.BICres$cp1.BIC.3cp<cp1.BICres$cp1.BIC.2cp & 
                                                  cp1.BICres$cp1.BIC.3cp<cp1.BICres$cp1.BIC.4cp &
                                                  (cp1.LRT.3cp.tau[,2]-cp1.LRT.3cp.tau[,1]<1 |
                                                     cp1.LRT.3cp.tau[,3]-cp1.LRT.3cp.tau[,2]<1 | 
                                                     cp1.LRT.3cp.tau[,1]<2 | cp1.LRT.3cp.tau[,3]==20) &
                                                  (cp1.LRT.2cp.tau[,2]-cp1.LRT.2cp.tau[,1]<1 |
                                                     cp1.LRT.2cp.tau[,1]<2 | cp1.LRT.2cp.tau[,2]==20)) |
                                               #4cp that should be 1 cp
                                               (cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.1cp & 
                                                  cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.2cp & 
                                                  cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.3cp &
                                                  (cp1.LRT.4cp.tau[,4]-cp1.LRT.4cp.tau[,3]<1 | 
                                                     cp1.LRT.4cp.tau[,3]-cp1.LRT.4cp.tau[,2]<1 | 
                                                     cp1.LRT.4cp.tau[,2]-cp1.LRT.4cp.tau[,1]<1 | 
                                                     cp1.LRT.4cp.tau[,1]<2 | cp1.LRT.4cp.tau[,4]==20) &
                                                  (cp1.LRT.3cp.tau[,3]-cp1.LRT.3cp.tau[,2]<1 |
                                                     cp1.LRT.3cp.tau[,2]-cp1.LRT.3cp.tau[,1]<1 | 
                                                     cp1.LRT.3cp.tau[,1]<2 | cp1.LRT.3cp.tau[,3]==20) &
                                                  (cp1.LRT.2cp.tau[,2]-cp1.LRT.2cp.tau[,1]<1 |
                                                     cp1.LRT.2cp.tau[,1]<2 | cp1.LRT.2cp.tau[,2]==20)))]
    
    
    cp1.LRT.test2cp<-cp1.LRT.2cp.value[which((cp1.BICres$cp1.BIC.2cp<cp1.BICres$cp1.BIC.1cp & 
                                                cp1.BICres$cp1.BIC.2cp<cp1.BICres$cp1.BIC.3cp & 
                                                cp1.BICres$cp1.BIC.2cp<cp1.BICres$cp1.BIC.4cp &
                                                cp1.LRT.2cp.tau[,2]-cp1.LRT.2cp.tau[,1]>1 &
                                                cp1.LRT.2cp.tau[,1]>=2 & cp1.LRT.2cp.tau[,2]<20)|
                                               #3cp that should be 2cp 
                                               (cp1.BICres$cp1.BIC.3cp<cp1.BICres$cp1.BIC.1cp & 
                                                  cp1.BICres$cp1.BIC.3cp<cp1.BICres$cp1.BIC.2cp & 
                                                  cp1.BICres$cp1.BIC.3cp<cp1.BICres$cp1.BIC.4cp &
                                                  (cp1.LRT.3cp.tau[,2]-cp1.LRT.3cp.tau[,1]<1 |
                                                     cp1.LRT.3cp.tau[,3]-cp1.LRT.3cp.tau[,2]<1 | 
                                                     cp1.LRT.3cp.tau[,1]<2 | cp1.LRT.3cp.tau[,3]==20) &
                                                  cp1.LRT.2cp.tau[,2]-cp1.LRT.2cp.tau[,1]>1 & 
                                                  cp1.LRT.2cp.tau[,1]>=2 & cp1.LRT.2cp.tau[,2]<20)|
                                               #4cp that should be 2cp
                                               (cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.1cp & 
                                                  cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.2cp & 
                                                  cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.3cp &
                                                  (cp1.LRT.4cp.tau[,4]-cp1.LRT.4cp.tau[,3]<1 | 
                                                     cp1.LRT.4cp.tau[,3]-cp1.LRT.4cp.tau[,2]<1 | 
                                                     cp1.LRT.4cp.tau[,2]-cp1.LRT.4cp.tau[,1]<1 | 
                                                     cp1.LRT.4cp.tau[,1]<2 | cp1.LRT.4cp.tau[,4]==20) &
                                                  (cp1.LRT.3cp.tau[,3]-cp1.LRT.3cp.tau[,2]<1 |
                                                     cp1.LRT.3cp.tau[,2]-cp1.LRT.3cp.tau[,1]<1 | 
                                                     cp1.LRT.3cp.tau[,1]<2 | cp1.LRT.3cp.tau[,3]==20) &
                                                  cp1.LRT.2cp.tau[,2]-cp1.LRT.2cp.tau[,1]>1 & 
                                                  cp1.LRT.2cp.tau[,1]>=2 & cp1.LRT.2cp.tau[,2]<20))]
    
    
    cp1.LRT.test3cp<-cp1.LRT.3cp.value[which((cp1.BICres$cp1.BIC.3cp<cp1.BICres$cp1.BIC.1cp & 
                                                cp1.BICres$cp1.BIC.3cp<cp1.BICres$cp1.BIC.2cp & 
                                                cp1.BICres$cp1.BIC.3cp<cp1.BICres$cp1.BIC.4cp &
                                                cp1.LRT.3cp.tau[,2]-cp1.LRT.3cp.tau[,1]>1 &
                                                cp1.LRT.3cp.tau[,3]-cp1.LRT.3cp.tau[,2]>1 & 
                                                cp1.LRT.3cp.tau[,1] >=2 & cp1.LRT.3cp.tau[,3]<20)|
                                               #4cp should be 3cp
                                               (cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.1cp & 
                                                  cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.2cp & 
                                                  cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.3cp &
                                                  #not actually 4 change-points
                                                  (cp1.LRT.4cp.tau[,4]-cp1.LRT.4cp.tau[,3]<1 | 
                                                     cp1.LRT.4cp.tau[,3]-cp1.LRT.4cp.tau[,2]<1 | 
                                                     cp1.LRT.4cp.tau[,2]-cp1.LRT.4cp.tau[,1]<1 | 
                                                     cp1.LRT.4cp.tau[,1] <2 | cp1.LRT.4cp.tau[,4]==20) &
                                                  #but is 3 distinct change-points
                                                  cp1.LRT.3cp.tau[,3]-cp1.LRT.3cp.tau[,2]>1 &
                                                  cp1.LRT.3cp.tau[,2]-cp1.LRT.3cp.tau[,1]>1 & 
                                                  cp1.LRT.3cp.tau[,1] >=2 & cp1.LRT.3cp.tau[,3]<20))]
    
    
    cp1.LRT.test4cp<-cp1.LRT.4cp.value[which(cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.1cp & 
                                               cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.2cp & 
                                               cp1.BICres$cp1.BIC.4cp<cp1.BICres$cp1.BIC.3cp &
                                               cp1.LRT.4cp.tau[,4]-cp1.LRT.4cp.tau[,3]>1 & 
                                               cp1.LRT.4cp.tau[,3]-cp1.LRT.4cp.tau[,2]>1 & 
                                               cp1.LRT.4cp.tau[,2]-cp1.LRT.4cp.tau[,1]>1 & 
                                               cp1.LRT.4cp.tau[,1] >=2 & cp1.LRT.4cp.tau[,4]<20)]
    
    ##alpha
    #power
    cp1.power.alpha<-sum(cp1.LRT.test1cp<probs1cp.alpha, cp1.LRT.test2cp<probs2cp.alpha, cp1.LRT.test3cp<probs3cp.alpha, 
                         cp1.LRT.test4cp<probs4cp.alpha)/nsim 
    #sensitivity
    cp1.sensitivity.alpha<-sum(cp1.LRT.test1cp<probs1cp.alpha)/nsim 
    ##alpha spending
    #power
    cp1.power.spending<-sum(cp1.LRT.test1cp<probs1cp.spending, cp1.LRT.test2cp<probs2cp.spending, cp1.LRT.test3cp<probs3cp.spending, 
                            cp1.LRT.test4cp<probs4cp.spending)/nsim 
    #sensitivity
    cp1.sensitivity.spending<-sum(cp1.LRT.test1cp<probs1cp.spending)/nsim 
    
    powerstats <- c(cp1.power.alpha, cp1.sensitivity.alpha, 
                    cp1.power.spending, cp1.sensitivity.spending)
    names(powerstats) <- c("Power", "Sensitivity", 
                           "Power (Sp)", "Sensitivity (Sp)")
    return(powerstats)
  }
  
  if(cp == 2) {
    #### 3cp ####
    cp2.BICres<-data.frame(cp2.BIC.1cp = resdta$BIC1cp,
                           cp2.BIC.2cp = resdta$BIC2cp, 
                           cp2.BIC.3cp = resdta$BIC3cp, 
                           cp2.BIC.4cp = resdta$BIC4cp)
    
    cp2.LRT.1cp.value <- resdta$LRT1cp
    cp2.LRT.2cp.value <- resdta$LRT2cp
    cp2.LRT.3cp.value <- resdta$LRT3cp
    cp2.LRT.4cp.value <- resdta$LRT4cp
    cp2.LRT.1cp.tau <- resdta$tau_1cp
    cp2.LRT.2cp.tau <- cbind(resdta$tau1_2cp, resdta$tau2_2cp)
    cp2.LRT.3cp.tau <- cbind(resdta$tau1_3cp, resdta$tau2_3cp, resdta$tau3_3cp)
    cp2.LRT.4cp.tau <- cbind(resdta$tau1_4cp, resdta$tau2_4cp, resdta$tau3_4cp, resdta$tau4_4cp)
    
    ## separate into testing datasets
    cp2.LRT.test1cp<-cp2.LRT.1cp.value[which((cp2.BICres$cp2.BIC.1cp<cp2.BICres$cp2.BIC.2cp & 
                                                cp2.BICres$cp2.BIC.1cp<cp2.BICres$cp2.BIC.3cp & 
                                                cp2.BICres$cp2.BIC.1cp<cp2.BICres$cp2.BIC.4cp) |
                                               #2cp that should be 1 cp
                                               (cp2.BICres$cp2.BIC.2cp<cp2.BICres$cp2.BIC.1cp & 
                                                  cp2.BICres$cp2.BIC.2cp<cp2.BICres$cp2.BIC.3cp & 
                                                  cp2.BICres$cp2.BIC.2cp<cp2.BICres$cp2.BIC.4cp &
                                                  (cp2.LRT.2cp.tau[,2]-cp2.LRT.2cp.tau[,1]<1 |
                                                     cp2.LRT.2cp.tau[,1]<2 | cp2.LRT.2cp.tau[,2]==20)) |
                                               #3cp that should be 1cp
                                               (cp2.BICres$cp2.BIC.3cp<cp2.BICres$cp2.BIC.1cp & 
                                                  cp2.BICres$cp2.BIC.3cp<cp2.BICres$cp2.BIC.2cp & 
                                                  cp2.BICres$cp2.BIC.3cp<cp2.BICres$cp2.BIC.4cp &
                                                  (cp2.LRT.3cp.tau[,2]-cp2.LRT.3cp.tau[,1]<1 |
                                                     cp2.LRT.3cp.tau[,3]-cp2.LRT.3cp.tau[,2]<1 | 
                                                     cp2.LRT.3cp.tau[,1]<2 | cp2.LRT.3cp.tau[,3]==20) &
                                                  (cp2.LRT.2cp.tau[,2]-cp2.LRT.2cp.tau[,1]<1 |
                                                     cp2.LRT.2cp.tau[,1]<2 | cp2.LRT.2cp.tau[,2]==20)) |
                                               #4cp that should be 1 cp
                                               (cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.1cp & 
                                                  cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.2cp & 
                                                  cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.3cp &
                                                  (cp2.LRT.4cp.tau[,4]-cp2.LRT.4cp.tau[,3]<1 | 
                                                     cp2.LRT.4cp.tau[,3]-cp2.LRT.4cp.tau[,2]<1 | 
                                                     cp2.LRT.4cp.tau[,2]-cp2.LRT.4cp.tau[,1]<1 | 
                                                     cp2.LRT.4cp.tau[,1]<2 | cp2.LRT.4cp.tau[,4]==20) &
                                                  (cp2.LRT.3cp.tau[,3]-cp2.LRT.3cp.tau[,2]<1 |
                                                     cp2.LRT.3cp.tau[,2]-cp2.LRT.3cp.tau[,1]<1 | 
                                                     cp2.LRT.3cp.tau[,1]<2 | cp2.LRT.3cp.tau[,3]==20) &
                                                  (cp2.LRT.2cp.tau[,2]-cp2.LRT.2cp.tau[,1]<1 |
                                                     cp2.LRT.2cp.tau[,1]<2 | cp2.LRT.2cp.tau[,2]==20)))]
    
    
    cp2.LRT.test2cp<-cp2.LRT.2cp.value[which((cp2.BICres$cp2.BIC.2cp<cp2.BICres$cp2.BIC.1cp & 
                                                cp2.BICres$cp2.BIC.2cp<cp2.BICres$cp2.BIC.3cp & 
                                                cp2.BICres$cp2.BIC.2cp<cp2.BICres$cp2.BIC.4cp &
                                                cp2.LRT.2cp.tau[,2]-cp2.LRT.2cp.tau[,1]>1 & 
                                                cp2.LRT.2cp.tau[,1]>=2 & cp2.LRT.2cp.tau[,2]<20)|
                                               #3cp that should be 2cp 
                                               (cp2.BICres$cp2.BIC.3cp<cp2.BICres$cp2.BIC.1cp & 
                                                  cp2.BICres$cp2.BIC.3cp<cp2.BICres$cp2.BIC.2cp & 
                                                  cp2.BICres$cp2.BIC.3cp<cp2.BICres$cp2.BIC.4cp &
                                                  (cp2.LRT.3cp.tau[,2]-cp2.LRT.3cp.tau[,1]<1 |
                                                     cp2.LRT.3cp.tau[,3]-cp2.LRT.3cp.tau[,2]<1 | 
                                                     cp2.LRT.3cp.tau[,1]<2 | cp2.LRT.3cp.tau[,3]==20) &
                                                  cp2.LRT.2cp.tau[,2]-cp2.LRT.2cp.tau[,1]>1 & 
                                                  cp2.LRT.2cp.tau[,1]>=2 & cp2.LRT.2cp.tau[,2]<20)|
                                               #4cp that should be 2cp
                                               (cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.1cp & 
                                                  cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.2cp & 
                                                  cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.3cp &
                                                  (cp2.LRT.4cp.tau[,4]-cp2.LRT.4cp.tau[,3]<1 | 
                                                     cp2.LRT.4cp.tau[,3]-cp2.LRT.4cp.tau[,2]<1 | 
                                                     cp2.LRT.4cp.tau[,2]-cp2.LRT.4cp.tau[,1]<1 | 
                                                     cp2.LRT.4cp.tau[,1]<2 | cp2.LRT.4cp.tau[,4]==20) &
                                                  (cp2.LRT.3cp.tau[,3]-cp2.LRT.3cp.tau[,2]<1 |
                                                     cp2.LRT.3cp.tau[,2]-cp2.LRT.3cp.tau[,1]<1 | 
                                                     cp2.LRT.3cp.tau[,1]<2 | cp2.LRT.3cp.tau[,3]==20) &
                                                  cp2.LRT.2cp.tau[,2]-cp2.LRT.2cp.tau[,1]>1 & 
                                                  cp2.LRT.2cp.tau[,1]>=2 & cp2.LRT.2cp.tau[,2]<20))]
    
    
    cp2.LRT.test3cp<-cp2.LRT.3cp.value[which((cp2.BICres$cp2.BIC.3cp<cp2.BICres$cp2.BIC.1cp & 
                                                cp2.BICres$cp2.BIC.3cp<cp2.BICres$cp2.BIC.2cp & 
                                                cp2.BICres$cp2.BIC.3cp<cp2.BICres$cp2.BIC.4cp &
                                                cp2.LRT.3cp.tau[,2]-cp2.LRT.3cp.tau[,1]>1 &
                                                cp2.LRT.3cp.tau[,3]-cp2.LRT.3cp.tau[,2]>1 & 
                                                cp2.LRT.3cp.tau[,1] >=2 & cp2.LRT.3cp.tau[,3]<20)|
                                               #4cp should be 3cp
                                               (cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.1cp & 
                                                  cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.2cp & 
                                                  cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.3cp &
                                                  #not actually 4 change-points
                                                  (cp2.LRT.4cp.tau[,4]-cp2.LRT.4cp.tau[,3]<1 | 
                                                     cp2.LRT.4cp.tau[,3]-cp2.LRT.4cp.tau[,2]<1 | 
                                                     cp2.LRT.4cp.tau[,2]-cp2.LRT.4cp.tau[,1]<1 | 
                                                     cp2.LRT.4cp.tau[,1] <2 | cp2.LRT.4cp.tau[,4]==20) &
                                                  #but is 3 distinct change-points
                                                  cp2.LRT.3cp.tau[,3]-cp2.LRT.3cp.tau[,2]>1 &
                                                  cp2.LRT.3cp.tau[,2]-cp2.LRT.3cp.tau[,1]>1 & 
                                                  cp2.LRT.3cp.tau[,1] >=2 & cp2.LRT.3cp.tau[,3]<20))]
    
    
    cp2.LRT.test4cp<-cp2.LRT.4cp.value[which(cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.1cp & 
                                               cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.2cp & 
                                               cp2.BICres$cp2.BIC.4cp<cp2.BICres$cp2.BIC.3cp &
                                               cp2.LRT.4cp.tau[,4]-cp2.LRT.4cp.tau[,3]>1 & 
                                               cp2.LRT.4cp.tau[,3]-cp2.LRT.4cp.tau[,2]>1 & 
                                               cp2.LRT.4cp.tau[,2]-cp2.LRT.4cp.tau[,1]>1 & 
                                               cp2.LRT.4cp.tau[,1] >=2 & cp2.LRT.4cp.tau[,4]<20)]
    
    ## run for 2cp results
    # alpha
    #power
    cp2.power.alpha<-sum(cp2.LRT.test1cp<probs1cp.alpha, cp2.LRT.test2cp<probs2cp.alpha, 
                         cp2.LRT.test3cp<probs3cp.alpha, cp2.LRT.test4cp<probs4cp.alpha)/nsim 
    
    #sensitivity
    cp2.sensitivity.alpha<-sum(cp2.LRT.test2cp<probs2cp.alpha)/nsim 
    
    #alpha spending
    #power
    cp2.power.spending<-sum(cp2.LRT.test1cp<probs1cp.spending, cp2.LRT.test2cp<probs2cp.spending, 
                            cp2.LRT.test3cp<probs3cp.spending, cp2.LRT.test4cp<probs4cp.spending)/nsim 
    
    #sensitivity
    cp2.sensitivity.spending<-sum(cp2.LRT.test2cp<probs2cp.spending)/nsim 
    
    powerstats <- c(cp2.power.alpha, cp2.sensitivity.alpha, 
                    cp2.power.spending, cp2.sensitivity.spending)
    names(powerstats) <- c("Power", "Sensitivity", 
                           "Power (Sp)", "Sensitivity (Sp)")
    return(powerstats)
  }
  
  if(cp == 3) {
    #### 4cp ####
    cp3.BICres<-data.frame(cp3.BIC.1cp = resdta$BIC1cp,
                           cp3.BIC.2cp = resdta$BIC2cp, 
                           cp3.BIC.3cp = resdta$BIC3cp,
                           cp3.BIC.4cp = resdta$BIC4cp)
    
    cp3.LRT.1cp.value <- resdta$LRT1cp
    cp3.LRT.2cp.value <- resdta$LRT2cp
    cp3.LRT.3cp.value <- resdta$LRT3cp
    cp3.LRT.4cp.value <- resdta$LRT4cp
    cp3.LRT.1cp.tau <- resdta$tau_1cp
    cp3.LRT.2cp.tau <- cbind(resdta$tau1_2cp, resdta$tau2_2cp)
    cp3.LRT.3cp.tau <- cbind(resdta$tau1_3cp, resdta$tau2_3cp, resdta$tau3_3cp)
    cp3.LRT.4cp.tau <- cbind(resdta$tau1_4cp, resdta$tau2_4cp, resdta$tau3_4cp, resdta$tau4_4cp)
    
    ## separate into testing datasets
    cp3.LRT.test1cp<-cp3.LRT.1cp.value[which((cp3.BICres$cp3.BIC.1cp<cp3.BICres$cp3.BIC.2cp & 
                                                cp3.BICres$cp3.BIC.1cp<cp3.BICres$cp3.BIC.3cp & 
                                                cp3.BICres$cp3.BIC.1cp<cp3.BICres$cp3.BIC.4cp) |
                                               #2cp that should be 1 cp
                                               (cp3.BICres$cp3.BIC.2cp<cp3.BICres$cp3.BIC.1cp & 
                                                  cp3.BICres$cp3.BIC.2cp<cp3.BICres$cp3.BIC.3cp & 
                                                  cp3.BICres$cp3.BIC.2cp<cp3.BICres$cp3.BIC.4cp &
                                                  (cp3.LRT.2cp.tau[,2]-cp3.LRT.2cp.tau[,1]<1 |
                                                     cp3.LRT.2cp.tau[,1]<2 | cp3.LRT.2cp.tau[,2]==20)) |
                                               #3cp that should be 1cp
                                               (cp3.BICres$cp3.BIC.3cp<cp3.BICres$cp3.BIC.1cp & 
                                                  cp3.BICres$cp3.BIC.3cp<cp3.BICres$cp3.BIC.2cp & 
                                                  cp3.BICres$cp3.BIC.3cp<cp3.BICres$cp3.BIC.4cp &
                                                  (cp3.LRT.3cp.tau[,2]-cp3.LRT.3cp.tau[,1]<1 |
                                                     cp3.LRT.3cp.tau[,3]-cp3.LRT.3cp.tau[,2]<1 | 
                                                     cp3.LRT.3cp.tau[,1]<2 | cp3.LRT.3cp.tau[,3]==20) &
                                                  (cp3.LRT.2cp.tau[,2]-cp3.LRT.2cp.tau[,1]<1 |
                                                     cp3.LRT.2cp.tau[,1]<2 | cp3.LRT.2cp.tau[,2]==20)) |
                                               #4cp that should be 1 cp
                                               (cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.1cp & 
                                                  cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.2cp & 
                                                  cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.3cp &
                                                  (cp3.LRT.4cp.tau[,4]-cp3.LRT.4cp.tau[,3]<1 | 
                                                     cp3.LRT.4cp.tau[,3]-cp3.LRT.4cp.tau[,2]<1 | 
                                                     cp3.LRT.4cp.tau[,2]-cp3.LRT.4cp.tau[,1]<1 | 
                                                     cp3.LRT.4cp.tau[,1]<2 | cp3.LRT.4cp.tau[,4]==20) &
                                                  (cp3.LRT.3cp.tau[,3]-cp3.LRT.3cp.tau[,2]<1 |
                                                     cp3.LRT.3cp.tau[,2]-cp3.LRT.3cp.tau[,1]<1 | 
                                                     cp3.LRT.3cp.tau[,1]<2 | cp3.LRT.3cp.tau[,3]==20) &
                                                  (cp3.LRT.2cp.tau[,2]-cp3.LRT.2cp.tau[,1]<1 |
                                                     cp3.LRT.2cp.tau[,1]<2 | cp3.LRT.2cp.tau[,2]==20)))]
    
    
    cp3.LRT.test2cp<-cp3.LRT.2cp.value[which((cp3.BICres$cp3.BIC.2cp<cp3.BICres$cp3.BIC.1cp & 
                                                cp3.BICres$cp3.BIC.2cp<cp3.BICres$cp3.BIC.3cp & 
                                                cp3.BICres$cp3.BIC.2cp<cp3.BICres$cp3.BIC.4cp &
                                                cp3.LRT.2cp.tau[,2]-cp3.LRT.2cp.tau[,1]>1 & 
                                                cp3.LRT.2cp.tau[,1]>=2 & cp3.LRT.2cp.tau[,2]<20)|
                                               #3cp that should be 2cp 
                                               (cp3.BICres$cp3.BIC.3cp<cp3.BICres$cp3.BIC.1cp & 
                                                  cp3.BICres$cp3.BIC.3cp<cp3.BICres$cp3.BIC.2cp & 
                                                  cp3.BICres$cp3.BIC.3cp<cp3.BICres$cp3.BIC.4cp &
                                                  (cp3.LRT.3cp.tau[,2]-cp3.LRT.3cp.tau[,1]<1 |
                                                     cp3.LRT.3cp.tau[,3]-cp3.LRT.3cp.tau[,2]<1 | 
                                                     cp3.LRT.3cp.tau[,1]<2 | cp3.LRT.3cp.tau[,3]==20) &
                                                  cp3.LRT.2cp.tau[,2]-cp3.LRT.2cp.tau[,1]>1 & 
                                                  cp3.LRT.2cp.tau[,1]>=2 & cp3.LRT.2cp.tau[,2]<20)|
                                               #4cp that should be 2cp
                                               (cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.1cp & 
                                                  cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.2cp & 
                                                  cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.3cp &
                                                  (cp3.LRT.4cp.tau[,4]-cp3.LRT.4cp.tau[,3]<1 | 
                                                     cp3.LRT.4cp.tau[,3]-cp3.LRT.4cp.tau[,2]<1 | 
                                                     cp3.LRT.4cp.tau[,2]-cp3.LRT.4cp.tau[,1]<1 | 
                                                     cp3.LRT.4cp.tau[,1]<2 | cp3.LRT.4cp.tau[,4]==20) &
                                                  (cp3.LRT.3cp.tau[,3]-cp3.LRT.3cp.tau[,2]<1 |
                                                     cp3.LRT.3cp.tau[,2]-cp3.LRT.3cp.tau[,1]<1 | 
                                                     cp3.LRT.3cp.tau[,1]<2 | cp3.LRT.3cp.tau[,3]==20) &
                                                  cp3.LRT.2cp.tau[,2]-cp3.LRT.2cp.tau[,1]>1 & 
                                                  cp3.LRT.2cp.tau[,1]>=2 & cp3.LRT.2cp.tau[,2]<20))]
    
    
    cp3.LRT.test3cp<-cp3.LRT.3cp.value[which((cp3.BICres$cp3.BIC.3cp<cp3.BICres$cp3.BIC.1cp & 
                                                cp3.BICres$cp3.BIC.3cp<cp3.BICres$cp3.BIC.2cp & 
                                                cp3.BICres$cp3.BIC.3cp<cp3.BICres$cp3.BIC.4cp &
                                                cp3.LRT.3cp.tau[,2]-cp3.LRT.3cp.tau[,1]>1 &
                                                cp3.LRT.3cp.tau[,3]-cp3.LRT.3cp.tau[,2]>1 & 
                                                cp3.LRT.3cp.tau[,1] >=2 & cp3.LRT.3cp.tau[,3]<20)|
                                               #4cp should be 3cp
                                               (cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.1cp & 
                                                  cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.2cp & 
                                                  cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.3cp &
                                                  #not actually 4 change-points
                                                  (cp3.LRT.4cp.tau[,4]-cp3.LRT.4cp.tau[,3]<1 | 
                                                     cp3.LRT.4cp.tau[,3]-cp3.LRT.4cp.tau[,2]<1 | 
                                                     cp3.LRT.4cp.tau[,2]-cp3.LRT.4cp.tau[,1]<1 | 
                                                     cp3.LRT.4cp.tau[,1] <2 | cp3.LRT.4cp.tau[,4]==20) &
                                                  #but is 3 distinct change-points
                                                  cp3.LRT.3cp.tau[,3]-cp3.LRT.3cp.tau[,2]>1 &
                                                  cp3.LRT.3cp.tau[,2]-cp3.LRT.3cp.tau[,1]>1 & 
                                                  cp3.LRT.3cp.tau[,1] >=2 & cp3.LRT.3cp.tau[,3]<20))]
    
    
    cp3.LRT.test4cp<-cp3.LRT.4cp.value[which(cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.1cp & 
                                               cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.2cp & 
                                               cp3.BICres$cp3.BIC.4cp<cp3.BICres$cp3.BIC.3cp &
                                               cp3.LRT.4cp.tau[,4]-cp3.LRT.4cp.tau[,3]>1 & 
                                               cp3.LRT.4cp.tau[,3]-cp3.LRT.4cp.tau[,2]>1 & 
                                               cp3.LRT.4cp.tau[,2]-cp3.LRT.4cp.tau[,1]>1 & 
                                               cp3.LRT.4cp.tau[,1] >=2 & cp3.LRT.4cp.tau[,4]<20)]
    
    ## run for 3cp results
    #alpha
    #power
    cp3.power.alpha<-sum(cp3.LRT.test1cp<probs1cp.alpha, cp3.LRT.test2cp<probs2cp.alpha, 
                         cp3.LRT.test3cp<probs3cp.alpha, cp3.LRT.test4cp<probs4cp.alpha)/nsim 
    
    #sensitivity
    cp3.sensitivity.alpha<-sum(cp3.LRT.test3cp<probs3cp.alpha)/nsim 
    
    #alpha spending
    cp3.power.spending<-sum(cp3.LRT.test1cp<probs1cp.spending, cp3.LRT.test2cp<probs2cp.spending, 
                            cp3.LRT.test3cp<probs3cp.spending, cp3.LRT.test4cp<probs4cp.spending)/nsim 
    
    #sensitivity
    cp3.sensitivity.spending<-sum(cp3.LRT.test3cp<probs3cp.spending)/nsim 
    
    powerstats <- c(cp3.power.alpha, cp3.sensitivity.alpha, 
                    cp3.power.spending, cp3.sensitivity.spending)
    
    names(powerstats) <- c("Power", "Sensitivity", 
                           "Power (Sp)", "Sensitivity (Sp)")
    return(powerstats)
  }
}


save(weib.1cp.nll, weib.1cp.LRT, weib.1cp.BIC,
     weib.2cp.nll, weib.2cp.LRT, weib.2cp.BIC,
     weib.3cp.nll, weib.3cp.LRT, weib.3cp.BIC,
     weib.4cp.nll, weib.4cp.LRT, weib.4cp.BIC,
     weib_invcdf, 
     exp.1cp.nll, exp.1cp.LRT, exp.1cp.BIC,
     exp.2cp.nll, exp.2cp.LRT, exp.2cp.BIC,
     exp.3cp.nll, exp.3cp.LRT, exp.3cp.BIC,
     exp.4cp.nll, exp.4cp.LRT, exp.4cp.BIC,
     exp_invcdf,
     is.error, test_cp, results,
     file = "C:/Users/Camille/Documents/D word/4 Writing the thing/simulation 2019 11 20/dwordfns20191206.RData")
