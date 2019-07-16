## Functions for applying test for multiple
## change-points within the Weibull change-point
## hazard model

##### overall functions ##### 

## par[1] is gamma, par[2] is theta0
nll <- function(par, dta) {
  ll <- log(par[2]) * sum(dta$censor) + 
          (par[1] - 1) * sum(dta$censor * log(dta$time)) - 
          (par[2]/par[1]) * sum(dta$time^(par[1]))
  return(-ll)
}

##### functions for one change-point #####

# MLE log-likelihood
nll.1cp <- function(tau, gamma = 2) {
  g <- gamma
  G <- sum(dta$censor)
  T <- (1/g) * sum(dta$time^g)
  G1 <- sum(dta$time < tau)
  G2 <- G - G1
  T2 <- (1/g) * sum((dta$time^g - tau^g) * 
                      (dta$time >= tau))
  T1 <- T - T2
  # using MLEs for thetas
  ll <- G1 * log(G1/T1) + G2 * log(G2/T2) - G1 - 
        G2 + (g - 1) * sum(dta$censor * log(dta$time))
  return(-ll)
}

# LRT
LRT.1cp <- function(tau, gamma = 2) {
  g <- gamma
  G <- sum(dta$censor)
  T <- (1/g) * sum(dta$time^g)
  G1 <- sum(dta$time < tau)
  G2 <- G - G1
  T2 <- (1/g) * sum((dta$time^g - tau^g) * 
                      (dta$time >= tau))
  T1 <- T - T2
  lrt <- G1 * log((G1/T1) * (T/G)) + 
          G2 * log((G2/T2) * (T/G))
  return(-lrt)  #need to use minimum so that we can compare to optim fn
}


# BIC
BIC.1cp <- function(tau, gamma = 2) {
  3 * log(n) + 2 * nll.1cp(tau, gamma)
}

# full nll
cp1.nll.full <- function(par, tau, dta) {
  theta1 <- par[1]
  theta2 <- par[2]
  g <- par[3]
  ll <- (g - 1) * sum(dta$censor * log(dta$time)) + 
        log(theta1) * sum((dta$time < tau)) + 
        log(theta2) * sum((tau <= dta$time)* dta$censor) - 
        (theta1/g) * sum((dta$time^g) * (dta$time < tau) + 
                           (tau^g) * (dta$time >= tau)) - 
    (theta2/g) * sum((dta$time^g - tau^g) * 
                       (dta$time >= tau))
  return(-ll)
}

##### functions for two change-points #####

# MLE log-likelihood
nll.2cp <- function(tau, gamma = 2) {
  g <- gamma
  G <- sum(dta$censor)
  T <- (1/g) * sum(dta$time^g)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum(dta$censor * (dta$time >= tau[2]))
  T1 <- (1/g) * sum((dta$time^g) * (dta$time < tau[1]) + 
                      (tau[1]^g) * (tau[1] <= dta$time))
  T2 <- (1/g) * sum((dta$time^g - tau[1]^g) * (tau[1] <= 
                                                 dta$time) * (dta$time < tau[2]) + (tau[2]^g - 
                                                                                      tau[1]^g) * (dta$time >= tau[2]))
  T3 <- (1/g) * sum((dta$time^g - tau[2]^g) * (dta$time >= 
                                                 tau[2]))
  ll <- G1 * log(G1/T1) + G2 * log(G2/T2) + G3 * 
    log(G3/T3) - G1 - G2 - G3 + (g - 1) * sum(dta$censor * 
                                                log(dta$time))
  return(-ll)
}

# LRT
LRT.2cp <- function(tau, gamma = 2) {
  g <- gamma
  d <- dta$censor
  G <- sum(dta$censor)
  T <- (1/g) * sum(dta$time^g)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum(dta$censor * (dta$time >= tau[2]))
  T1 <- (1/g) * sum((dta$time^g) * (dta$time < tau[1]) + 
                      (tau[1]^g) * (tau[1] <= dta$time))
  T2 <- (1/g) * sum((dta$time^g - tau[1]^g) * (tau[1] <= 
                                                 dta$time) * (dta$time < tau[2]) + (tau[2]^g - 
                                                                                      tau[1]^g) * (dta$time >= tau[2]))
  T3 <- (1/g) * sum((dta$time^g - tau[2]^g) * (dta$time >= 
                                                 tau[2]))
  lrt <- G1 * log((G1/T1) * (T/G)) + G2 * log((G2/T2) * 
                                                (T/G)) + G3 * log((G3/T3) * (T/G))
  return(-lrt)
}

# BIC
BIC.2cp <- function(tau, gamma = 2) {
  5 * log(n) + 2 * nll.2cp(tau, gamma)
}

# full nll
cp2.nll.full <- function(par, tau, dta) {
  theta1 <- par[1]
  theta2 <- par[2]
  theta3 <- par[3]
  g <- par[4]
  ll <- (g - 1) * sum(dta$censor * log(dta$time)) + 
    log(theta1) * sum((dta$time < tau[1])) + log(theta2) * 
    sum((tau[1] <= dta$time) * (dta$time < tau[2])) + 
    log(theta3) * sum((dta$time >= tau[2]) * dta$censor) - 
    (theta1/g) * sum((dta$time^g) * (dta$time < 
                                       tau[1]) + (tau[1]^g) * (dta$time >= tau[1])) - 
    (theta2/g) * sum((dta$time^g - tau[1]^g) * 
                       (dta$time >= tau[1]) * (dta$time < tau[2]) + 
                       (tau[2]^g - tau[1]^g) * (dta$time >= tau[2])) - 
    (theta3/g) * sum((dta$time^g - tau[2]^g) * 
                       (dta$time >= tau[2]))
  return(-ll)
}

##### functions for three change-points #####

# MLE log-likelihood
nll.3cp <- function(tau, gamma = 2) {
  g <- gamma
  G <- sum(dta$censor)
  T <- (1/g) * sum(dta$time^g)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum((tau[2] <= dta$time) * (dta$time < tau[3]))
  G4 <- sum((dta$time >= tau[3]) * dta$censor)
  T1 <- (1/g) * sum((dta$time^g) * (dta$time < tau[1]) + 
                      (tau[1]^g) * (tau[1] <= dta$time))
  T2 <- (1/g) * sum((dta$time^g - tau[1]^g) * (tau[1] <= 
                                                 dta$time) * (dta$time < tau[2]) + (tau[2]^g - 
                                                                                      tau[1]^g) * (dta$time >= tau[2]))
  T3 <- (1/g) * sum((dta$time^g - tau[2]^g) * (tau[2] <= 
                                                 dta$time) * (dta$time < tau[3]) + (tau[3]^g - 
                                                                                      tau[2]^g) * (dta$time >= tau[3]))
  T4 <- (1/g) * sum((dta$time^g - tau[3]^g) * (dta$time >= 
                                                 tau[3]))
  # using MLEs for thetas
  ll <- G1 * log(G1/T1) + G2 * log(G2/T2) + G3 * 
    log(G3/T3) + G4 * log(G4/T4) - G1 - G2 - G3 - 
    G4 + (g - 1) * sum(dta$censor * log(dta$time))
  return(-ll)
}

# LRT
LRT.3cp <- function(tau, gamma = 2) {
  g <- gamma
  G <- sum(dta$censor)
  T <- (1/g) * sum(dta$time^g)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum((tau[2] <= dta$time) * (dta$time < tau[3]))
  G4 <- sum((dta$time >= tau[3]) * dta$censor)
  T1 <- (1/g) * sum((dta$time^g) * (dta$time < tau[1]) + 
                      (tau[1]^g) * (tau[1] <= dta$time))
  T2 <- (1/g) * sum((dta$time^g - tau[1]^g) * (tau[1] <= 
                                                 dta$time) * (dta$time < tau[2]) + (tau[2]^g - 
                                                                                      tau[1]^g) * (dta$time >= tau[2]))
  T3 <- (1/g) * sum((dta$time^g - tau[2]^g) * (tau[2] <= 
                                                 dta$time) * (dta$time < tau[3]) + (tau[3]^g - 
                                                                                      tau[2]^g) * (dta$time >= tau[3]))
  T4 <- (1/g) * sum((dta$time^g - tau[3]^g) * (dta$time >= 
                                                 tau[3]))
  lrt <- G1 * log((G1/T1) * (T/G)) + G2 * log((G2/T2) * 
                                                (T/G)) + G3 * log((G3/T3) * (T/G)) + G4 * log((G4/T4) * 
                                                                                                (T/G))
  return(-lrt)
}

# BIC
BIC.3cp <- function(tau, gamma = 2) {
  7 * log(n) + 2 * nll.3cp(tau, gamma)
}

# full nll
cp3.nll.full <- function(par, tau, dta) {
  theta1 <- par[1]
  theta2 <- par[2]
  theta3 <- par[3]
  theta4 <- par[4]
  g <- par[5]
  ll <- (g - 1) * sum(dta$censor * log(dta$time)) + 
    log(theta1) * sum((dta$time < tau[1])) + log(theta2) * 
    sum((tau[1] <= dta$time) * (dta$time < tau[2])) + 
    log(theta3) * sum((tau[2] <= dta$time) * (dta$time < 
                                                tau[3])) + log(theta4) * sum((dta$time >= 
                                                                                tau[3]) * dta$censor) - (theta1/g) * sum((dta$time^g) * 
                                                                                                                           (dta$time < tau[1]) + (tau[1]^g) * (dta$time >= 
                                                                                                                                                                 tau[1])) - (theta2/g) * sum((dta$time^g - tau[1]^g) * 
                                                                                                                                                                                               (dta$time >= tau[1]) * (dta$time < tau[2]) + 
                                                                                                                                                                                               (tau[2]^g - tau[1]^g) * (dta$time >= tau[2])) - 
    (theta3/g) * sum((dta$time^g - tau[2]^g) * 
                       (dta$time >= tau[2]) * (dta$time < tau[3]) + 
                       (tau[3]^g - tau[2]^g) * (dta$time >= tau[3])) - 
    (theta4/g) * sum((dta$time^g - tau[3]^g) * 
                       (dta$time >= tau[3]))
  return(-ll)
}


##### functions for four change-points #####

# MLE log-likelihood
nll.4cp <- function(tau, gamma = 2) {
  g <- gamma
  G <- sum(dta$censor)
  T <- (1/g) * sum(dta$time^g)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum((tau[2] <= dta$time) * (dta$time < tau[3]))
  G4 <- sum((tau[3] <= dta$time) * (dta$time < tau[4]))
  G5 <- sum((dta$time >= tau[4]) * dta$censor)
  T1 <- (1/g) * sum((dta$time^g) * (dta$time < tau[1]) + 
                      (tau[1]^g) * (tau[1] <= dta$time))
  T2 <- (1/g) * sum((dta$time^g - tau[1]^g) * (tau[1] <= 
                                                 dta$time) * (dta$time < tau[2]) + (tau[2]^g - 
                                                                                      tau[1]^g) * (dta$time >= tau[2]))
  T3 <- (1/g) * sum((dta$time^g - tau[2]^g) * (tau[2] <= 
                                                 dta$time) * (dta$time < tau[3]) + (tau[3]^g - 
                                                                                      tau[2]^g) * (dta$time >= tau[3]))
  T4 <- (1/g) * sum((dta$time^g - tau[3]^g) * (tau[3] <= 
                                                 dta$time) * (dta$time < tau[4]) + (tau[4]^g - 
                                                                                      tau[3]^g) * (dta$time >= tau[4]))
  T5 <- (1/g) * sum((dta$time^g - tau[4]^g) * (dta$time >= 
                                                 tau[4]))
  ll <- G1 * log(G1/T1) + G2 * log(G2/T2) + G3 * 
    log(G3/T3) + G4 * log(G4/T4) + G5 * log(G5/T5) - 
    G1 - G2 - G3 - G4 - G5 + (g - 1) * sum(dta$censor * 
                                             log(dta$time))
  return(-ll)
}

# LRT
LRT.4cp <- function(tau, gamma = 2) {
  g <- gamma
  G <- sum(dta$censor)
  T <- (1/g) * sum(dta$time^g)
  G1 <- sum(dta$time < tau[1])
  G2 <- sum((tau[1] <= dta$time) * (dta$time < tau[2]))
  G3 <- sum((tau[2] <= dta$time) * (dta$time < tau[3]))
  G4 <- sum((tau[3] <= dta$time) * (dta$time < tau[4]))
  G5 <- sum((dta$time >= tau[4]) * dta$censor)
  T1 <- (1/g) * sum((dta$time^g) * (dta$time < tau[1]) + 
                      (tau[1]^g) * (tau[1] <= dta$time))
  T2 <- (1/g) * sum((dta$time^g - tau[1]^g) * (tau[1] <= 
                                                 dta$time) * (dta$time < tau[2]) + (tau[2]^g - 
                                                                                      tau[1]^g) * (dta$time >= tau[2]))
  T3 <- (1/g) * sum((dta$time^g - tau[2]^g) * (tau[2] <= 
                                                 dta$time) * (dta$time < tau[3]) + (tau[3]^g - 
                                                                                      tau[2]^g) * (dta$time >= tau[3]))
  T4 <- (1/g) * sum((dta$time^g - tau[3]^g) * (tau[3] <= 
                                                 dta$time) * (dta$time < tau[4]) + (tau[4]^g - 
                                                                                      tau[3]^g) * (dta$time >= tau[4]))
  T5 <- (1/g) * sum((dta$time^g - tau[4]^g) * (dta$time >= 
                                                 tau[4]))
  lrt <- G1 * log((G1/T1) * (T/G)) + G2 * log((G2/T2) * 
                                                (T/G)) + G3 * log((G3/T3) * (T/G)) + G4 * log((G4/T4) * 
                                                                                                (T/G)) + G5 * log((G5/T5) * (T/G))
  return(-lrt)
}

# BIC
BIC.4cp <- function(tau, gamma = 2) {
  9 * log(n) + 2 * nll.4cp(tau, gamma)
}

# full nll
cp4.nll.full <- function(par, tau, dta) {
  theta1 <- par[1]
  theta2 <- par[2]
  theta3 <- par[3]
  theta4 <- par[4]
  theta5 <- par[5]
  g <- par[6]
  ll <- (g - 1) * sum(dta$censor * log(dta$time)) + 
    log(theta1) * sum((dta$time < tau[1])) + log(theta2) * 
    sum((tau[1] <= dta$time) * (dta$time < tau[2])) + 
    log(theta3) * sum((tau[2] <= dta$time) * (dta$time < 
                                                tau[3])) + log(theta4) * sum((tau[3] <= 
                                                                                dta$time) * (dta$time < tau[4])) + log(theta5) * 
    sum((dta$time >= tau[4]) * dta$censor) - (theta1/g) * 
    sum((dta$time^g) * (dta$time < tau[1]) + (tau[1]^g) * 
          (dta$time >= tau[1])) - (theta2/g) * sum((dta$time^g - 
                                                      tau[1]^g) * (dta$time >= tau[1]) * (dta$time < 
                                                                                            tau[2]) + (tau[2]^g - tau[1]^g) * (dta$time >= 
                                                                                                                                 tau[2])) - (theta3/g) * sum((dta$time^g - tau[2]^g) * 
                                                                                                                                                               (dta$time >= tau[2]) * (dta$time < tau[3]) + 
                                                                                                                                                               (tau[3]^g - tau[2]^g) * (dta$time >= tau[3])) - 
    (theta4/g) * sum((dta$time^g - tau[3]^g) * 
                       (dta$time >= tau[3]) * (dta$time < tau[4]) + 
                       (tau[4]^g - tau[3]^g) * (dta$time >= tau[4])) - 
    (theta5/g) * sum((dta$time^g - tau[4]^g) * 
                       (dta$time >= tau[4]))
  return(-ll)
}