rm(list =ls())

# Required packages and source files
library(Rtwalk)
library(knitr)
ql <- function(p) quantile(p, 0.025)
qu <- function(p) quantile(p, 0.975)

######################################################################################
# Original Data: Wilson's growth data
# https://bscheng.com/2014/05/07/modeling-logistic-growth-data-in-r/
######################################################################################

mass <- c(6.25,10,20,23,26,27.6,29.8,31.6,37.2,41.2,48.7,54,54,63,66,72,72.2,
          76,75) #Wilson's mass in pounds
days.since.birth <- c(31,62,93,99,107,113,121,127,148,161,180,214,221,307,
                      452,482,923, 955,1308) #days since Wilson's birth

plot(days.since.birth, mass, xlab = "time (days)", ylab = "mass", 
     cex.axis = 1.5, cex.lab = 1.5, pch = 19) 

######################################################################################
# Data preparation
######################################################################################

# log mass
y <- log(mass)
# Time in years
time <- days.since.birth/365.24


######################################################################################
# Gompertz model
######################################################################################

# Gompertz growth
gompertzode <- function(t, lambda, kappa,N0, log = FALSE){
  ll <-  log(kappa) + (log(N0) - log(kappa))*exp(-lambda*t) 
  if (log)  return(ll)
  else return(exp(ll))
}

# negative log-likelihood function
loglikG <- function(par){
  lambda <- exp(par[1])
  kappa <- exp(par[2])
  N0 <- exp(par[3])
  sigma <- exp(par[4])
  
  loglik <- sum( dnorm( y -gompertzode(time, lambda, kappa, N0, log = TRUE)  , 
                        mean = 0, sd = sigma, log = TRUE) )
  
  return(-loglik)
  
}


######################################################################################
# Maximum likelihood estimation
######################################################################################

# Optimisation step
OPTG <- nlminb(start = rep(0.1,4), objective = loglikG)

# Maximum likelihood estimate
MLEG <- exp(OPTG$par)

kable(MLEG, digits = 3)


# Fitted Gompertz model
fitgompertz <- Vectorize(function(t) gompertzode(t, MLEG[1], MLEG[2], MLEG[3], log = FALSE))

plot(time, mass, xlab = "time (years)", ylab = "mass", main = "Gompertz model",
     cex.axis = 1.5, cex.lab = 1.5, pch =19) 
curve(fitgompertz, 0,4, add=TRUE, lwd = 2, col = "blue")

#calculate BIC

bic <- 2*loglikG(par = OPTG$par)+3*log(19)
## m=3 is the number of parameters estimated in the model, does this m include sigma?
#bic = −2⋅ln(likelihood of the model)+k⋅ln(n)
# but here lgolikeG is negative log likelihood


##################################################################################################
# Profile likelihood of the parameters
# See Chapter 3 of https://github.com/FJRubio67/StatisticalInference
##################################################################################################
# Proile likelihood function for parameter "ind"
prof.likG <- function(par1, ind){
  
  tempf <- function(par){
    tempv <- rep(0,p)
    tempv <- replace(x = tempv, c(1:p)[-ind] , par)
    tempv <- replace(x = tempv, ind , par1)
    out0 <- loglikG(tempv)
    return(out0)
  } 
  
  out <-  -nlminb(OPTG$par[-ind],tempf, control = list(iter.max = 10000))$objective + MG
  
  return(exp(out))
}


prof_indG <- Vectorize(function(par) prof.likG(log(par),indprof))


# Required quantities
# Number of parameters
p <- 4
# (minus) Maximum value of the log-likelihood
MG <- OPTG$objective


# Profile likelihoods

# Profile likelihood of Parameter 1
indprof <- 1
curve(prof_indG,3.5,5 , n = 200, lwd = 2, xlab = expression(lambda), ylab = "Profile Likelihood",
      cex.axis = 1.5, cex.lab = 1.5)

# Profile likelihood of Parameter 2
indprof <- 2
curve(prof_indG,65,80 , n = 200, lwd = 2, xlab = expression(lambda), ylab = "Profile Likelihood",
      cex.axis = 1.5, cex.lab = 1.5)

# Profile likelihood of Parameter 3
indprof <- 3
curve(prof_indG,1,3, n = 200, lwd = 2, xlab = expression(kappa), ylab = "Profile Likelihood",
      cex.axis = 1.5, cex.lab = 1.5)

# Profile likelihood of Parameter 4
indprof <- 4
curve(prof_indG,0.02,0.15 , n = 200, lwd = 2, xlab = expression(sigma), ylab = "Profile Likelihood",
      cex.axis = 1.5, cex.lab = 1.5)


######################################################################################
# Bayesian estimation
######################################################################################

# negative log-posterior function
logpostG <- function(par){
  lambda <- exp(par[1])
  kappa <- exp(par[2])
  N0 <- exp(par[3])
  sigma <- exp(par[4])
  
  loglik <- sum( dnorm( y - gompertzode(time, lambda, kappa, N0, log = TRUE)  , 
                        mean = 0, sd = sigma, log = TRUE) )
  
  # logprior <- dgamma(lambda, shape = 2, scale = 2) - 
  #             dgamma(kappa, shape = 2, scale = 2) -
  #             dgamma(N0, shape = 2, scale = 2) -
  #             dgamma(sigma, shape = 2, scale = 2) 
  
  logprior <- 0
  
  logjacob <- sum(par)
  
  return(-loglik - logprior - logjacob)
  
}

# Support
SupportG <- function(x) {   TRUE }

# Random initial points
X0G <- function(x) { OPTG$par + runif(4,-0.01,0.01) }

# twalk for analytic solution
set.seed(1234)
infoG <- Runtwalk( dim=4,  Tr=110000,  Obj=logpostG, Supp=SupportG, 
                   x0=X0G(), xp0=X0G(),PlotLogPost = FALSE) 

# Posterior sample after burn-in and thinning
ind=seq(10000,110000,100) 

# Summaries of the posterior samples
summG <- apply(exp(infoG$output[ind,]),2,summary)
colnames(summG) <- c("lambda","kappa","h_0","sigma")
kable(summG, digits = 3)
sd <- apply(exp(infoG$output[ind,]),2,sd)
mean <- apply(infoG$output[ind,],2,mean) 
#bic
bic_bayesian <- 2*loglikG(par = mean)+3*log(19) #bic= -43.76159



# KDEs of the posterior samples
lambdap <- exp(infoG$output[,1][ind])
kappap <- exp(infoG$output[,2][ind])
N0p <- exp(infoG$output[,3][ind])
sigmap <- exp(infoG$output[,4][ind])




plot(density(lambdap), main = "", xlab = expression(lambda), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
plot(density(kappap), main = "", xlab = expression(kappa), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
plot(density(N0p), main = "", xlab = expression(N[0]), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)
plot(density(sigmap), main = "", xlab = expression(sigma), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)

######################################################################################
# Posterior envelopes and mean growth curve
######################################################################################

# Creating the posterior envelopes
tvec <- seq(0,4,by = 0.01)
ntvec <- length(tvec)

# Logistic
NCIG <- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    NCIG[j,k ] <- gompertzode( tvec[k],lambdap[j], kappap[j], N0p[j]) 
  }
} 

NG <-  colMeans(NCIG)

NCIGL <- apply(NCIG, 2, ql)
NCIGU <- apply(NCIG, 2, qu)


# Gompertz
plot(tvec,  NG, type = "l", ylim = c(0,100), xlab = "Time", ylab = "Mass", main = "Gompertz model",
     cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 2)
points(tvec,  NCIGL, col = "gray", type = "l")
points(tvec,  NCIGU, col = "gray", type = "l")
polygon(c(tvec, rev(tvec)), c(NCIGL[order(tvec)], rev(NCIGU[order(tvec)])),
        col = "gray", border = NA)
points(tvec,  NG,type = "l", col = "black", lwd = 2, lty =2)
points(time, mass, xlab = "time (years)", ylab = "mass", 
       cex.axis = 1.5, cex.lab = 1.5, pch = 19)




