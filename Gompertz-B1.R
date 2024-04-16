# Required packages and source files
library(Rtwalk)
library(knitr)
ql <- function(p) quantile(p, 0.025)
qu <- function(p) quantile(p, 0.975)


########### data 
library(growthcurver)
head(growthdata[,c("time","B1")],1) #B1 well

plot(growthdata[,c("time","B1")])

index <- which(growthdata[,1] >= 7 & growthdata[,1] <= 24)

mass <- as.numeric(as.vector(growthdata[index,'B1'])) #mass here means the absorbance readings
time <- as.numeric(as.vector(growthdata[index,1])) #Time in years
time <- time - min(time)
plot(time, mass)

N0 <- mass[1]

#prepare the data
y <- log(mass)  #since we have taken the log of the likelihood and the ODE before


########## Gompertz growth
gompertzode <- function(t, lambda, kappa, log = FALSE){
  ll <-  log(kappa) + (log(N0) - log(kappa))*exp(-lambda*t) 
  if (log)  return(ll)
  else return(exp(ll))
}

# negative log-likelihood function
loglikG <- function(par){
  lambda <- exp(par[1])
  kappa <- exp(par[2])
  sigma <- exp(par[3])
  
  loglik <- sum( dnorm( y -gompertzode(time, lambda, kappa, log = TRUE)  , 
                        mean = 0, sd = sigma, log = TRUE) )
  
  return(-loglik)
  
}


######################################################################################
# Maximum likelihood estimation
######################################################################################

# Optimisation step
OPTG <- nlminb(start = rep(-1,3), objective = loglikG)

# Maximum likelihood estimate
MLEG <- exp(OPTG$par)

kable(MLEG, digits = 3)


# Fitted Gompertz model
fitgompertz <- Vectorize(function(t) gompertzode(t, MLEG[1], MLEG[2], log = FALSE))

plot(time, mass, ylim=c(0,0.5),xlab = "time (years)", ylab = "absorbance readings", main = "Gompertz model",
     cex.axis = 1.5, cex.lab = 1.5, pch =19) 
curve(fitgompertz, add=TRUE, lwd = 2, col = "blue")

#calculate BIC

bic <- 2*loglikG(par = OPTG$par)+2*log(103)
## m=2 is the number of parameters estimated in the model, does this m include sigma?
#bic = −2⋅ln(likelihood of the model)+k⋅ln(n)
# but here lgolikeG is negative log likelihood



######################################################################################
# Bayesian estimation
######################################################################################

# negative log-posterior function
logpostG <- function(par){
  lambda <- exp(par[1])
  kappa <- exp(par[2])
  sigma <- exp(par[3])
  
  loglik <- sum( dnorm( y - gompertzode(time, lambda, kappa, log = TRUE)  , 
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
X0G <- function(x) { OPTG$par + runif(3,-0.01,0.01) }

# twalk for analytic solution
set.seed(1234)
infoG <- Runtwalk( dim=3,  Tr=110000,  Obj=logpostG, Supp=SupportG, 
                   x0=X0G(), xp0=X0G(),PlotLogPost = FALSE) 

# Posterior sample after burn-in and thinning
ind=seq(10000,110000,100) 

# Summaries of the posterior samples
summG <- apply(exp(infoG$output[ind,]),2,summary)
colnames(summG) <- c("lambda","kappa","sigma")
kable(summG, digits = 3)
sd <- apply(exp(infoG$output[ind,]),2,sd)
mean <- apply(infoG$output[ind,],2,mean) 
#bic
bic_bayesian <- 2*loglikG(par = mean)+2*log(103) #bic



# KDEs of the posterior samples
lambdap <- exp(infoG$output[,1][ind])
kappap <- exp(infoG$output[,2][ind])
sigmap <- exp(infoG$output[,3][ind])





######################################################################################
# Posterior envelopes and mean growth curve
######################################################################################

# Creating the posterior envelopes
tvec <- seq(0,17,by = 0.01)
ntvec <- length(tvec)

# Logistic
NCIG <- matrix(0, ncol = ntvec, nrow = length(ind))

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    NCIG[j,k ] <- gompertzode( tvec[k],lambdap[j], kappap[j]) 
  }
} 

NG <-  colMeans(NCIG)

NCIGL <- apply(NCIG, 2, ql)
NCIGU <- apply(NCIG, 2, qu)


# Gompertz
plot(tvec,  NG, type = "l", xlab = "Time", ylab = "absorbance readings", main = "Gompertz model",
     cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 2)
points(tvec,  NCIGL, col = "gray", type = "l")
points(tvec,  NCIGU, col = "gray", type = "l")
polygon(c(tvec, rev(tvec)), c(NCIGL[order(tvec)], rev(NCIGU[order(tvec)])),
        col = "gray", border = NA)
points(tvec,  NG,type = "l", col = "black", lwd = 2, lty =2)
points(time, mass, xlab = "time (years)", ylab = "mass", 
       cex.axis = 1.5, cex.lab = 1.5, pch = 19)
