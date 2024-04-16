
######### understanding the code of inference on logistic model
library(Rtwalk)
library(knitr)

#define the logistic model:(it takes the equation in the logistic html, but divide by e^lambda*t
#on both numerator and denominator, and then take a log)
logisode <- function(t, lambda, kappa,N0, log = FALSE){
  ll <-  log(kappa) + log(N0) - log( (kappa-N0)*exp(-lambda*t) + N0)
  if (log)  return(ll)
  else return(exp(ll))
}

# define the negative log-likelihood function
#y-logisode(time, lambda, kappa, N0, log = TRUE) is the same as yi−y(λ,κ,h0) in the html
#dnorm means calculate the density of yi−y(λ,κ,h0) using normal distribution with mean=0 and 
#sd=sigma, all these are the same in html
loglikL <- function(par){
  lambda <- exp(par[1])
  kappa <- exp(par[2])
  N0 <- exp(par[3])
  sigma <- exp(par[4])
  
  loglik <- sum( dnorm( y - logisode(time, lambda, kappa, N0, log = TRUE)  , 
                        mean = 0, sd = sigma, log = TRUE) )
  
  return(-loglik)
  
}



#insert the data and plot it 
mass <- c(6.25,10,20,23,26,27.6,29.8,31.6,37.2,41.2,48.7,54,54,63,66,72,72.2,
          76,75)

days.since.birth <- c(31,62,93,99,107,113,121,127,148,161,180,214,221,307,
                      452,482,923, 955,1308)

plot(days.since.birth, mass, xlab = "time (days)", ylab = "mass", 
     cex.axis = 1.5, cex.lab = 1.5, pch = 19) 


#prepare the data
y <- log(mass)  #since we have taken the log of the likelihood and the ODE before
time <- days.since.birth/365.24 #Time in years


############## maximum likelihood ##################
OPTL <- nlminb(start = rep(0.1,4), objective = loglikL)
#rep(0.1,4), vector length 4 because we have 4 parameter
#nlminb use to find the optimal
MLEL <- exp(OPTL$par) #we have the log-likelihood, so we take exp back
kable(MLEL, digits = 3)
#kable()formats the MLEs in a table, here we require the knitr package

# Fitted the estimates into logistic model(now you have the maximum likelihood estimates, you plug back to the
# logistic model)
fitlogis <- Vectorize(function(t) logisode(t, MLEL[1], MLEL[2], MLEL[3], log = FALSE))
#create an anonymous function being vectorized, this function takes a single argument t 
#and calculates the value of the logistic model at time t using the ML estimates 
#vectorization refers to the process of applying a function to each element of a vector 
plot(time, mass, xlab = "time (years)", ylab = "mass", 
     cex.axis = 1.5, cex.lab = 1.5, pch = 19) #plot the real observations points
curve(fitlogis, 0,4, add=TRUE, lwd = 2, col = "blue") 
#plot the curve, to see how the logsitic model with MLE works?

#bic
bic <- 2*loglikL(par = OPTL$par)+3*log(19)


################ profile likelihood of the parameters (lecture note 3.4, page 7.2)#######
# Required quantities
# Number of parameters (we have 4 unkonwn parameter)
p <- 4
# (minus) Maximum value of the log-likelihood
ML <- OPTL$objective

# Profile likelihood function for parameter "ind"
#par1 is the parameter you are interest in 
#ind is the index of the parameter of interest
prof.likL <- function(par1, ind){
  
  tempf <- function(par){
    tempv <- rep(0,p)
    tempv <- replace(x = tempv, c(1:p)[-ind] , par)
    tempv <- replace(x = tempv, ind , par1)
    out0 <- loglikL(tempv)
    return(out0)
  } 
  
  out <-  -nlminb(OPTL$par[-ind],tempf, control = list(iter.max = 10000))$objective + ML
  
  return(exp(out))
}


prof_indL <- Vectorize(function(par) prof.likL(log(par),indprof)) 
#indprof is the index of the parameter of interest

# Profile likelihood of Parameter 1
indprof <- 1
curve(prof_indL,6,10 , n = 200, lwd = 2, xlab = expression(lambda), ylab = "Profile Likelihood",
      cex.axis = 1.5, cex.lab = 1.5)
#n is the number of points, 6-10 is the range(xlimit),

# Profile likelihood of Parameter 2
indprof <- 2
curve(prof_indL,60,80 , n = 200, lwd = 2, xlab = expression(lambda), ylab = "Profile Likelihood",
      cex.axis = 1.5, cex.lab = 1.5)

# Profile likelihood of Parameter 3
indprof <- 3
curve(prof_indL,2,5 , n = 200, lwd = 2, xlab = expression(kappa), ylab = "Profile Likelihood",
      cex.axis = 1.5, cex.lab = 1.5)

# Profile likelihood of Parameter 4
indprof <- 4
curve(prof_indL,0.02,0.15 , n = 200, lwd = 2, xlab = expression(sigma), ylab = "Profile Likelihood",
      cex.axis = 1.5, cex.lab = 1.5)



################ Bayesian Inference ##########
#Here we use the Markov chain Monte Carlo (MCMC) method for Bayesian inference

# Analytic solution of Logistic ODE is define at very beginning, the logisode function
# Support 
SupportL <- function(x) {   TRUE }
#The support function specifies the valid region of the parameter space, it always returns TRUE
#which implies that the entire parameter space is considered valid.

# Random initial points
X0L <- function(x) { OPTL$par + runif(4,-0.01,0.01) }
#this generate random initial points base on the MLE, It adds small random noise to MLE

#define the log of the posterior function 
logpostL <- function(par){
  lambda <- exp(par[1])  #in log scale, use exp to transform them back
  kappa <- exp(par[2])
  N0 <- exp(par[3])
  sigma <- exp(par[4])
  
  loglik <- sum( dnorm( y - logisode(time, lambda, kappa, N0, log = TRUE)  , 
                        mean = 0, sd = sigma, log = TRUE) ) 
  #this is the log-likelihood, up to here, the code is the same as loglikL() above
  
  #logprior <- dexp(lambda, rate=1) - 
  #           dunif(lambda, min = 0, max = 200) -
  #           dgamma(kappa, shape = 2, scale = 2) -
  #           dunif(sigma, min = 0, max = 1) 
  
  logprior <- 0 
  #this set to 0 meaning that the prior is uniform(constant) across the parameter space
  
  logjacob <- sum(par)
  #it calculates the logarithm of the Jacobian determinant (?)
  
  return(-loglik - logprior - logjacob)
  #return the negative log-posterior
}


# twalk for analytic solution
set.seed(1234)
infoL <- Runtwalk( dim=4,  Tr=110000,  Obj=logpostL, Supp=SupportL, 
                   x0=X0L(), xp0=X0L(),PlotLogPost = FALSE) 
#use the twalk algorithm, this is a type of MCMC 
#Tr is the number of iterations

# Posterior sample after burn-in and thinning
ind=seq(10000,110000,100) 
#After running the MCMC algorithm, it extracts posterior samples from the MCMC chain. 
#It selects every 100th iteration starting from the 10,000th iteration

# Summaries of the posterior samples
summL <- apply(exp(infoL$output[ind,]),2,summary) 
#find summary statistics for each parameter based on the posterior samples
#it provides an insight into the posterior distribution of the parameters.
colnames(summL) <- c("lambda","kappa","h_0","sigma")
kable(summL, digits = 3)
sd <- apply(exp(infoL$output[ind,]),2,sd) 
mean <- apply(infoL$output[ind,],2,mean) 
#bic
bic_bayesian <- 2*loglikL(par = mean)+3*log(19) # BIC= -38.87437


# KDEs of the posterior samples
lambdap <- exp(infoL$output[,1][ind])
#extracts the posterior samples for the parameter lambda from the MCMC output
kappap <- exp(infoL$output[,2][ind])
N0p <- exp(infoL$output[,3][ind])
sigmap <- exp(infoL$output[,4][ind])


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
tvec <- seq(0,4,by = 0.01) #time points, used for plotting the growth curve
ntvec <- length(tvec) #number of time points

# Logistic
NCIL <- matrix(0, ncol = ntvec, nrow = length(ind)) 
#empty matrix with ncol=number of time points, nrow=number of posterior samples

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    NCIL[j,k ] <- logisode( tvec[k],lambdap[j], kappap[j], N0p[j]) 
  }
} 
#filled the emty NCIL matrix with the value of logistic growth model for each 
#combination of posterior samples and time points 

NL <-  colMeans(NCIL) #mean grwoth curve 

#NCILL <- apply(NCIL, 2, ql)
#NCILU <- apply(NCIL, 2, qu)
# "2" means column, ql and qu are functions compute the lower and upper quantiles
#calculates the quantiles of the posterior distribution for each time point
NCILL <- apply(NCIL, 2, function(x) quantile(x, probs = 0.025))
NCILU <- apply(NCIL, 2, function(x) quantile(x, probs = 0.975))



# Logistic
plot(tvec,  NL, type = "l", ylim = c(0,100), xlab = "Time", ylab = "Mass", main = "Logistic model",
     cex.axis = 1.5, cex.lab = 1.5, lwd =2, lty = 2)
points(tvec,  NCILL, col = "gray", type = "l")
points(tvec,  NCILU, col = "gray", type = "l")
#add the upper and lower quantile as the grey lines
polygon(c(tvec, rev(tvec)), c(NCILL[order(tvec)], rev(NCILU[order(tvec)])),
        col = "gray", border = NA)
#filled the area with a gray polygon to represent the posterior envelopes
points(tvec,  NL,type = "l", col = "black", lwd = 2, lty =2)
#Adds the mean growth curve(becasue when you filled with grey polygon, the dash line we plot
# before would be covered, so we plot the line again)
points(time, mass, xlab = "time (years)", ylab = "mass", 
       cex.axis = 1.5, cex.lab = 1.5, pch = 19) 
#Adds the observed data points




  
  
  