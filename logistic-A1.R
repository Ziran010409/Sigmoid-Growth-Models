
######### understanding the code of inference on logistic model
library(Rtwalk)
library(knitr)

#define the logistic model:(it takes the equation in the logistic html, but divide by e^lambda*t
#on both numerator and denominator, and then take a log)
# Logistic growth
logisode <- function(t, lambda, kappa,log = FALSE){
  ll <-  log(kappa) + log(N0) - log( (kappa-N0)*exp(-lambda*t) + N0)
  if (log)  return(ll)
  else return(exp(ll))
}

# negative log-likelihood function
loglikL <- function(par){
  lambda <- exp(par[1])
  kappa <- exp(par[2])
  sigma <- exp(par[3])
  
  loglik <- sum( dnorm( y - logisode(time, lambda, kappa, log = TRUE)  , 
                        mean = 0, sd = sigma, log = TRUE) )
  
  return(-loglik)
  
}

# define the negative log-likelihood function
#y-logisode(time, lambda, kappa, N0, log = TRUE) is the same as yi−y(λ,κ,h0) in the html
#dnorm means calculate the density of yi−y(λ,κ,h0) using normal distribution with mean=0 and 
#sd=sigma, all these are the same in html



############# prepare dataset
library(growthcurver)
head(growthdata[,c(1,2)],1) #A1 well

plot(growthdata[,c(1,2)])

index <- which(growthdata[,1] >= 6 & growthdata[,1] <= 24)

mass <- as.numeric(as.vector(growthdata[index,2])) #mass here means the absorbance readings
time <- as.numeric(as.vector(growthdata[index,1])) #Time in years
time <- time - min(time)
plot(time,mass)

N0 <- mass[1]

#prepare the data
y <- log(mass)  #since we have taken the log of the likelihood and the ODE before



############## maximum likelihood ##################
OPTL <- nlminb(start = rep(-1,3), objective = loglikL)
#rep(0.1,4), vector length 4 because we have 4 parameter
#nlminb use to find the optimal
MLEL <- exp(OPTL$par) #we have the log-likelihood, so we take exp back
kable(MLEL, digits = 3)
#kable()formats the MLEs in a table, here we require the knitr package

# Fitted the estimates into logistic model(now you have the maximum likelihood estimates, you plug back to the
# logistic model)
fitlogis <- Vectorize(function(t) logisode(t, MLEL[1], MLEL[2], log = FALSE))
#create an anonymous function being vectorized, this function takes a single argument t 
#and calculates the value of the logistic model at time t using the ML estimates 
#vectorization refers to the process of applying a function to each element of a vector 
plot(time, mass, xlab = "time (years)", ylab = "absorbance readings", 
     cex.axis = 1.5, cex.lab = 1.5, pch = 19, ylim = c(0.0,0.5)) #plot the real observations points
curve(fitlogis, 0, max(time), add=TRUE, lwd = 2, col = "blue") 
#plot the curve, to see how the logsitic model with MLE works?

#bic
bic <- 2*loglikL(par = OPTL$par)+2*log(109)





################ Bayesian Inference ##########
#Here we use the Markov chain Monte Carlo (MCMC) method for Bayesian inference

# Analytic solution of Logistic ODE is define at very beginning, the logisode function
# Support 
SupportL <- function(x) {   TRUE }
#The support function specifies the valid region of the parameter space, it always returns TRUE
#which implies that the entire parameter space is considered valid.

# Random initial points
X0L <- function(x) { OPTL$par + runif(3,-0.01,0.01) }
#this generate random initial points base on the MLE, It adds small random noise to MLE

#define the log of the posterior function 
logpostL <- function(par){
  lambda <- exp(par[1])  #in log scale, use exp to transform them back
  kappa <- exp(par[2])
  sigma <- exp(par[3])
  
  loglik <- sum( dnorm( y - logisode(time, lambda, kappa, log = TRUE)  , 
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
infoL <- Runtwalk( dim=3,  Tr=110000,  Obj=logpostL, Supp=SupportL, 
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
colnames(summL) <- c("lambda","kappa","sigma")
kable(summL, digits = 3)
sd <- apply(exp(infoL$output[ind,]),2,sd) 
mean <- apply(infoL$output[ind,],2,mean) 
#bic
bic_bayesian <- 2*loglikL(par = mean)+2*log(109) 


# KDEs of the posterior samples
lambdap <- exp(infoL$output[,1][ind])
#extracts the posterior samples for the parameter lambda from the MCMC output
kappap <- exp(infoL$output[,2][ind])
sigmap <- exp(infoL$output[,3][ind])


plot(density(lambdap), main = "", xlab = expression(lambda), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)

plot(density(kappap), main = "", xlab = expression(kappa), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)



plot(density(sigmap), main = "", xlab = expression(sigma), ylab = "Density",
     cex.axis = 1.5, cex.lab = 1.5, lwd = 2)


######################################################################################
# Posterior envelopes and mean growth curve
######################################################################################

# Creating the posterior envelopes
tvec <- seq(0,18,by = 0.01) #time points, used for plotting the growth curve
ntvec <- length(tvec) #number of time points

# Logistic
NCIL <- matrix(0, ncol = ntvec, nrow = length(ind)) 
#empty matrix with ncol=number of time points, nrow=number of posterior samples

for(j in 1:length(ind)){
  for(k in 1:ntvec){
    NCIL[j,k ] <- logisode( tvec[k],lambdap[j], kappap[j]) 
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
plot(tvec,  NL, type = "l", xlab = "Time", ylab = "absorbance readings", main = "Logistic model",
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



