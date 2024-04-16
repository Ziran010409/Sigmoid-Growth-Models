library(Rtwalk)
library(knitr)
library(deSolve)

############# prepare dataset
library(growthcurver)
head(growthdata[,c("time","B1")],1) #B1 well

plot(growthdata[,c("time","B1")])

index <- which(growthdata[,1] >= 7 & growthdata[,1] <= 24)

mass <- as.numeric(as.vector(growthdata[index,'B1'])) #mass here means the absorbance readings
time <- as.numeric(as.vector(growthdata[index,1])) #Time in years
time <- time - min(time)

plot(time, mass)
N0 <- mass[1]

y <- log(mass) 


# define the ODE of Gompertz-Richards model, t is necessary though not used in the calculation
gompertz_richardsode <- function(t, y, parms){
  with(as.list(c(y,parms)),{
    dN <- lambda* N *((1 - (N / kappa) ^ beta) * log(kappa / N))
    return(list(dN))
  })
}


# Define the negative log-likelihood function for Richards' growth model
loglikGR <- function(par){
  lambda <- exp(par[1])
  kappa <- exp(par[2])
  beta <- exp(par[3])
  sigma <- exp(par[4])
  
  #find numerical solutions
  num_sol <- ode(y=c(N=mass[1]),times = time, func = gompertz_richardsode, method = "lsode",
                 parms = c(lambda=lambda,kappa=kappa,beta=beta))[,"N"]
  
  #define log-likelihood
  loglik <- sum( dnorm( y - log(num_sol), 
                        mean = 0, sd = sigma, log = TRUE) )
  
  return(-loglik)
}

###################################
# Maximum Likelihood Estimation (MLE)
################################### 
OPTGR <- nlminb(start = rep(-0.2,4), objective = loglikGR)
MLEGR <- exp(OPTGR$par)
kable(MLEGR, digits = 3)

fitGR <- Vectorize(function(t) {
  sol <- ode(y = c(N = mass[1]), times = c(0, t), func = gompertz_richardsode, 
             parms = c(lambda = MLEGR[1], kappa = MLEGR[2], beta = MLEGR[3], N0 = mass[1]))
  N_t <- sol[2, "N"] 
  return(N_t)
})

plot(time, mass, xlab = "time (years)", ylab = "absorbance readings", 
     main = "Gompertz-Richards' model",
     cex.axis = 1.5, cex.lab = 1.5, pch =19) 
curve(fitGR, add=TRUE,lwd = 2, col = "blue")


bic <- 2*loglikGR(par = OPTGR$par)+2*log(103) # bic





####################################
#  Bayesian Inference for Gompertz-Richards' growth model
####################################

###################
SupportGR <- function(x) {TRUE }
X0GR <- function(x) { OPTGR$par + runif(4,-0.1,0.1) }

#posterior function
logpostGR <- function(par){
  lambda <- exp(par[1])
  kappa <- exp(par[2])
  beta <- exp(par[3])
  sigma <- exp(par[4])
  
  #find numerical solutions
  num_sol <- ode(y=c(N=mass[1]),times = time, func = gompertz_richardsode, method = "lsode",
                 parms = c(lambda=lambda,kappa=kappa,beta=beta))[,"N"]
  
  #define log-likelihood
  loglik <- sum( dnorm( y - log(num_sol), 
                        mean = 0, sd = sigma, log = TRUE) )
  
  logprior <- sum(dnorm(par,0,1,log=TRUE)) #use N(0,1) prior
  
  logjacob <- sum(par)
  
  return(-loglik - logprior - logjacob)
}

# twalk 
set.seed(1234)
infoGR <- Runtwalk( dim=4,  Tr=110000,  Obj=logpostGR, Supp=SupportGR, 
                    x0=X0GR(), xp0=X0GR(),PlotLogPost = TRUE) 


# Posterior sample after burn-in and thinning
ind=seq(10000,110000,100) 

# Summaries of the posterior samples
summGR <- apply(exp(infoGR$output[ind,]),2,summary)
colnames(summGR) <- c("lambda","kappa","beta","sigma")
kable(summGR, digits = 3)
sd <- apply(exp(infoGR$output[ind,]),2,sd)
mean <- apply(infoGR$output[ind,],2,mean)
#bic
bic_bayesian <- 2*loglikGR(par = mean)+2*log(103)



lambdap <- exp(infoGR$output[,1][ind])
kappap <- exp(infoGR$output[,2][ind])
betap <- exp(infoGR$output[,3][ind])
sigmap <- exp(infoGR$output[,4][ind])


#################################
# Posterior envelopes and mean growth curve
#################################

# Creating the posterior envelopes
tvec <- seq(0,17,by = 0.01) #time points, used for plotting the growth curve
ntvec <- length(tvec) #number of time points

# Gompertz-Richards' model
NCIGR <- matrix(0, ncol = ntvec, nrow = length(ind))  # 151*401 matrix
#empty matrix with ncol=number of time points, nrow=number of posterior samples

#It takes quite a long time to run because you need to run the ODE solver so many times 
for(j in 1:length(ind)){
  for(k in 1:ntvec){
    NCIGR[j,k] <- ode(y=c(N=mass[1]),times =tvec, func = gompertz_richardsode, method = "lsode",
                      parms = c(lambda=lambdap[j],kappa=kappap[j],beta=betap[j]))[,"N"][k]
  }
} 
#filled the empty NCIGR matrix with the value of Gompertz-Richards' growth model for each 
#combination of posterior samples and time points 

NGR <-  colMeans(NCIGR) #mean grwoth curve 

#calculates the quantiles of the posterior distribution for each time point
NCIGRL <- apply(NCIGR, 2, function(x) quantile(x, probs = 0.025))
NCIGRU <- apply(NCIGR, 2, function(x) quantile(x, probs = 0.975))

#Gompertz-Richards' model
plot(tvec,  NGR, type = "l", xlab = "Time", ylab = "absorbance readings", 
     main = "Gompertz-Richards model",cex.main=0.5,
     cex.axis = 0.5, cex.lab = 0.5, lwd =1, lty = 1)
points(tvec,  NCIGRL, col = "gray", type = "l")
points(tvec,  NCIGRU, col = "gray", type = "l")
#add the upper and lower quantile as the grey lines
polygon(c(tvec, rev(tvec)),c(NCIGRL[order(tvec)], rev(NCIGRU[order(tvec)])),
        col = "gray", border = NA)
#filled the area with a gray polygon to represent the posterior envelopes
points(tvec,  NGR,type = "l", col = "black", lwd = 2, lty =2)
#Adds the mean growth curve(becasue when you filled with grey polygon, the dash line we plot
# before would be covered, so we plot the line again)
points(time, mass, xlab = "time (years)", ylab = "mass", 
       cex.axis = 1.5, cex.lab = 1.5, pch = 19) 
#Adds the observed data points

