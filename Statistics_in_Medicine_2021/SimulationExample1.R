### This simulation study was carried using R version 3.5.0.
### Simulation Study Example: letting prior pi = 0.05, number of feature p = 500, sample size n = 300, 
### number of ordinal levels K = 3, non-zero coefficients |beta| = ln(3).
rm(list=ls())
require("HDInterval")
require("MASS")
require("Matrix")
require("ordinalgmifs")
require("ordinalNet")
require("dclone")
### Prameter set up
set.seed(12345)
a = 0.1 # a and b are hyperparameters in the Gamma prior for lambda: lambda ~ Gamma(a, b)
b = 0.1
N=300 # define sample size
P=500 # define number of covariates
p = 0.05 # calculate prior probability for gamma
Nsim = 30 # define number of replications
Nzero = 10 # define number of nonzero coefficients
epsilon = 0.1 # epsilon is a small positive number for hypothesis testing using bayes factor 
cl = makePSOCKcluster(3) # Create a Parallel Socket Cluster for parallel MCMC chains

### AR1 correlated covariates in the block diagonal
blocksize = 50 # define block size
rho = 0.5 # define correlation between adjacent covariates for AR1 
CovBlock = diag(blocksize) # create block
CovBlock = rho^abs(row(CovBlock)-col(CovBlock)) # build AR1 correlation structure in each block
sigma = bdiag(CovBlock,CovBlock,CovBlock,CovBlock,CovBlock,CovBlock,CovBlock,CovBlock,CovBlock,CovBlock) # this is for P=500, add two extra lines for P=1500

### Define value of coefficients
true.beta = c(rep(0,P))
blockindex = sample(1:blocksize,Nzero,replace=F) # random select the non-zero coefficient at each block
true.index = blockindex + seq(from = 0, to = P-blocksize, by = blocksize) # index of non-sero coefficient when P = 500 
# true.index = blockindex + seq(from = 0, to = P-3*blocksize, by = 3*blocksize) # index of non-sero coefficient when P = 1500 
true.beta[true.index] = c(rep(-log(3), ceiling(Nzero/2)), rep(log(3), Nzero - ceiling(Nzero/2)))

### Proposed Bayesian Penalized Cumulative Logit Model Fixing pi = 0.05
model <- function(){
  for(i in 1:N){
    mu[i] <- inprod(X[i,],betgma[]) 
    logit(Q[i,1]) <- alpha[1]-mu[i]
    p[i,1] <- Q[i,1]
    for(j in 2:(k-1)){
      logit(Q[i,j]) <- alpha[j]-mu[i]
      p[i,j] <- Q[i,j] - Q[i,j-1]
    }
    p[i,k] <- 1 - Q[i,(k-1)]
    y[i,] ~ dcat(p[i,1:k]) 
  }
  for (b in 1:P){ 
    beta[b] ~ ddexp(0,lambda)
  }             
  lambda ~ dgamma(0.1,0.1)
  
  for (j in 1:P){
    betgma[j] <- beta[j] * gamma[j]
    gamma[j] ~ dbern(0.05)            
  }
  for (r in 1:(k-1)) {
    alpha0[r] ~ dnorm(0, 1.0E-1)
  }
  alpha <- sort(alpha0)
}

### JAGS Set-up
adaptSteps <- 5000              #number of steps to "tune" the samplers
burnInSteps <- 5000             #number of steps to "burn-in" the samplers
nChains <- 3                    #number of chains to run
numSavedSteps <- 9999           #total number of steps in chains to save
thinSteps <- 3                 	#number of steps to "thin" (1=keep every step)
nIter <- ceiling((numSavedSteps*thinSteps )/nChains) 	#steps per chain
set.seed(10001)

### Initialize vector/matrix/array to save results
betgmamatrix = array(NA, c(numSavedSteps, P, Nsim))
gammamatrix = array(NA, c(numSavedSteps, P, Nsim))
gamma.TP = numeric()
gamma.FP = numeric()
gammaBF.TP = numeric()
gammaBF.FP = numeric()
bgTP = numeric()
bgFP = numeric()
HPDbg.TP = numeric()
HPDbg.FP = numeric()
ETbg.TP = numeric()
ETbg.FP = numeric()
gmifsAIC.TP = numeric()
gmifsAIC.FP = numeric()
gmifsBIC.TP = numeric()
gmifsBIC.FP = numeric()
OrdinalNetAIC.TP = numeric()
OrdinalNetAIC.FP = numeric()
OrdinalNetBIC.TP = numeric()
OrdinalNetBIC.FP = numeric()
### For loop for Nsim replications
for (j in 1:Nsim){
  ### generate & standardize data 
  X = mvrnorm(N,mu=rep(0,P),Sigma=sigma)
  meanx = apply(X, 2, mean)
  sdx = apply(X, 2, sd)
  for (l in 1:P) {
    for (k in 1:N) {
      X[k,l] <- (X[k,l]-meanx[l])/sdx[l]
    }
  }
  error = rlogis(N, location = 0, scale = 1)
  z <- X%*%true.beta + error
  y <- z
  quantile(z, probs = c(1/3, 2/3))
  y[z < quantile(z,1/3)] <- 1
  y[z >= quantile(z,1/3) & z < quantile(z,2/3)] <- 2
  y[z >= quantile(z,2/3)] <- 3
  k<-length(unique(y))
  
  #initialize and adapt the model
  jagsdf <- list(y=y,X=X, N=dim(X)[1], P =dim(X)[2], k=k) 
  params <- c("gamma","betgma")
  inits1 <- list("alpha0" = sort(c(rnorm(1,log(1/2),.5), rnorm(1,log(2/1),.5))), "beta" = rep(0,P), "pi" = rep(0,P))
  inits2 <- list("alpha0" = c(log(1/2), log(2/1)), "beta" = rep(0,P), "pi" = rep(0,P))
  inits3 <- list("alpha0" = sort(c(rnorm(1,log(1/2),.5), rnorm(1,log(2/1),.5))), "beta" = rep(0,P), "pi" = rep(0,P))
  
  parJagsModel(cl, name = "test", file = model, 
               data = jagsdf, inits = list(inits1, inits2, inits3),
               n.chains = nChains, n.adapt=adaptSteps, quiet=FALSE)
  # Burn-in the algorithm
  print(paste0("Replication Number: ", j))
  cat( "Burning in the MCMC chain...\n")
  parUpdate(cl, "test", n.iter=burnInSteps)
  # Run MCMC algorithm
  cat("Sampling final MCMC chain...\n" )
  codaSamples <- parCodaSamples(cl, "test", variable.names=params, n.iter=nIter, thin=thinSteps)
  
  ### Retrieve posterior samples
  mcmcChain <- as.matrix(codaSamples)
  betgmamatrix[ , ,j] = mcmcChain[,grep("betgma",dimnames(mcmcChain)[[2]])]
  gammamatrix[ , ,j] = mcmcChain[,grep("gamma",dimnames(mcmcChain)[[2]])]
  
  ### select variable using top 10 gamma mean's
  gammamean = apply(gammamatrix[ , ,j], 2,mean)
  gamma10.index = sort(gammamean, decreasing = T, index.return =T)$ix[1:10]
  gamma.TP[j] = length(intersect(gamma10.index, true.index))
  gamma.FP[j] = length(gamma10.index) - length(intersect(gamma10.index, true.index))
  
  ### select variable using bayes factor based on gamma
  gammasum = apply(gammamatrix[ , ,j], 2,sum)
  PostOdds = gammasum/(dim(gammamatrix)[1]-gammasum)
  PriorOdds = p/(1-p)
  gamma.index = (1:P)[PostOdds/PriorOdds > 5] # let threshold = 5 for bayes factor
  if(length(gamma.index)>0){
    gammaBF.TP[j] = length(intersect(gamma.index, true.index))
    gammaBF.FP[j] = length(gamma.index) - length(intersect(gamma.index, true.index))
  }
  
  ### select variable using bayes factor based on gammaXbeta
  BetgmaPost1 = apply(betgmamatrix[ , ,j], 2, function(x) sum(abs(x) > epsilon))
  BetgmaPost2 = apply(betgmamatrix[ , ,j], 2, function(x) sum(abs(x) <= epsilon))
  PostOdds = BetgmaPost1/BetgmaPost2
  PriorOdds = (p*b^a*gamma(a+1))/(p*gamma(a+1)*(b+epsilon)^a-p*gamma(a+1)*b^a+(1-p)*a*gamma(a)*(b+epsilon)^a)
  BetgmaFactor = PostOdds/PriorOdds # calculate bayes factor
  BFbeta = (1:P)[BetgmaFactor>5] # let threshold = 5 for bayes factor
  if(length(BFbeta)>0){
    bgTP[j] = length(intersect(BFbeta, true.index)) 
    bgFP[j] = length(BFbeta) - length(intersect(BFbeta, true.index))
  }
  
  ### select variable using 95% equal-tailed credible interval of gammaXbeta
  e = (1:P)[!((apply(betgmamatrix[ , ,j],2,quantile,0.025) <= 0) & (0 <= apply(betgmamatrix[ , ,j],2,quantile,0.975)))]
  if(length(e)>0){
    ETbg.TP[j] = length(intersect(e, true.index))
    ETbg.FP[j] = length(e) - length(intersect(e, true.index))
  }
  
  ### select variable using 95% HPD credible interval of gammaXbeta
  e = (1:P)[!((apply(betgmamatrix[ , ,j],2,hdi,0.95)[1,] <= 0) & (0 <= apply(betgmamatrix[ , ,j],2,hdi,0.95)[2,]))]
  if(length(e)>0){
    HPDbg.TP[j] = length(intersect(e, true.index))
    HPDbg.FP[j] = length(e) - length(intersect(e, true.index))
  }  
  
  ##### ordinalgmifs method 
  mydata = as.data.frame(cbind(y, X))
  colnames(mydata) = c("y", sprintf("V%s",seq(1:P)))
  mydata$y = factor(mydata$y, ordered = T)
  fitGMIFS = ordinalgmifs(y ~ 1, x = mydata[,-1], data = mydata, tol=0.01, probability.model = "Cumulative", link = "logit")
  GMIFS.AIC = which(coef(fitGMIFS, model.select = "AIC")!=0)[-c(1:2)] - 2 # index for variables selected using ordinalgmifs with minimum AIC
  gmifsAIC.TP[j] = length(intersect(GMIFS.AIC, true.index)) # Number of true positives using ordinalgmifs with minimum AIC
  gmifsAIC.FP[j] = length(GMIFS.AIC) - length(intersect(GMIFS.AIC, true.index)) # Number of false positives using ordinalgmifs with minimum AIC
  GMIFS.BIC = which(coef(fitGMIFS, model.select = "BIC")!=0)[-c(1:2)] - 2 # index for variables selected using ordinalgmifs with minimum BIC
  gmifsBIC.TP[j] = length(intersect(GMIFS.BIC, true.index))  # Number of true positives using using ordinalgmifs with minimum BIC
  gmifsBIC.FP[j] = length(GMIFS.BIC) - length(intersect(GMIFS.BIC, true.index))  # Number of false positives using ordinalgmifs with minimum BIC
  
  ##### ordinalNet method 
  y = factor(y, ordered = T)
  FitOrdinalNet <- ordinalNet(X, y, family="cumulative", link="logit",
                              parallelTerms=T, nonparallelTerms=F)
  OrdinalNet.AIC = as.numeric(which(coef(FitOrdinalNet, criteria = "aic")!=0))[-c(1:2)] - 2 # index for variables selected using ordinalNet with minimum AIC
  OrdinalNetAIC.TP[j] = length(intersect(OrdinalNet.AIC, true.index))  # Number of true positives using ordinalNet with minimum AIC
  OrdinalNetAIC.FP[j] = length(OrdinalNet.AIC) - length(intersect(OrdinalNet.AIC, true.index)) # Number of false positives using ordinalNet with minimum AIC   
  OrdinalNet.BIC = as.numeric(which(coef(FitOrdinalNet, criteria = "bic")!=0))[-c(1:2)] - 2 # index for variables selected using ordinalNet with minimum BIC
  OrdinalNetBIC.TP[j] = length(intersect(OrdinalNet.BIC, true.index)) # Number of true positives using using ordinalNet with minimum BIC
  OrdinalNetBIC.FP[j] = length(OrdinalNet.BIC) - length(intersect(OrdinalNet.BIC, true.index))  # Number of false positives using ordinalNet with minimum BIC
}
save.image("SimulationExample1.RData")

