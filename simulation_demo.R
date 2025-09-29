library('CPAT')
library('cointReg')
library('sde')
library('MASS')
library('stats')
library('pracma')
library('fGarch')

# load confidence_interval_functions.R
source("confidence_interval_functions.R")

######################################
##### data generating process [IID] 
 
dgp_homo<-function(beta, N, delta, thetastar){ 
  y=c(rep(NA,N))
  x=matrix(1,N,2)
  x[,2] = rnorm(N, mean = 0, sd = 1)
  e = rnorm(N, mean = 0, sd = 1)

  kN = floor(thetastar*N)
  y[1:kN] = x[1:kN,]%*%beta+e[1:kN]
  y[(kN+1):N] = x[(kN+1):N,]%*%(beta+delta)+e[(kN+1):N]

  return(list(x = x,y = y))
}


dgp_hetero0<-function(beta, N, delta, thetastar){ 
  y=c(rep(NA,N))
  v = randn(N,1)
  sigma = c(rep(NA,N))
  kN = floor(thetastar*N)

  omegas = c(3, 0.5)
  e = c(rep(NA,N))
  e[1:kN] = rnorm( kN, mean = 0, sd = omegas[1])
  e[(kN+1):N] = rnorm( (N-kN), mean = 0, sd = omegas[2])

  x = matrix(1,N,2)
  x[,2] = rnorm(N, mean = 0, sd = 1)

  y[1:kN] = x[1:kN,]%*%beta+e[1:kN]
  y[(kN+1):N] = x[(kN+1):N,]%*%(beta+delta)+e[(kN+1):N]

  return(list(x = x,y = y,e = e))
}



####### simulation demo
set.seeds(123)

reps = 20 # change the replications to 1000-2000 times
N1 = 300

# parameter setting
h_band=ceiling(sqrt(N1))
kappa = 1/2
alpha = 0.95


# location of change point k_1
thetastar = 0.5 
zeta = sqrt( thetastar*3^2 + (1-thetastar)*0.5^2 ) # for comparable signal-to-noise ratio
# size of change
delta = 3*(N1^(-1/5))*zeta
# model coefficient
beta = as.matrix(c(1,1),2,1)

emp_coverage = matrix(0,reps,1)  
cilength = k_est = matrix(0,reps,1) 
confiband = matrix(NA, reps, 2)

for (i in 1:reps){
    dt = dgp_hetero0(beta, N1, delta, thetastar)
    x = dt$x
    y = dt$y

    khat = confi_21(y,x,kappa)
    kstar = N1*thetastar
    # if(khat>= (N1*0.95)){khat = N1*0.95} # truncation at 0.05
    # if(khat<= (N1*0.05)){khat = N1*0.05} 
    khat <- pmin(pmax(khat, N1 * 0.05), N1 * 0.95) # truncation at 0.05
    deltahat = ols_e(y[1:khat],x[1:khat,]) - ols_e(y[(khat+1):N1],x[(khat+1):N1,])
    
    varl=sigma_step_smooth(y,x,khat, deltahat, h_band)
    a_minus = varl[1]
    a_plus = varl[2]
    
    lim = lim_sim_hete(kappa=kappa, theta=khat/N1, N=20000, nrep=5000, a_minus, a_plus )
    sigmastar = sigma_standard(y,x, khat, deltahat, a_minus, a_plus )
    quan_lim = quantile(lim, c((1-alpha)/2, 1-(1-alpha)/2))
    rev_eng  = quan_lim*sigmastar/( euc_norm(deltahat) )
    rev_eng[2] = max(0, rev_eng[2])
    confi = khat - rev_eng
    confi = sort(confi)

    k_est[i,1] = khat
    confiband[i, 1:2] = c(confi[1],confi[2])
    
    confi[1] <- max(confi[1], 1)
    confi[2] <- min(confi[2], N1)
    cilength[i,1] = confi[2]-confi[1]

    if(kstar>=confi[1] & kstar<=confi[2]){
        emp_coverage[i,1] = 1 
    }
}

# empirical coverage
sum(emp_coverage)/reps
# length of coverage
mean(cilength)/N1
