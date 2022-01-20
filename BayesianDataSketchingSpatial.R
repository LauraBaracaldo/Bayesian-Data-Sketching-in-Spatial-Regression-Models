######################################################################
################ Data analysis with static and varying coefficients

###Packages
library(mvtnorm)
library(Rcpp)
library(Matrix)
library(splines)
sourceCpp("test.cpp")  

############################################
###### Data set
############################################
Data<- data_cleaned
N <-n <- dim(Data)[1]
Y<- log(Data.train[, c("NDVI")] + 1)   #Response
X<- cbind(rep(1, n), Data.train$Non_Vegetated_or_Builtup_Lands )  #Covariate
X_D <- t(bdiag(split(X,1:nrow(X))))
S <- cbind(Data.train$scaled_x, Data.train$scaled_y) ##Coordinates
p<- dim(X)[2]


#####################################################################
### (1) Geostatistical sketching
#####################################################################
m<- 10*ceiling(sqrt(n))
Phi<- matrix(rnorm(m*n, 0,  sqrt(1/n)),m,n)  ##Compressing matrix
Y.tilde <- eigenMatMult(Phi, as.matrix(Y)) 
X_D.tilde <- Phi%*%X_D
X.tilde <- Phi%*%X


#####################################################################
### (2)  B-splines (Tensor product of marginal B-splines)
#####################################################################

# Create B-splines
deg   = 3  ##Polynomial degree
N  <-  n
ndx   =  (100 *N^(1/5))^(1/2)
#ndx= 12
xr_1    = max(S[,1])
xl_1    = min(S[,1])
xmax_1  = xr_1 + 0.01 * (xr_1 - xl_1)
xmin_1  = xl_1 - 0.01 * (xr_1 - xl_1)
dt_1    = (xmax_1 - xmin_1) / ndx
knots_1 = seq(xmin_1 - deg * dt_1, xmax_1 + deg * dt_1, by = dt_1)

xr_2    = max(S[,2])
xl_2    = min(S[,2])
xmax_2  = xr_2 + 0.01 * (xr_2 - xl_2)
xmin_2  = xl_2 - 0.01 * (xr_2 - xl_2)
dt_2    = (xmax_2 - xmin_2) / ndx
knots_2 = seq(xmin_2 - deg * dt_2, xmax_2 + deg * dt_2, by = dt_2)

B1     = splineDesign(knots = knots_1, x = S[,1], ord = deg + 1, derivs = 0,outer.ok = TRUE)
B2     = splineDesign(knots = knots_2, x = S[,2], ord = deg + 1, derivs = 0,outer.ok = TRUE)

#### Tensor product of marginal b-splines
B = matrix(nrow = N, ncol = dim(B1)[2]*dim(B2)[2])
for (i in 1:(dim(B)[1])) 
{B[i,] = kronecker(B1[i,], B2[i,])}
B<- Matrix(B, sparse = TRUE) 
diagp<-Matrix(diag((p+1)), sparse=TRUE) 
####Matrix B of tensor product of B-splines 

BB<-kronecker(diagp ,t(B[1,]))
for(k in 2:length(Y))
{BB <- rbind(BB, kronecker(diagp ,t(B[k,])) )
  print(k)}


XB<- X_D.tilde%*%BB
Y.tilde<- as.matrix(Y.tilde)
XB<- Matrix(XB, sparse = TRUE) 
X.tilde<- Matrix(X.tilde, sparse = TRUE) 
X<- Matrix(X, sparse = TRUE) 
#####################################################################
###### Full conditionals
#####################################################################
sample.Gamma<- function(XB, sigma2, tau2, Y.tilde, H, p,m, X.tilde, beta){
  
  
  Delta<-   bdiag((tau2[1])*diag(H))
  for(k in 2:(p+1))
  {
    Delta<- bdiag(Delta, (tau2[k])*diag(H))
  }
  
  g1<- rmvnormCPP(rep(0, H*(p+1)), as.matrix(Delta))
  g2 <- rmvnormCPP(rep(0, m), diag(m))
  g3<- (1/(sqrt(sigma2)))*eigenMatMult(XB,g1) + g2
  g4 <- eigenMatMult(chol2inv(chol((1/sigma2)*eigenMatMult(XB, eigenMatMult(as.matrix(Delta), t(XB)))  + diag(m))) , ((Y.tilde - X.tilde%*%beta) /sqrt(sigma2) - g3))
  Gamma <- g1 + (1/sqrt(sigma2))*eigenMatMult(as.matrix(Delta), eigenMatMult(t(XB), g4))
  
  
  return(Gamma)
}
sample.beta<- function(X.tilde, sigma2, Y.tilde, XB, Gamma, sig2.beta, p, m){
  
  Sigma<- chol2inv(chol( (1/sigma2)*eigenMatMult(t(X.tilde), X.tilde)     + (1/sig2.beta)*diag(p+1)))
  mu<-  (1/sigma2)*eigenMatMult(as.matrix(Sigma), eigenMatMult(t(X.tilde),   (Y.tilde - XB%*%Gamma) )) 
  beta<- rmvnormCPP(mu, as.matrix(Sigma))
  return(beta)
}
sample.sig2<- function(a.sig, b.sig, M, Y.tilde, XB, Gamma, X.tilde, beta){
  YY <- Y.tilde - eigenMatMult(as.matrix(XB), as.matrix(Gamma)) - eigenMatMult(as.matrix(X.tilde), as.matrix(beta))
  sigma2<- 1/(rgamma(1,a.sig + M/2, b.sig + 1/2*(eigenMatMult( t(YY), YY))))
  return(sigma2)}
sample.tau2 <- function(a.tau, b.tau, H, Gamma, p){ 
  tau2<- rep(0, (p+1))
  for(h in 1:(1+p))
  { Gamma_h <- as.matrix(Gamma[((h-1)*H +1):(h*H)])
  tau2[h]<- 1/(rgamma(1, a.tau + H/2,b.tau + 1/2*eigenMatMult(t(Gamma_h), Gamma_h)))}
  return(tau2)
}

#####################################################################
## (3) MCMC algorithm
#####################################################################
H<- (dim(BB)[2])/(p+1)
M<- length(Y.tilde)

##Initial values
Gamma<- rnorm((p+1)*H, 0, 1)
sigma2<- 1/rgamma(1, 1, 1)
tau2<- 1/rgamma((p+1), 1, 1)
beta<- rnorm(p+1)




##Hyperparameters
a.tau <- 2
b.tau <- 0.1
a.sig<- 2
b.sig<- 0.1
sig2.beta<- 100

# Setup number of runs and variables where output will be stored.
maxit <- 2000
burn<-1000
Gamma.out <- matrix(NA,maxit-burn,(H*(p+1)))
tau2.out <- matrix(NA,maxit-burn,(p+1))
sig2.out<-  matrix(NA,maxit-burn,1)
beta.out <- matrix(NA,maxit-burn,((p+1)))

# MCMC loop


start_time <- Sys.time()
for(s in 1:maxit){
  
  Gamma    <- sample.Gamma(XB, sigma2, tau2, Y.tilde, H, p,M, X.tilde, beta)  
  beta<- sample.beta(X.tilde, sigma2, Y.tilde, XB, Gamma, sig2.beta, p, m)
  sigma2<- sample.sig2(a.sig, b.sig, M, Y.tilde, XB, Gamma, X.tilde, beta)
  tau2<- sample.tau2(a.tau, b.tau, H, Gamma, p)
  
  
  
  print(s)
  
  
  if((s/500)==floor(s/500)){
    print(s)
  }
  if(s>burn){
    Gamma.out[s-burn,] <- Gamma
    beta.out[s-burn,] <- beta
    tau2.out[s-burn,] <- tau2
    sig2.out[s-burn] <- sigma2
    
    
  }
}
end_time <- Sys.time()
time1<- end_time-start_time

#####################################################################
#### Estimation of varying coeff.and prediction of response variable

BB<-Matrix(BB, sparse = TRUE) 
y.hat <- array(NA, dim=c(maxit-burn, N))
for(s in 1:(maxit-burn))
{
  Beta.hat.tem <- (BB%*%Gamma.out[s,])
  mean.temp<- as.vector(X_D%*%Beta.hat.tem + X%*%beta.out[s,])
  y.hat[s,] <- rnorm(n,  mean.temp , sqrt(sig2.out[s]) )
  print(s)
}

#####################################################################
#### (4) Prediction of response and varying coefficients
#####################################################################

####################################################
## B-splines basis functions on new sample
####################################################
B1.new     = splineDesign(knots = knots_1, x = S.new[,1], ord = deg + 1, derivs = 0,outer.ok = TRUE)
B2.new     = splineDesign(knots = knots_2, x = S.new[,2], ord = deg + 1, derivs = 0,outer.ok = TRUE)


#### Tensor product of marginal b-splines
B.new = matrix(nrow = n.new, ncol = dim(B1.new)[2]*dim(B2.new)[2])
for (i in 1:(dim(B.new)[1])) 
{
  B.new[i,] = kronecker(B1.new[i,], B2.new[i,])
}



####Matrix B of tensor product of B-splines
B.new<-Matrix(B.new, sparse = TRUE) 
p<- 1
diagp<-Matrix(diag((p+1)), sparse=TRUE) 


BB.new<-kronecker(diagp,t(B.new[1,]))
for(k in 2:length(Y.new))
{
  BB.new <- rbind(BB.new, kronecker(diagp,t(B.new[k,])) )
  print(k)
}

#####################################################################
#### Sampling from predictive distribution
#####################################################################
BB.new<-Matrix(BB.new, sparse = TRUE) 
Y.new.pred<- matrix(NA, maxit-burn, n.new)
Beta.new.hat <- matrix(NA, maxit-burn, n.new*(p+1) )

for(s in 1:(maxit-burn)){
   B.n<- as.matrix(BB.new%*%Gamma.out[s,])
  mu <- X_D.new%*%B.n + X.new%*%beta.out[s,]
  Y.new.pred[s,]<- rnorm(n.new,  as.vector(mu), sqrt(sig2.out[s]) )
  print(s)
  }












