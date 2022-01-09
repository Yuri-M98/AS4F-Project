
#Dataset and Time Series
require(quantmod)
getSymbols("AMT", from="2016-01-01", to="2021-07-21")
colnames(AMT)
attr(AMT, "src")
plot(AMT, main=" American Tower Inc. Stock", col = "red" )
S <- AMT [,"AMT.Close"]
lineChart(S, up.col= "Black", name="AMT Closing Price", theme = "white")
X <- diff(log(S))



### Change point analysis (from 2016 to 2021)
library(tseries)
S <- get.hist.quote("AMT", start = "2016-01-01", end = "2021-07-21")
chartSeries(S, up.col= "Black",TA=c(addVo(), addBBands()), theme="white")
S <- S$Close
install.packages("sde")
require(sde)
cpoint(S)
addVLine = function(dtlist) plot(addTA(xts(rep(TRUE,NROW(dtlist)),dtlist),on=1,col="red"))
addVLine(cpoint(S)$tau0)

## all log-returns with change-point line
S <- as.numeric(S)
n <- length(S)
X <- log(S[-1]/S[-n])
plot(X, type = "l", main = "AMT Stock log returns")
abline(v = cpoint(X)$tau0, col = "red")

## log ret CP onwards
getSymbols("AMT", from="2020-02-07", to="2021-07-21")
S <- AMT$AMT.Close
X <- diff(log(S))
plot(X)

## Density plot and QQplot from the change-point onwards
library(quantmod)
getSymbols("AMT", from="2020-02-07", to="2021-07-21")
S <- AMT$AMT.Close
lineChart(S, name="AMT",theme="white", up.col= "Black")

## Construction of the log-returns
X <- diff(log(S))
plot(X)
head(X)
X <- na.omit(X)
head(X)

## Density plot
plot(density(X), lwd=2, main="AMT Stock Density Plot")
f <- function(u) dnorm(u, mean=mean(X), sd=sd(X))
curve( f, -0.1, 0.1, add=TRUE, col="red",lwd=2)

## QQplot
qqnorm(X, main="AMT Stock Q-Q plot")
qqline(X, col="red",lwd=2)

----------

## Estimation of the B&S parameters
  #historical volatility
require(quantmod)
getSymbols("AMT", from="2020-02-07", to="2021-07-21")
attr(AMT, "src")
S <- AMT$AMT.Close
plot(S)
X <- na.omit(diff(log(S)))
plot(X)
Delta <- 1/252
alpha.hat <- mean(X, na.rm=TRUE)/Delta
sigma.hat <- sqrt(var(X, na.rm=TRUE)/Delta)
mu.hat <- alpha.hat + 0.5*sigma.hat^2
sigma.hat
mu.hat

##implied volatility
library(fOptions)
S0 <- as.numeric(tail(AMT$AMT.Close, n=1))
K <- 250
T <- 22 * Delta
r <- 0.0149
p<- 35
sigma.imp <- GBSVolatility(p, "c", S = S0, X = K, Time = T, r = r, b = r)
sigma.imp

#max likelihood estimation
set.seed (123)
library ("stats4")
x <- na.omit(diff(log(S)))
log.lik <- function (mu = 0.1922518, sigma =  0.3861056) -sum ( dnorm (x, mean = mu ,
                                                         sd = sigma , log = TRUE ))
fit <- mle (log.lik)
fit
summary(fit)

------------------------------------------------

#### Simulation part ####
library(tseries)
library(quantmod)
library(sde)

# B&S estimates from 2020-02-07 till 2021-07-21
SS <- get.hist.quote("AMT", start = "2020-02-07", end = "2021-07-21")
SS <- SS$Close
XX <- diff(log(SS))
str (SS)
Delta<- 1/252
alpha.hat2 <- mean(XX,na.rm=TRUE)/Delta
sigma.hat2 <- sqrt(var(XX,na.rm=TRUE)/Delta)
mu.hat2<- alpha.hat2 + 0.5*sigma.hat2^2
mu.hat2
sigma.hat2


## Data
getSymbols("AMT", from="2020-02-07", to="2021-07-21")   # starting price
S0 <- AMT[,"AMT.Close"]       # closing prices for the simulated period

## GBM Simulation
set.seed(123)
str(S0)                 # to work out correctly the length of the period
n <- 252                # set the length of the period
mu <- mu.hat2
sigma <- sigma.hat2
Ssim <- numeric(n)
Ssim[1] <- S0
for(i in 2:n)
  Ssim[i] <- Ssim[i-1] + mu*Ssim[i-1]*Delta + sigma*Ssim[i-1]*sqrt(Delta)*rnorm(1)
t <- (0:(n-1))*Delta

S0<- 283
#####B&S
#Call formula
call.price <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1, 
                       K = 250.00) {
  d2 <- (log(x/K) + (r - 0.5 * sigma^2) * (T - t))/(sigma * 
                                                      sqrt(T - t))
  d1 <- d2 + sigma * sqrt(T - t)
  x * pnorm(d1) - K * exp(-r * (T - t)) * pnorm(d2)
}

# PUT formula

put.price <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1, 
                      K = 250.00) {
  d2 <- (log(x/K) + (r - 0.5 * sigma^2) * (T - t))/(sigma * 
                                                      sqrt(T - t))
  d1 <- d2 + sigma * sqrt(T - t)
  K * exp(-r * (T - t)) * pnorm(-d2) - x * pnorm(-d1)
}

C <- call.price(x = S0, t = 0, T = T, r = r, K = K, sigma = sigma)
C

P <- put.price(x = S0, t = 0, T = T, r = r, K = K, sigma = sigma)
P


#parity
35.41098-283+250.00*exp(1)^-(0.01490*22/252)

## Rescalability
a<-4
GBSOption(TypeFlag = "c", S=a*S0, X=a*K, Time=T, r=r, b=r, sigma=sigma)@price
a*C

## Monte Carlo Pricing (includes antithetic technique)
MCPrice <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1,
                    M = 1000, f) {
  h <- function(m) {
    u <- rnorm(m/2)
    tmp <- c(x * exp((r - 0.5 * sigma^2) * (T - t) + sigma *
                       sqrt(T - t) * u), x * exp((r - 0.5 * sigma^2) * (T -
                                                                          t) + sigma * sqrt(T - t) * (-u)))
    mean(sapply(tmp, function(xx) f(xx)))
  }
  p <- h(M)
  p * exp(-r * (T - t))
}

# AMT Informationa & MC pricing (put)
S0 <- 283
K <- 250
r <- 0.0149
T <- 22/252
sigma <- 0.3861056
GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)@price
f <- function(x) max(0, K-x)
set.seed(123)
M <- 1000
MCPrice(x = S0, t = 0, T = T, r = r, sigma, M = M, f = f)
set.seed(123)
M <- 50000
MCPrice(x = S0, t = 0, T = T, r = r, sigma, M = M, f = f)
set.seed(123)
M <- 1e+06
MCPrice(x = S0, t = 0, T = T, r = r, sigma, M = M, f = f)


# Speed of convergence
set.seed(123)
m <- c(10, 50, 100, 150, 200, 250, 500, 1000)
p1 <- NULL
err <- NULL
nM <- length(m)
repl <- 100
mat <- matrix(, repl, nM)
for (k in 1:nM) {
  tmp <- numeric(repl)
  for (i in 1:repl) tmp[i] <- MCPrice(x = S0, t = 0, T = T,
                                      r = r, sigma, M = m[k], f = f)
  mat[, k] <- tmp
  p1 <- c(p1, mean(tmp))
  err <- c(err, sd(tmp))
}
colnames(mat) <- m
p0 <- GBSOption(TypeFlag = "p", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)@price
minP <- min(p1 - err)
maxP <- max(p1 + err)
plot(m, p1, type = "n", ylim = c(minP, maxP), axes = F, ylab = "MC price", xlab = "MC replications")
lines(m, p1 + err, col = "blue")
lines(m, p1 - err, col = "blue")
axis(2, p0, "B&S price")
axis(1, m)
boxplot(mat, add = TRUE, at = m, boxwex = 15, col = "orange", axes = F)
points(m, p1, col = "blue", lwd = 3, lty = 3)
abline(h = p0, lty = 2, col = "red", lwd = 3)

#FFT pricing CALL
FFTcall.price <- function(phi, S0, K, r, T, alpha = 1, N = 2^12, eta = 0.25) {
  m <- r - log(phi(-(0+1i)))
  phi.tilde <- function(u) (phi(u) * exp((0+1i) * u * m))^T
  psi <- function(v) exp(-r * T) * phi.tilde((v - (alpha + 
                                                     1) * (0+1i)))/(alpha^2 + alpha - v^2 + (0+1i) * (2 * 
                                                                                                        alpha + 1) * v)
  lambda <- (2 * pi)/(N * eta)
  b <- 1/2 * N * lambda
  ku <- -b + lambda * (0:(N - 1))
  v <- eta * (0:(N - 1))
  tmp <- exp((0+1i) * b * v) * psi(v) * eta * (3 + (-1)^(1:N) - 
                                                 ((1:N) - 1 == 0))/3
  ft <- fft(tmp)
  res <- exp(-alpha * ku) * ft/pi
  inter <- spline(ku, Re(res), xout = log(K/S0))
  return(inter$y * S0)
}


phiBS <- function(u) exp((0+1i) * u * (mu - 0.5 * sigma^2) - 0.5 * sigma^2 * u^2)

S0 <- 283
K <- 250
r <- 0.0149
T <- 22/252
sigma <- 0.3861056
mu <- 0.1922518 
require(fOptions)
GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r,   sigma = sigma)@price
FFTcall.price(phiBS, S0 = S0, K = K, r = r, T = T)
-----------------------------------------------------------------------------------
 

### American options: Explicit finite-difference method


AmericanPutExp <- function(Smin=0, Smax,  T=22/252, N=10, M=10, K, r=0.0149, sigma=0.3861056){
  Dt = T/N 
  DS = (Smax-Smin)/M
  t <- seq(0, T, by =Dt) 
  S <- seq(Smin, Smax, by=DS)
  A <- function(j) (-0.5*r*j*Dt + 0.5*sigma^2*j^2*Dt)/(1+r*Dt) 
  B <- function(j) (1-sigma^2*j^2*Dt)/(1+r*Dt) 
  C <- function(j) (0.5*r*j*Dt + 0.5*sigma^2*j^2*Dt)/(1+r*Dt)
  P <- matrix(, M+1, N+1)
  colnames(P) <- round(t,2)
  rownames(P) <- round(rev(S),2)
  P[M+1, ] <- K   # C(,j=0) = K
  P[1,] <- 0   # C(,j=M) = 0
  P[,N+1] <- sapply(rev(S), function(x) max(K-x,0))
  optTime <- matrix(FALSE, M+1, N+1)
  optTime[M+1,] <- TRUE
  optTime[which(P[,N+1]>0),N+1] <- TRUE
  
  for(i in (N-1):0){
    for(j in 1:(M-1)){
      J <- M+1-j
      I <- i+1
      P[J,I] <- A(j)*P[J+1,I+1] + B(j)*P[J,I+1] + C(j)*P[J-1,I+1]
      if(P[J,I] < P[J,N+1])
        optTime[J,I] <- TRUE
    }
  }
  colnames(optTime) <- colnames(P)
  rownames(optTime) <- rownames(P)
  ans <- list(P=P, t=t, S=S, optTime=optTime,N=N,M=M)
  class(ans) <- "AmericanPut"
  return(invisible(ans))
}


### plot method for American put

plot.AmericanPut <- function( obj ){
  plot(range(obj$t),range(obj$S),type="n",axes=F,xlab="t", ylab="S")
  axis(1,obj$t,obj$t)
  axis(2,obj$S,obj$S)
  abline(v = obj$t, h = obj$S, col = "darkgray", lty = "dotted")
  for(i in 0:obj$N){
    for(j in 0:obj$M){
      J <- obj$M+1-j
      I <- i+1
      cl <- "grey"; 
      if(obj$optTime[J,I])
        cl <- "black"
      text(obj$t[i+1],obj$S[j+1], round(obj$P[J,I],2),cex=0.75, col=cl)
    }
  }
  DS <- mean(obj$S[1:2])
  y <- as.numeric(apply(obj$optTime,2, function(x) which(x)[1]))
  lines(obj$t, obj$S[obj$M+2-y]+DS, lty=2)
}


put <- AmericanPutExp(Smax =  400, sigma = 0.3861056 , K = 250)
round(put$P,2)
par(mar=c(3,3,1,1))
plot(put)


S0 <- 283
myval <- round(put$P[which(rownames(put$P)==S0),1],2)
myval

### instability

put.bad <- AmericanPutExp(Smax = 400, sigma = 0.3861056, K = 250, M=15)
round(put.bad$P,2)



### American options: finite-difference implicit method

AmericanPutImp <- function( Smin=0, Smax,  T=22/252, N=10, M=10, K, r=0.0149, sigma=0.3861056){
  Dt = T/N 
  DS = (Smax-Smin)/M
  t <- seq(0, T, by =Dt) 
  S <- seq(Smin, Smax, by=DS)

  A <- function(j) 0.5*r*j*Dt - 0.5*sigma^2*j^2*Dt 
  B <- function(j) 1+sigma^2*j^2*Dt+r*Dt
  C <- function(j) -0.5*r*j*Dt - 0.5*sigma^2*j^2*Dt
  
  a <- sapply(0:M, A)
  b <- sapply(0:M, B)
  c <- sapply(0:M, C)
  
  P <- matrix(, M+1, N+1)
  colnames(P) <- round(t,2)
  rownames(P) <- round(rev(S),2)
  
  P[M+1, ] <- K   # C(,j=0) = K
  P[1,] <- 0   # C(,j=M) = 0
  P[,N+1] <- sapply(rev(S), function(x) max(K-x,0))
  
  AA <- matrix(0, M-1, M-1)
  for(j in 1:(M-1)){
    if(j>1) AA[j,j-1] <- A(j)
    if(j<M) AA[j,j] <- B(j)
    if(j<M-1) AA[j,j+1] <- C(j)
  }
  
  optTime <- matrix(FALSE, M+1, N+1)
  for(i in (N-1):0){
    I <- i+1
    bb <- P[M:2,I+1]
    bb[1] <- bb[1]-A(1)*P[M+1-0,I+1]
    bb[M-1] <- bb[M-1]-C(M-1)*P[M+1-M,I+1] 
    P[M:2,I] <- solve(AA,bb)
    idx <- which(P[,I] < P[,N+1])
    P[idx,I] <- P[idx,N+1] 
    optTime[idx, I] <- TRUE
  }
  optTime[M+1,] <- TRUE 
  optTime[which(P[,N+1]>0),N+1] <- TRUE
  colnames(optTime) <- colnames(P)
  rownames(optTime) <- rownames(P)
  ans <- list(P=P, t=t, S=S, optTime=optTime,N=N,M=M)
  class(ans) <- "AmericanPut"
  return(invisible(ans))
}



put <- AmericanPutImp(Smax = 400, sigma = 0.3861056, K = 250)
round(put$P,2)


par(mar=c(3,3,1,1))
plot(put)


# Barone-Adesi Wiley approx
install.packages("fOptions")
library(fOptions)
T <- 22/252
sigma= 0.3861056
r=0.0149
S0 <- 283
K <- 250
BAWAmericanApproxOption("p", S=S0, X=K, Time=T, r=r, b=r, sigma=sigma)@price

put <- AmericanPutImp(0,400,T=T, K=K,r=r,sigma=sigma,M=500,N=00)
put$P["283",1]
?fOptions
#include it in the code!


# Broadie and Glasserman Monte Carlo method

simTree <- function(b,d, S0, sigma, T, r){
  tot <- sum(b^(1:(d-1)))
  S <- numeric(tot+1) 
  S[1] <- S0
  dt <- T/d
  for(i in 0:(tot - b^(d-1))){
    for(j in 1:b){
      S[i*b+j +1] <- S[i+1]*exp((r - 0.5*sigma^2)*dt + sigma*sqrt(dt)*rnorm(1))
    }
  }
  S
}



upperBG <- function(S, b, d, f){
  tot <- sum(b^(1:(d-1)))
  start <- tot - b^(d-1) +1
  end <- tot +1
  P <- S
  P[start:end] <- f(S[start:end])
  tot1 <- sum(b^(1:(d-2)))
  for(i in tot1:0){
    m <- mean(P[i*b+1:b+1])
    v <- f(S[i+1])
    P[i+1] <- max(v,m)
  }
  P
}

lowerBG <- function(S, b, d, f){
  tot <- sum(b^(1:(d-1)))
  start <- tot - b^(d-1) +1
  end <- tot +1
  p <- S 
  p[start:end] <- f(S[start:end])
  tot1 <- sum(b^(1:(d-2)))
  
  m <- numeric(b)
  for(i in tot1:0){
    v <- f(S[i+1])
    for(j in 1:b){
      m[j] <- mean(p[i*b+(1:b)[-j]+1])
      m[j] <- ifelse( v>m[j], v, p[i*b+(1:b)[j]+1])
    }
    p[i+1] <- mean(m)
  }
  p
}



lowerBG(S, b,d,f)
upperBG(S, b,d,f)


# example of use
set.seed(123)
b <- 5
d <- 5
K <- 250
f <- function(x) sapply(x, function(x) max(K-x,0)) 
T <- 22/252
r <- 0.0149
sigma <- 0.3861056
S0 <- 283

low <- 0
upp <- 0
M <- 1000
for(i in 1:M){
  S <- simTree(b,d, S0, sigma, T, r)
  low <- low + lowerBG(S, b,d,f)[1]
  upp <- upp + upperBG(S, b,d,f)[1]
}
low/M
upp/M




### least squares method

LSM <- function(n, d, S0, K, sigma, r, T){
  s0 <- S0/K
  dt <- T/d
  z <- rnorm(n)
  s.t <- s0*exp((r-1/2*sigma^2)*T+sigma*z*(T^0.5))
  s.t[(n+1):(2*n)] <- s0*exp((r-1/2*sigma^2)*T-sigma*z*(T^0.5))
  CC <- pmax(1-s.t, 0)
  payoffeu <- exp(-r*T)*(CC[1:n]+CC[(n+1):(2*n)])/2*K
  euprice <- mean(payoffeu)
  
  for(k in (d-1):1){
    z <- rnorm(n)
    mean <- (log(s0) + k*log(s.t[1:n]))/(k+1)
    vol <- (k*dt/(k+1))^0.5*z
    s.t.1 <- exp(mean+sigma*vol)
    mean <- (log(s0) + k*log( s.t[(n+1):(2*n)] )) / ( k + 1 )
    s.t.1[(n+1):(2*n)] <- exp(mean-sigma*vol)
    CE <- pmax(1-s.t.1,0)
    idx<-(1:(2*n))[CE>0]
    discountedCC<- CC[idx]*exp(-r*dt)
    basis1 <- exp(-s.t.1[idx]/2)
    basis2 <- basis1*(1-s.t.1[idx])
    basis3 <- basis1*(1-2*s.t.1[idx]+(s.t.1[idx]^2)/2)
    
    p <- lm(discountedCC ~ basis1+basis2+basis3)$coefficients
    estimatedCC <- p[1]+p[2]*basis1+p[3]*basis2+p[4]*basis3
    EF <- rep(0, 2*n)
    EF[idx] <- (CE[idx]>estimatedCC)
    CC <- (EF == 0)*CC*exp(-r*dt)+(EF == 1)*CE
    s.t <- s.t.1
  }
  
  payoff <- exp(-r*dt)*(CC[1:n]+CC[(n+1):(2*n)])/2
  usprice <- mean(payoff*K)
  error <- 1.96*sd(payoff*K)/sqrt(n)
  earlyex <- usprice-euprice
  data.frame(usprice, error, euprice)
}


### comparison

S0 <- 283
K <- 250
T <- 22/252
r <- 0.0149
sigma <- 0.3861056
LSM(10000, 3, K, S0, sigma, r, T)
require(fOptions)
BSAmericanApproxOption("p", S0,K,T,r,r,sigma)@price
BAWAmericanApproxOption("p",S0, K, T,r,r, sigma)@price
put <- AmericanPutImp(0,100,T=T, K=K,r=r,sigma=sigma,M=100,N=100)
put$P["283",1]


# Double Monte Carlo experiment
set.seed(123)
b <- 3
d <- 3
K <- 230
f <- function(x) sapply(x, function(x) max(K-x,0))
T <- 22/252
r <- 0.0149
sigma <- 0.3861056
S0 <- 283

low <- 0
upp <- 0
M <- 1000
for(i in 1:M){
  S <- simTree(b,d, S0, sigma, T, r)
  low <- low + lowerBG(S, b,d,f)[1]
  upp <- upp + upperBG(S, b,d,f)[1]
}
low/M
upp/M


# increasing the number of branches at each node
set.seed(123)
b <- 5
d <- 3
K <- 250
f <- function(x) sapply(x, function(x) max(K-x,0))
T <-22/252
r <- 0.0149
sigma <- 0.3861056
S0 <- 283

low <- 0
upp <- 0
M <- 1000
for(i in 1:M){
  S <- simTree(b,d, S0, sigma, T, r)
  low <- low + lowerBG(S, b,d,f)[1]
  upp <- upp + upperBG(S, b,d,f)[1]
}
low/M
upp/M

#too high value



############################################################################################
#Lévy PROCESSES

#Packages
install.packages("quantmod")
install.packages("moments")
install.packages("VarianceGamma")
install.packages("viridis")
install.packages("seqHMM")
install.packages('fBasics')

#Obtaining the dataset -> S&P500 
library(quantmod)   
library(moments)
library(seqHMM)
require(fBasics)
getSymbols("AMT", from="2020-02-07", to="2021-07-21")
attr(AMT, "src")
S <- AMT [,"AMT.Close"]
l_ret <- diff(log(S))
l_ret <- na.omit(l_ret)

#VARIANCE-GAMMA MODEL

#Parameters fitting
library(VarianceGamma)
library(viridis)

vgFit(l_ret, method = "nlm")
vgFit(l_ret, method = "BFGS")
vgFit(l_ret)#estimated VG parameters on the sample nelder mead
str(vgFit(l_ret))
vg_param <- vgFit(l_ret)$param


c <- as.numeric(vg_param[1])
sigma <- as.numeric(vg_param[2])
theta <- as.numeric(vg_param[3])
nu <- as.numeric(vg_param[4])

T <- 22/252 #option maturity
N <- 100    #number of steps for each path
r <- 0.0149    #arbitrary risk-free rate
nsim <- 1000  #number of simulated path


#Variance Gamma Function
VG=function(sigma, nu, theta, Tf, N, r) { 
  a=1/nu 
  b=1/nu 
  h=T/N
  t=(0:N)*T/N 
  X=rep(0, N+1) 
  I=rep(0,N) 
  X[1]=0 
  for(i in 1:N) { 
    I[i]=rgamma(1,a*h,b) 
    X[i+1]=X[i] + theta*I[i]+sigma*sqrt(I[i])*rnorm(1)
  }
  return((X)) }


set.seed(123)
VG_paths<-matrix(nrow = nsim, ncol=N+1)  #fill the matrix with random paths that follow
for(i in 1:nsim){                        #the function VG just created
  VG_paths[i,]<-VG(sigma, nu, theta, T, N, r)
}

VG_paths



#plot the Monte Carlo Simulation
library("viridis")

colori=viridis(nsim)
plot(VG_paths[1,], col=0, type="l", ylim = c(min(VG_paths),max(VG_paths)), 
     main = "Monte Carlo Simulation for VG returns", sub = "100 steps, 10 paths", 
     xlab = "Time", ylab = "VG returns")
for(i in 2:nsim){
  lines(VG_paths[i,], col=colori[i], lwd = 2);
}


##TESTS (both graphical and not) OF DISTRIBUTIONAL ASSUMPTIONS

#QQplot

l_ret.s <- sort(as.numeric(l_ret)) 
p <- ppoints(length(l_ret.s)) #plotting position??? 
VG.q <- qvg(p, vgC=c, sigma=sigma, theta=theta, nu=nu) #compute the quantile 
plot(VG.q, l_ret.s, main = "Variance-Gamma Q-Q Plot", 
     xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", abline(0, 1, col = 2, lwd=3)) 



#good result, linear

#Density comparison
#kernel density and VG overlayed (Gaussian kernel, Silverman's rule of thumb)
plot(density(l_ret[-1,]), type = "l", lwd = 2, lty = 3, col = "red",
     xlim= c(-0.03,0.03), ylim=c(0,120), main ="", xlab ="", ylab = "")
legend ("topright", inset = .02, c("Kernel", "VG"),
        col=c("red","black"), lwd=c(2,1), lty=c(3,1), cex = 0.8, bty = "n")
points(seq(min(l_ret[-1,]), max(l_ret[-1,]), length.out=500), 
       dvg(seq(min(l_ret[-1,]), max(l_ret[-1,]), length.out=500),
           mean(l_ret[-1,]), sd(l_ret[-1,])), type="l", col="black")
#better fitting than the log-normal case

#gridplot <- seq(min(l_ret[-1,]), max(l_ret[-1,]), length.out=length(l_ret[-1,]))
#plot(density(l_ret[-1,]), type = "l", lwd = 2, lty = 3, col = "coral2",
#     xlim= c(-0.03,0.03), ylim=c(0,120), main ="", xlab ="", ylab = "")
#legend ("topright", inset = .02, c("Kernel", "VG"),
#        col=c("coral2","seagreen3"), lwd=c(2,1), lty=c(3,1), cex = 0.8, bty = "n")
#points(gridplot, dvg(gridplot, param = c(c, sigma, theta, nu)),
#       type="l", col="seagreen3")
##better code (we include the VG parameters), worse result



#Tests

#H0 = The data is consistent with a specified reference distribution.
#H1 = The data is NOT consistent with a specified reference distribution

#Chi^2 test
chisq.test(l_ret.s, VG.q)
#high p-value (0.24) -> we can't reject the null hypotesis


#K-S test
ks.test(as.numeric(l_ret), rvg(length(as.numeric(l_ret)), 
                               param = c(c, sigma, theta, nu)))
#high p-value (0.10) -> we can't reject the null hypotesis



#summary statistics
final_retVG<-VG_paths[,N]
basicStats(final_retVG)
hist(final_retVG) #not much disclosing when nsim is small 



#####################################################################################



#FROM VARIANCE-GAMMA RETURNS TO STOCK PRICES

S0 <- as.numeric(tail(AMT$AMT.Close, n=1)) 
S0


r <- 0.0149 

#function for stock price with VG returns
VGexp=function(sigma, nu, theta, T, N, r, S0) { 
  a=1/nu 
  b=1/nu 
  h=T/N
  t=(0:N)*T/N 
  X=rep(0, N+1) 
  I=rep(0,N) 
  X[1]=S0 
  for(i in 1:N) { 
    I[i]=rgamma(1,a*h,b) 
    X[i+1]=X[i]*exp(r*t+theta*I[i]+sigma*sqrt(I[i])*rnorm(1))
  }
  return(X)}


set.seed(123)
VGexp_paths<-matrix(nrow = nsim, ncol=N+1)
for(i in 1:nsim){
  VGexp_paths[i,]<-VGexp(sigma, nu, theta, T, N, r, S0)
}

VGexp_paths


#plot MCS
plot(VGexp_paths[1,], col=0, type="l", ylim = c(min(VGexp_paths),max(VGexp_paths)), 
     main = "MC Simlation for VG stock prices", sub = "100 steps, 10 paths", 
     xlab = "Time", ylab = "S&P 500")
for(i in 2:nsim){
  lines(VGexp_paths[i,], col=colori[i], lwd = 2);
  
}


#statistics on final prices
final_pricesVG<-VGexp_paths[,N]

#risk neutral transform
rn_final_pricesVG<-S0*(final_pricesVG)*(exp(r*T)/(mean(final_pricesVG)))
rn_final_pricesVG

basicStats(rn_final_pricesVG)




##CALL OPTION PRICING 
K <- 250
payoff_VG <- pmax(rn_final_pricesVG - K, 0)
optprice_VG <- mean(payoff_VG)*exp(-r*T)
optprice_VG


#############################################################
#######Pricing with FFT under the VG process
# VG process
#Nelder-Mead
theta <- 0.0006035
nu <- 1.3570967
r <- 0.0149 
sigma <- 0.0230158
T <- 22/252
K <- 250
S <- 283
alpha <- 1.65
phiVG <- function(u) {
  omega <- (1/nu) * (log(1 - theta * nu - sigma^2 * nu/2))
  tmp <- 1 - (0+1i) * theta * nu * u + 0.5 * sigma^2 * u^2 * nu
  tmp <- tmp^(-1/nu)
  exp((0+1i) * u * log(S0) + u * (r + omega) * (0+1i)) * tmp
}

FFTcall.price <- function(phi, S0, K, r, T, alpha = 1, N = 2^12, eta = 0.25) {
  m <- r - log(phi(-(0+1i)))
  phi.tilde <- function(u) (phi(u) * exp((0+1i) * u * m))^T
  psi <- function(v) exp(-r * T) * phi.tilde((v - (alpha + 
                                                     1) * (0+1i)))/(alpha^2 + alpha - v^2 + (0+1i) * (2 * 
                                                                                                        alpha + 1) * v)
  lambda <- (2 * pi)/(N * eta)
  b <- 1/2 * N * lambda
  ku <- -b + lambda * (0:(N - 1))
  v <- eta * (0:(N - 1))
  tmp <- exp((0+1i) * b * v) * psi(v) * eta * (3 + (-1)^(1:N) - 
                                                 ((1:N) - 1 == 0))/3
  ft <- fft(tmp)
  res <- exp(-alpha * ku) * ft/pi
  inter <- spline(ku, Re(res), xout = log(K/S0))
  return(inter$y * S0)
}

FFTcall.price(phiVG, S0 = S0, K = K, r = r, T = T)

#######################################################################################


#MEIXNER MODEL

install.packages("fBasics")
install.packages("Runuran")

require('fBasics')
require('Runuran')

#Moments: mean, variance, skewness, kurtosis
x   <-mean(l_ret, na.rm = TRUE)
y   <-sd(l_ret, na.rm = TRUE)
z <-as.numeric(skewness(l_ret, na.rm = TRUE))
w <-as.numeric(kurtosis(l_ret, na.rm = TRUE))


#Mom: estimates parameters m, a, b, d as functions of the moments
m <- x-((z*sqrt(y))/(w-(z^2)-3)) 
a <- sqrt(y*(2*w-3*(z^2)-6)) 
b <- 2*atan(-sqrt((z^2)/(2*w-3*(z^2)-6))) 
d <- 1/(w-(z^2)-3) 

m
a
b
d

#risk neutral transformation 1
#Esscher transform: Meixner(a, a*theta + b, d, m) distribution
theta <- -1/a * (b + 2 * atan((-cos(a/2)+ exp((m-r)/2*d))/sin(a/2)))
b <- a*theta+b

#risk neutral transformation 2
#final prices
rn_final_pricesMX<-S0*(final_pricesMX)*(exp(r*Tf)/(mean(final_pricesMX)))
rn_final_pricesMX


#risk neutral transformation 3
# modifying m
m <- r -2 *d*log(cos(b/2)/cos((a+b)/2))


library(fBasics) 
library(Runuran)

#Meixner function
MX=function(a, b, d, m, N) {
  distr <- udmeixner(a, b, d, m) #meixner distribution
  gen <- pinvd.new(distr) #Polynomial interpolation of INVerse CDF
  rdmMXgen <- ur(gen,N) #randomly draws N objects from gen (from a Meixner distr)
  h=T/N
  X=rep(0,N+1)
  for (i in 1:N){
    X[i+1]=X[1]+rdmMXgen[i]
  }
  return(X)
}

MX(a, b, d, m, N)


set.seed(123)
MX_paths<-matrix(nrow = nsim, ncol=N+1) #fill the matrix with random paths that follow
for(i in 1:nsim){                        #the function MX just created
  MX_paths[i,]<-MX(a,b,d,m,N)
}


#plot the Monte Carlo Simulation
plot(MX_paths[1,], col=0, type="l", ylim = c(min(MX_paths),max(MX_paths)), 
     main = "Monte Carlo Simulation for Meixner returns", sub = "100 steps, 10 paths", 
     xlab = "Time", ylab = "MXNR returns")
for(i in 2:nsim){
  lines(MX_paths[i,], col=colori[i], lwd = 2);
}



#QQplot
MX.q <- uq(pinvd.new(udmeixner(a, b, d, m)), p) #compute the quantile

plot(MX.q, l_ret.s, main = "Meixner Q-Q Plot",
     xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(MX.q, 
       lwd = 2,      # spessore
       col = "red"   # colore
)
#good result, linear

#summary statistics
final_retMX<-MX_paths[,N]
basicStats(final_retMX)
hist(final_retMX) #not much disclosing when nsim is small 


#####################################################################################



#FROM MEIXNER RETURNS TO STOCK PRICES

#function for stock price with Meixner returns
MXexp=function(a, b, d, m, N, T, r, S0) {
  distr <- udmeixner(a, b, d, m) #meiner distribution
  gen <- pinvd.new(distr) #Polynomial interpolation of INVerse CDF
  generazioni <- ur(gen,N) #randomly draws N objects from gen (from a Meixner distr)
  h=T/N
  t=(0:N)*T/N
  X=rep(0,N+1)
  X[1]=S0
  for (i in 1:N){
    X[i+1]=X[1]*exp(r*t+generazioni[i])
  }
  return(X)
}

MXexp(a, b, d, m, N, T, r, S0)


set.seed(123)
MXexp_paths<-matrix(nrow = nsim, ncol=N+1)
for(i in 1:nsim){
  MXexp_paths[i,]<-MXexp(a,b,d,m,100,T,r,S0) #vengono tutte le linee uguali perché MX non varia!
}

MXexp_paths

#statistics on final prices
final_pricesMX<-MXexp_paths[,N]


########payoff if you use the Esscher transform 

payoff_MX <- pmax(final_pricesMX - K, 0)

optprice_MX <- mean(payoff_MX)*exp(-r*T)

optprice_MX

########payoff rn transformation of final prices

rn_payoff_MX <- pmax(rn_final_pricesMX - K, 0)

optprice_MX <- mean(rn_payoff_MX)*exp(-r*T)

optprice_MX

----------------------------------------------------------------------------
#CHANGE POINT DETECTION
#CCI
library(tseries)
CCI <- get.hist.quote("CCI", start = "2016-01-01", end = "2021-07-21")
chartSeries(CCI, up.col= "Black",TA=c(addVo(), addBBands()), theme="white")
CCI <- CCI$Close
require(sde)
cpoint(CCI)

#prices CP line
addVLine = function(dtlist) plot(addTA(xts(rep(TRUE,NROW(dtlist)),dtlist),on=1,col="red"))
addVLine(cpoint(S)$tau0)

#logret CP line CCI
CCI <- as.numeric(CCI)
n <- length(CCI)
CCI <- log(CCI[-1]/CCI[-n])
plot(CCI, type = "l", main = "CCI Stock log returns")
abline(v = cpoint(CCI)$tau0, col = "red")

-----------------------------------------------------------------------------

library(Quandl)
Quandl.api_key('-t1iySvb5iK76Eg6M_aC')
require(quantmod)

#AMT
getSymbols("AMT", from="2020-02-07", to="2021-07-21")
AMT<- AMT$AMT.Close
ts(AMT)
plot.ts(AMT, main="timeseries")
X <- na.omit(diff(log(AMT)))
plot(X, main="AMT 1yr return")
plot(density(X), type="h", lwd=2, main="AMT density plot")
f <- function(u) dnorm(u, mean=mean(X), sd=sd(X)) 
curve( f, -0.1, 0.1, add=TRUE, col="red", lwd=2, type = "l")
qqnorm(X, main="AMT Q-Q plot")
qqline(X, col="green",lwd=2)

#CCI
getSymbols("CCI", from="2020-02-07", to="2021-07-21")
CCI<- CCI$CCI.Close
ts(CCI)
plot.ts(CCI, main="timeseries")
Y <- na.omit(diff(log(CCI)))
plot(Y, main="CCI 1yr return")
plot(density(Y), type="h", lwd=2, main="CCI density plot")
f <- function(u) dnorm(u, mean=mean(Y), sd=sd(Y)) 
curve( f, -0.1, 0.1, add=TRUE, col="red", lwd=2, type = "l")
qqnorm(Y, main="CCI Q-Q plot")
qqline(Y, col="green",lwd=2)


#Volatilities
str(X)## to work out correctly the length of the period
Delta<- 1/252
alpha.hat2 <- mean(X,na.rm=TRUE)/Delta
sigma.hat.AMT <- sqrt(var(X,na.rm=TRUE)/Delta)
mu.hat2<- alpha.hat2 + 0.5*sigma.hat.AMT^2
mu.hat2
sigma.hat.AMT

str(Y)## to work out correctly the length of the period
Delta<- 1/252
alpha.hat3 <- mean(Y,na.rm=TRUE)/Delta
sigma.hat.CCI <- sqrt(var(Y,na.rm=TRUE)/Delta)
mu.hat3 <- alpha.hat3 + 0.5*sigma.hat.CCI^2
mu.hat3
sigma.hat.CCI

volatilities <- c(sigma.hat.AMT, sigma.hat.CCI)
volatilities


###################################COPULAS###################################
require(copula)
require(qrng)
require(quantmod)
X <- diff(log(AMT))
A <- c(X)
W <- diff(log(CCI))
B <- c(W)
M <- cbind(A, B)
corKendall(as.matrix(M))
#alternatively cor(as.matrix(M), method="kendall")


#### MULTIASSET ####
#Rho Lin Correlation
# R script for creating GBM paths of correlated daily asset prices. In this 
# example, we do this for automotive listed companies as: FCA, BMW and Daimler Group. 
####

GBM <- function(N, sigma, mu, S0, Wt = NULL) {
  # Creates a single asset path of daily prices using Geometric Brownian Motion. 
  # One year is 252 days since that is about how many trading days are in any
  # given year.
  #
  # Args:
  #   N: Number of days in the path.
  #   sigma: Volatility or standard deviation of daily continuously compounded 
  #          returns.
  #   mu: Drift or average daily continuously compounded returns. 
  #   S0: The initial price of the asset. 
  #   Wt: The cumulative Brownian motion of the model. This can be supplied or 
  #       left as NULL. In the case that it is NULL, a vector will be provided.
  #       If you include this argument, it must be a vector of length N of the 
  #       cumulative sum of a random variable to work properly. 
  #
  # Returns:
  #   A vector of length N containing the asset prices generated by the specified
  #   GBM. 
  if (is.null(Wt)) {
    Wt <- cumsum(rnorm(N, 0, 1))
  }
  t <- (1:N)/252
  p1 <- (mu - 0.5*(sigma*sigma)) * t
  p2 <- sigma * Wt
  St = S0 * exp(p1 + p2)
  return(St)
}

CorrelatedGBM <- function(N, S0, mu, sigma, cor.mat) {
  # Creates a matrix of correlated daily price paths using Geometric 
  # Brownian Motion. 
  #
  # Args: 
  #   N: Number of days in the path.
  #   mu: Drift or average daily continuously compounded returns.  
  #   sigma: Volatility or standard deviation of daily continuously compounded 
  #          returns. 
  #   S0: The initial price of the asset. 
  #   cor.mat: The correlation matrix of the daility contiuously compounded 
  #            returns. 
  #
  # Returns:
  #   A matrix of simulated daily price paths of length N having the same number
  #   of assets as in the mu and sigma vectors. Note that mu and sigma must have
  #   the same dimensions. 
  mu <- as.matrix(mu)
  sigma <- as.matrix(sigma)
  GBMs <- matrix(nrow = N, ncol = nrow(mu))
  Wt <- matrix(rnorm(N * nrow(mu), 0, 1), ncol = nrow(mu))
  Wt <- apply(Wt, 2, cumsum)
  chol.mat <- chol(cor.mat) # upper triangular cholesky decomposition
  Wt <- Wt %*% chol.mat   # key trick for creating correlated paths
  for (i in 1:nrow(mu)) {
    GBMs[,i] <- GBM(N, sigma[i], mu[i] , S0[i], Wt[, i])
  }
  return(GBMs)
}

GetPrices <- function(tickers, startDate='2020-02-07') {
  # Gets price information from Yahoo based on ticker information and the 
  # specified start date. 
  #
  # Args:
  #   tickers: A vector of the ticker symbols to download.
  #   startDate: The beginning date. Defaults to 1/2/1992.
  #
  # Returns:
  #   A zoo object with the price information. 
  prices <- get.hist.quote(instrument = tickers[1], start = startDate, 
                           quote = 'AdjClose')
  # download the rest of the prices
  for (tik in 2:length(tickers)) {
    tmp <- get.hist.quote(instrument = tickers[tik], 
                          start = start, quote = 'AdjClose')
    prices <- merge(prices, tmp)
  }
  return(prices)
}

require(zoo)
require(tseries)

set.seed(1994)
N <- 1 * 252  #252 trading days
t <- (1:N)/252
prices <- cbind.zoo(AMT, CCI)

# get the cc returns and vectors of average returns and volatiliy
returns.mat <- as.matrix(na.omit(diff(log(prices))))
mean.vec <- as.numeric(colMeans(returns.mat))
sigma.vec <- as.numeric(sqrt(apply(returns.mat, 2, var)))
prices.vec <- as.numeric(prices[nrow(prices)])
cor.mat <-as.matrix(cor(returns.mat))

# make correlated asset paths
paths <- CorrelatedGBM(N, prices.vec , mean.vec, sigma.vec, cor.mat)

# make a basic r plot 
colors <- c('black', 'darkred')
plot(t, paths[,1], type = 'l', ylim = c(0, max(paths)), xlab = "Year", 
     ylab = "Price", main = "Simulated Asset Prices", col = colors[1])
for (i in 2:ncol(paths)) {
  lines(t, paths[, i], col = colors[i])
}
legend(x = 2.9, y = 1725, c('AMT', 'CCI'), lty = c(1,1), col = colors, cex = 0.7)



##' @title Risk measures Value-at-Risk (VaR), Expected Shortfall (ES)
##'        and Capital Allocated on Aggregate of Losses
##' @param x (n, d)-matrix of losses
##' @param alpha confidence level
##' @return 5-vector of VaR, ES, Alloc(X.first), Alloc(X.mid), Alloc(X.last)
##' @author Mathieu Cambou, Christiane Lemieux, Marius Hofert
risk.measures <- function(x, alpha)
{
  if(!is.matrix(x)) x <- rbind(x)
  n <- nrow(x)
  d <- ncol(x)
  
  aloss  <- rowSums(x) # n-vector of aggregated losses
  VaR <- quantile(aloss, probs=alpha, names=FALSE) # VaR estimate
  l <- sum(xcd <- aloss > VaR)
  ES <- mean(aloss[xcd]) / (1-alpha) # ES estimate
  
  Alloc.first <- x[aloss>=VaR, 1] %*% rep(1/n, l) / (l/n) # capital allocated to X_1
  Alloc.mid   <- x[aloss>=VaR, floor(d/2)] %*% rep(1/n, l) / (l/n) # capital allocated to X_{floor(d/2)}
  Alloc.last  <- x[aloss>=VaR,d] %*% rep(1/n, l) / (l/n) # Capital Allocated to X_d
  
  ## return estimated risk measures
  c(VaR = VaR, ES = ES, Alloc.first = Alloc.first, Alloc.mid = Alloc.mid,
    Alloc.last = Alloc.last)
}



########################## 2 Case Study ########################################
### 2.1 Define parameters #####
n <- 1e5 # Monte Carlo sample size
d <- 2 # dimension

nu <- 3
volatilities
## Stochastic process parameters
sigma <- c(0.3861056, 0.3762671) #volatilities --> estimate my own sigma and put there 2 different sigma.
r <- 0.0149 #continuously compounded short rate --> rf rate used before
S0 <- c(238.75 ,  150.49 ) #initial stocks' levels
K <- 1.1 #option strike standardized (correspond to K*100)
N <- 1000 # option notional
T <- 22/252 # time horizon 


## Copulas
tau <- 0.642529 # Kendall's tau ----> corKendall(as.matrix(M)) 

#ALPHA.CONFIDENCE
alpha <- 0.99 # confidence level for VaR, ES



####FITTING####
#Clayton,Gumbel,Frank
loglikCopula(th.C,returns.mat, claytonCopula()) 
#returns the copula log-likelihood evaluated at the parameter (vector) param given the data u.
clay <- fitCopula(claytonCopula(dim = 2), returns.mat, method ="itau")
gumb <- fitCopula(gumbelCopula(dim = 2), returns.mat, method ="itau")
frank <- fitCopula(frankCopula(dim = 2), returns.mat, method ="itau")
summary(clay)
summary(gumb)
summary(frank)

#Frank###############################################################
family.F <- "Frank"
th.F <- iTau(getAcop(family.F), tau) # corresponding parameter
frank.cop <- onacopulaL(family.F, nacList=list(th.F, 1:d)) # Frank copula
### 2.2 Sampling 
set.seed(1994)
U.CDM  <- matrix(runif(n*d), ncol=d) # pseudo
U.C.CDM  <- cCopula(U.CDM,  cop=frank.cop, inverse=TRUE) # pseudo
### 2.3 Functional Calculation
erT <- exp(-r*T)
rm.C.CDM <- risk.measures(U.C.CDM, alpha)
### 2.4 Results
res <- array(dim=c(2,2,1), dimnames=list(type=c(paste0("VaR.", alpha), paste0("ES.", alpha)),
                                         copula=c("Frank", paste0("t", nu)),
                                         method=c("CDM")))
res[paste0("VaR.", alpha),,] <- matrix(c(rm.C.CDM[1]), ncol=1)
res[paste0("ES.", alpha),,]  <- matrix(c(rm.C.CDM[2]), ncol=1)
res

#Gumbel###############################################################
family.G <- "Gumbel"
th.G <- iTau(getAcop(family.G), tau) # corresponding parameter
gumbel.cop <- onacopulaL(family.G, nacList=list(th.G, 1:d)) # Gumbel copula
#gumbel.cop <- archmCopula(family.G, param=th.G, dim=d) # alternative cop.calc
### 2.2 Sampling 
set.seed(1994)
U.CDM  <- matrix(runif(n*d), ncol=d) # pseudo
U.C.CDM  <- cCopula(U.CDM,  cop=gumbel.cop, inverse=TRUE) # pseudo
### 2.3 Functional Calculation 
erT <- exp(-r*T)
rm.C.CDM <- risk.measures(U.C.CDM, alpha)
### 2.4 Results 
res <- array(dim=c(2,2,1), dimnames=list(type=c(paste0("VaR.", alpha), paste0("ES.", alpha)),
                                         copula=c("Gumbel", paste0("t", nu)),
                                         method=c("CDM")))
res[paste0("VaR.", alpha),,] <- matrix(c(rm.C.CDM[1]), ncol=1)
res[paste0("ES.", alpha),,]  <- matrix(c(rm.C.CDM[2]), ncol=1)
res

##Clayton###############################################################
family.C <- "Clayton"
th.C <- iTau(getAcop(family.C), tau) # corresponding parameter
clayton.cop <- onacopulaL(family.C, nacList=list(th.C, 1:d)) # Clayton copula
### 2.2 Sampling 
set.seed(1994)
U.CDM  <- matrix(runif(n*d), ncol=d) # pseudo
U.C.CDM  <- cCopula(U.CDM,  cop=clayton.cop, inverse=TRUE) # pseudo
### 2.3 Functional Calculation
erT <- exp(-r*T)
rm.C.CDM <- risk.measures(U.C.CDM, alpha)
### 2.4 Results 
res <- array(dim=c(2,2,1), dimnames=list(type=c(paste0("VaR.", alpha), paste0("ES.", alpha)),
                                         copula=c("Clayton", paste0("t", nu)),
                                         method=c("CDM")))
res[paste0("VaR.", alpha),,] <- matrix(c(rm.C.CDM[1]), ncol=1)
res[paste0("ES.", alpha),,]  <- matrix(c(rm.C.CDM[2]), ncol=1)
res


#TIME Copula calculation###############################################################
system.time(cCopula(U.CDM,  cop=frank.cop, inverse=TRUE))
system.time(cCopula(U.CDM,  cop=gumbel.cop, inverse=TRUE))
system.time(cCopula(U.CDM,  cop=clayton.cop, inverse=TRUE))

