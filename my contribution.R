require(quantmod)
getSymbols("AMT", from="2016-01-01", to="2021-07-21")
attr(AMT, "src")
S <- AMT [,"AMT.Close"]
X <- diff(log(S))
BAWValue
require("fOptions")
## Cox-Ross-Rubinstein Binomial Tree Option Model:
# call
CRRBinomialTreeOption(TypeFlag = "ce", S = 283, X = 250, 
                      Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)
CRRBinomialTreeOption(TypeFlag = "ca", S = 283, X = 250, 
                      Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)
# put 
CRRBinomialTreeOption(TypeFlag = "pe", S = 283, X = 250, 
                      Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)
CRRBinomialTreeOption(TypeFlag = "pa", S = 283, X = 250, 
                      Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)

##  Jarrow and Rudd (1983) Binomial Tree Option Model:
#call
JRBinomialTreeOption(TypeFlag = "ce", S = 283, X = 250, 
                      Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)
JRBinomialTreeOption(TypeFlag = "ca", S = 283, X = 250, 
                      Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)
#put
JRBinomialTreeOption(TypeFlag = "pe", S = 283, X = 250, 
                     Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)
JRBinomialTreeOption(TypeFlag = "pa", S = 283, X = 250, 
                     Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)

##  Tian Binomial Tree Option Model:
#call
TIANBinomialTreeOption(TypeFlag = "ce", S = 283, X = 250, 
                     Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)
TIANBinomialTreeOption(TypeFlag = "ca", S = 283, X = 250, 
                     Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)
#put
TIANBinomialTreeOption(TypeFlag = "pe", S = 283, X = 250, 
                     Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)
TIANBinomialTreeOption(TypeFlag = "pa", S = 283, X = 250, 
                     Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)
BAWValue
--------------------------------------------------------------------------------

## CRR - JR - TIAN Model Comparison plot, call US:  
par(mfrow = c(2, 1), cex = 0.9)
steps = 200
CRROptionValue =  JROptionValue = TIANOptionValue = 
  rep(NA, times = steps)
for (n in 3:steps) { 
  CRROptionValue[n] = CRRBinomialTreeOption(TypeFlag = "ca", S = 283, X = 250, 
                                            Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price
  JROptionValue[n] = JRBinomialTreeOption(TypeFlag = "ca", S = 283, X = 250, 
                                          Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price 
  TIANOptionValue[n] = TIANBinomialTreeOption(TypeFlag = "ca", S = 283, X = 250, 
                                              Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price 
}           
plot(CRROptionValue[3:steps], type = "l", col = "red", ylab = "Option Value")
lines(JROptionValue[3:steps], col = "green")
lines(TIANOptionValue[3:steps], col = "blue")
# Add Result from BAW Approximation:
BAWValue =  BAWAmericanApproxOption(TypeFlag = "c", S = 283, X = 250, 
                                    Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564)@price
abline(h = BAWValue, lty = 2)
title(main = "US call: Price Convergence (200 steps) ")
data.frame(CRROptionValue, JROptionValue, TIANOptionValue)
legend(150, 35.05, legend=c("CRR", "JR", "TIAN", "BAW"),
       col=c("red", "blue", "green", "black"), lty= c(1, 1, 1, 2), cex = 0.9)
## CRR - JR - TIAN Model Comparison plot, put US:  
steps = 200
CRROptionValue =  JROptionValue = TIANOptionValue = 
  rep(NA, times = steps)
for (n in 3:steps) { 
  CRROptionValue[n] = CRRBinomialTreeOption(TypeFlag = "pa", S = 283, X = 250, 
                                            Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price
  JROptionValue[n] = JRBinomialTreeOption(TypeFlag = "pa", S = 283, X = 250, 
                                          Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price 
  TIANOptionValue[n] = TIANBinomialTreeOption(TypeFlag = "pa", S = 283, X = 250, 
                                              Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price 
}           
plot(CRROptionValue[3:steps], type = "l", col = "red", ylab = "Option Value")
lines(JROptionValue[3:steps], col = "green")
lines(TIANOptionValue[3:steps], col = "blue")
# Add Result from BAW Approximation:
BAWValue =  BAWAmericanApproxOption(TypeFlag = "p", S = 283, X = 250, 
                                    Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564)@price
abline(h = BAWValue, lty = 2)
title(main = "US put: Price Convergence (200 steps)")
data.frame(CRROptionValue, JROptionValue, TIANOptionValue)
legend(150, 2.08, legend=c("CRR", "JR", "TIAN", "BAW"),
       col=c("red", "blue", "green", "black"), lty=c(1,1,1,2), cex = 0.9)

## CRR - JR - TIAN Model Comparison plot, call EU:  
par(mfrow = c(2, 1), cex = 0.9)
steps = 200
CRROptionValue =  JROptionValue = TIANOptionValue = 
  rep(NA, times = steps)
for (n in 3:steps) { 
  CRROptionValue[n] = CRRBinomialTreeOption(TypeFlag = "ce", S = 283, X = 250, 
                                            Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price
  JROptionValue[n] = JRBinomialTreeOption(TypeFlag = "ce", S = 283, X = 250, 
                                          Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price 
  TIANOptionValue[n] = TIANBinomialTreeOption(TypeFlag = "ce", S = 283, X = 250, 
                                              Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price 
}           
plot(CRROptionValue[3:steps], type = "l", col = "red", ylab = "Option Value")
lines(JROptionValue[3:steps], col = "green")
lines(TIANOptionValue[3:steps], col = "blue")
# Add Result from BAW Approximation:
abline(h = 35.41098, lty = 2) #B&S
abline(h = 35.0894, lty = 4, col = "orange")
title(main = "EU call: Price Convergence (200 steps)")
data.frame(CRROptionValue, JROptionValue, TIANOptionValue)
legend(150, 35.05, legend=c("CRR", "JR", "TIAN", "B&S", "Convergence Price"),
       col=c("red", "blue", "green", "black", "orange"), lty=c(1, 1, 1, 2, 4), cex = 0.9)

## CRR - JR - TIAN Model Comparison plot, put EU:  
steps = 200
CRROptionValue =  JROptionValue = TIANOptionValue = 
  rep(NA, times = steps)
for (n in 3:steps) { 
  CRROptionValue[n] = CRRBinomialTreeOption(TypeFlag = "pe", S = 283, X = 250, 
                                            Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price
  JROptionValue[n] = JRBinomialTreeOption(TypeFlag = "pe", S = 283, X = 250, 
                                          Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price 
  TIANOptionValue[n] = TIANBinomialTreeOption(TypeFlag = "pe", S = 283, X = 250, 
                                              Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = n)@price 
}           
plot(CRROptionValue[3:steps], type = "l", col = "red", ylab = "Option Value")
lines(JROptionValue[3:steps], col = "green")
lines(TIANOptionValue[3:steps], col = "blue")
# Add Result from BAW Approximation:
abline(h = 2.08599, lty = 2) #B&S
abline(h=2.1332, lty = 4, col="orange")
title(main = "EU put: Price Convergence (200 steps)")
data.frame(CRROptionValue, JROptionValue, TIANOptionValue)
legend(150, 2.08, legend=c("CRR", "JR", "TIAN", "B&S", "Convergence Price"),
       col=c("red", "blue", "green", "black", "orange"), lty=c(1,1,1,2,4), cex = 0.9)

------#Trinomial Tree

#   Write a function to compute the call and put price of an 
#   European or American style option using a trinomial tree
#   approach.

TrinomialTreeOption  = 
  function(AmeEurFlag, CallPutFlag, S, X, Time, r, b, sigma, n)
  {   # A function implemented by Diethelm Wuertz           
    
    # Description:
    #   Calculates option prices from the Trinomial tree model.
    
    # Arguments:
    #   AmeEurFlag - a character value, either "a" or "e" for 
    #       a European or American style option
    #   CallPutFlag - a character value, either "c" or "p" for 
    #       a call or put option
    #   S, X, Time, r, b, sigma - the usual option parameters
    #   n - an integer value, the depth of the tree
    
    # Value:
    #   Returns the price of the options.
    
    # Details:
    #   Trinomial trees in option pricing are similar to
    #   binomial trees. Trinomial trees can be used to 
    #   price both European and American options on a single 
    #   underlying asset.
    #   Because the asset price can move in three directions 
    #   from a given node, compared with only two in a binomial
    #   tree, the number of time steps can be reduced to attain
    #   the same accuracy as in the binomial tree. 
    
    # Reference:
    #   E.G Haug, The Complete Guide to Option Pricing Formulas
    #   Chapter 3.2
    
    # FUNCTION:
    
    # Settings:            
    OptionValue  =  rep(0, times=2*n+1)  
    
    # Call-Put Flag:
    if (CallPutFlag == "c") z  =  +1 
    if (CallPutFlag == "p") z  =  -1  
    
    # Time Interval: 
    dt  =  Time/n
    
    # Up-and-down jump sizes:
    u  =  exp(+sigma * sqrt(2*dt))
    d  =  exp(-sigma * sqrt(2*dt)) 
    
    # Probabilities of going up and down:  
    pu  =  ((exp(b * dt/2) - exp( -sigma * sqrt(dt/2))) / 
              (exp(sigma * sqrt(dt/2)) - exp(-sigma * sqrt(dt/2)))) ^ 2
    pd  =  (( exp(sigma * sqrt(dt/2)) - exp( b * dt/2)) / 
              (exp(sigma * sqrt(dt/2)) - exp(-sigma * sqrt(dt/2)))) ^ 2
    
    # Probability of staying at the same asset price level:
    pm  =  1 - pu - pd
    Df  =  exp(-r*dt)   
    for (i in 0:(2*n)) {
      OptionValue[i+1]  =  max(0, z*(S*u^max(i-n, 0) * 
                                       d^max(n*2-n-i, 0) - X))}
    for (j in (n-1):0) {
      for (i in 0:(j*2)) {
        # European Type:
        if (AmeEurFlag == "e") {
          OptionValue[i+1]  =  (
            pu * OptionValue[i+3] + 
              pm * OptionValue[i+2] + 
              pd * OptionValue[i+1]) * Df }
        # American Type:
        if (AmeEurFlag == "a") {
          OptionValue[i+1]  =  max((z*(S*u^max(i-j, 0) * 
                                         d ^ max(j*2-j-i, 0) - X)), (
                                           pu * OptionValue[i+3] + 
                                             pm * OptionValue[i+2] + 
                                             pd * OptionValue[i+1]) * Df) } } }
    TrinomialTree  =  OptionValue[1]
    
    # Return Value:
    TrinomialTree
  }

#EU CALL
TrinomialTreeOption(AmeEurFlag = "e", CallPutFlag = "c",S = 283, X = 250, 
                    Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)

#EU PUT
TrinomialTreeOption(AmeEurFlag = "e", CallPutFlag = "p",S = 283, X = 250, 
                    Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)

#US CALL
TrinomialTreeOption(AmeEurFlag = "a", CallPutFlag = "c",S = 283, X = 250, 
                    Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)

#US PUT
TrinomialTreeOption(AmeEurFlag = "a", CallPutFlag = "p",S = 283, X = 250, 
                    Time = 22/252, r = 0.0149, b = 0, sigma = 0.38610564, n = 5000)

