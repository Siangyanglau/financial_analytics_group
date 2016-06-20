
library(PerformanceAnalytics)
library(nloptr)
library(tseries)
# prepare the dataset removing the Funds of Funds fund indes
data("edhec")
edhec1<-edhec[,-13]

# N, mean and covariance monthly basis
N <- ncol(edhec1)
mu <- apply(edhec1,2,mean)
covMat <- cov(edhec1)

# create the stats table on yearly basis
stats <- cbind(mu*12*100, sqrt(diag(covMat))*sqrt(12)*100, mu/sqrt(diag(covMat))*sqrt(12))
colnames(stats) <- c("Ann. Returns (%)", "Volatility (%)", "Sharpe Ratio")

### equal weighted 
w_EW <- rep(1/N,N)
### inverse volatility
w_IV <- 1/sqrt(diag(covMat))/(sum(1/sqrt(diag(covMat))))




### Min-Var without weight constraints
eval_MinVar_RiskObjective <- function(w,covMat) {
  w <- as.matrix(w, col=1)
  vol <- sqrt( t(w) %*% covMat %*% w )[1]
  return(list("objective" = vol,
              "gradient"  = t(w) %*% covMat / vol ))
}

eval_MinVar_RiskConstraint <- function(w, covMat) {
  w <- as.matrix(w, col=1)
  vol <- sqrt( t(w) %*% covMat %*% w )[1]
  return(list("constraints" = t(w)%*%matrix(1,12)-1,
              "jacobian"    = matrix(1,1,12) ))
}



solve.MinVar <- function(covMat, constrained = TRUE) {
  wCap <- 100 # dummy max weight
  N <- ncol(covMat)
  initial <- 1/sqrt(diag(covMat))/(sum(1/sqrt(diag(covMat))))# assume initial is inverse risk
  
  upper <- rep(wCap,N)
  lower <- rep(ifelse(constrained,0,-wCap),N)
  
  result <- nloptr(x0 = initial,
                   eval_f = eval_MinVar_RiskObjective,
                   lb = lower,
                   ub = upper,
                   eval_g_eq = eval_MinVar_RiskConstraint,
                   opts = list("algorithm" = "NLOPT_LD_SLSQP",
                               "print_level"=0,
                               "maxeval" = 2000,
                               "check_derivatives" = FALSE,
                               "check_derivatives_print" = "all"),
                               covMat = covMat)
  
  #result$solution[abs(result$solution)<=0.005] = 0 # remove weights smaller than 0.5% in absolute value
  
  return(list(sol = result$solution,
              iter= result$iterations))
}

w_MinVar <- solve.MinVar(covMat = covMat)$sol













### ### Min-Var with weight constraints (add 3% and 20% constraints)
eval_MinVar_RiskObjective <- function(w,covMat) {
  w <- as.matrix(w, col=1)
  vol <- sqrt( t(w) %*% covMat %*% w )[1]
  return(list("objective" = vol,
              "gradient"  = t(w) %*% covMat / vol ))
}

eval_MinVar_RiskConstraint <- function(w, covMat) {
  w <- as.matrix(w, col=1)
  vol <- sqrt( t(w) %*% covMat %*% w )[1]
  return(list("constraints" = t(w)%*%matrix(1,12)-1,
              "jacobian"    = matrix(1,1,12) ))
}



solve.MinVar1 <- function(covMat, constrained = TRUE) {
  wCap <- 0.2 # dummy max weight
  N <- ncol(covMat)
    
  upper <- rep(wCap,N)
  lower <- rep(ifelse(constrained,0.03,-wCap),N)
  initial <- initial <- pmax(1/sqrt(diag(covMat))/(sum(1/sqrt(diag(covMat)))),lower)  # assume initail is inverse risk 

  result <- nloptr(x0 = initial,
                   eval_f = eval_MinVar_RiskObjective,
                   lb = lower,
                   ub = upper,
                   eval_g_eq = eval_MinVar_RiskConstraint,
                   opts = list("algorithm" = "NLOPT_LD_SLSQP",
                               "print_level"=0,
                               "maxeval" = 2000,
                               "check_derivatives" = FALSE,
                               "check_derivatives_print" = "all"),
                               covMat = covMat)
  
  #result$solution[abs(result$solution)<=0.005] = 0 # remove weights smaller than 0.5% in absolute value
  
  return(list(sol = result$solution,
              iter= result$iterations))
}

w_MinVar1 <- solve.MinVar1(covMat = covMat)$sol

"gradient" =t(sigma)/vol-(t(w)%*%covMat/vol)*(t(w)%*%sigma)[1]/vol[1]^2[1] ))}





### Max Diversification

eval_MaxDiv_RiskObjective <- function(w,covMat,sigma) {
  w <- as.matrix(w, col=1)
  sigma<-sqrt(diag(covMat))
  vol <- sqrt( t(w) %*% covMat %*% w )[1]
  div <- (t(w)%*%sigma)/vol
  return(list("objective" = -div,
  #### gradient of div need to check
              "gradient" =-t(sigma)/vol+(t(w)%*%covMat/vol)*(t(w)%*%sigma)[1]/vol[1]^2[1] ))}

eval_MaxDiv_RiskConstraint <- function(w, covMat,sigma) {
  w <- as.matrix(w, col=1)
  sigma<-sqrt(diag(covMat))
  vol <- sqrt( t(w) %*% covMat %*% w )[1]
  return(list("constraints" = t(w)%*%matrix(1,12)-1,
              "jacobian"    = matrix(1,1,12) ))
}



solve.MaxDiv <- function(covMat,sigma, constrained = TRUE) {
  wCap <- 100 # dummy max weight
  N <- ncol(covMat)
  sigma<-sqrt(diag(covMat))
  initial <- rep(1/N,N) # equally-weighted portfolio
  
  upper <- rep(wCap,N)
  lower <- rep(ifelse(constrained,0,-wCap),N)
  
  result <- nloptr(x0 = initial,
                   eval_f = eval_MaxDiv_RiskObjective,
                   lb = lower,
                   ub = upper,
                   eval_g_eq = eval_MaxDiv_RiskConstraint,
                   opts = list("algorithm" = "NLOPT_LD_SLSQP",
                               "print_level"=0,
                               "maxeval" = 2000,
                               "check_derivatives" = FALSE,
                               "check_derivatives_print" = "all"),
                               covMat = covMat, sigma=sigma)
  
  result$solution[abs(result$solution)<=0.005] = 0 # remove weights smaller than 0.5% in absolute value
  
  return(list(sol = result$solution,
              iter= result$iterations))
}

w_MaxDiv <- solve.MaxDiv(covMat = covMat,sigma=simga)$sol





### risk parity
eval_RP_RiskObjective <- function(w,covMat, tgtVol) {
  w <- as.matrix(w, col=1)
  rp<-sum(log(w))
  return(list("objective" = -rp,
              "gradient"  = -t(1/w) ))
}

eval_RP_RiskConstraint <- function(w, covMat, tgtVol) {
  w <- as.matrix(w, col=1)
  vol <- sqrt( t(w) %*% covMat %*% w )[1]
  return(list("constraints" = vol-tgtVol,
              "jacobian"    = t(w) %*% covMat / vol  ))
}



solve.RP <- function(covMat, tgtVol = 0.04, constrained = TRUE) {
  wCap <- 100 # dummy max weight
  N <- ncol(covMat)
  initial <- 1/sqrt(diag(covMat))/(sum(1/sqrt(diag(covMat))))# assume initial is inverse risk
  
  upper <- rep(wCap,N)
  lower <- rep(ifelse(constrained,0,-wCap),N)
  
  result <- nloptr(x0 = initial,
                   eval_f = eval_RP_RiskObjective,
                   lb = lower,
                   ub = upper,
                   eval_g_ineq = eval_RP_RiskConstraint,
                   opts = list("algorithm" = "NLOPT_LD_SLSQP",
                               "print_level"=0,
                               "maxeval" = 2000,
                               "check_derivatives" = FALSE,
                               "check_derivatives_print" = "all"),
                               covMat = covMat, tgtVol = tgtVol)
  
  #result$solution[abs(result$solution)<=0.005] = 0 # remove weights smaller than 0.5% in absolute value
  
  return(list(sol = result$solution,
              iter= result$iterations))
}

w_RP <- solve.RP(covMat = covMat,0.04/sqrt(12))$sol
w_RP_scaled <- w_RP/sum(w_RP)


### A matrix (12X7) with 12 assets and 7 weighting schemes
 
weights<-cbind(w_EW,w_IV,w_MinVar,w_MinVar1,w_MaxDiv,w_RP,w_RP_scaled)

# name the columns 
colnames(weights) <- c("EW","IV", "MinVar_unconstrained", "MinVar_constrained", "MD","RP","RP_scaled")

####barplot
library(ggplot2)
library(cowplot)
weigt<-data.frame(weights)
assets<-rownames(weigt)
ew<-ggplot(data=weigt, aes(x=assets, y=weigt$EW))+geom_bar(stat="identity")+ylim(0,1)
#text = element_text(size=20),                                                                                            axis.text.x = element_text(angle=90, vjust=1)) 
iv<-ggplot(data=weigt, aes(x=assets, y=weigt$IV))+geom_bar(stat="identity")+ylim(0,1)
mv_unc<-ggplot(data=weigt, aes(x=assets, y=weigt$MinVar_unconstrained))+geom_bar(stat="identity")+ylim(0,1)
mv_c<-ggplot(data=weigt, aes(x=assets, y=weigt$MinVar_constrained))+geom_bar(stat="identity")+ylim(0,1)
md<-ggplot(data=weigt, aes(x=assets, y=weigt$MD))+geom_bar(stat="identity")+ylim(0,1)
rp<-ggplot(data=weigt, aes(x=assets, y=weigt$RP))+geom_bar(stat="identity")+ylim(0,1)
rp_s<-ggplot(data=weigt, aes(x=assets, y=weigt$RP_scaled))+geom_bar(stat="identity")+ylim(0,1)

plot_grid(ew, iv, mv_unc,mv_c,md,rp,rp_s, labels=c("EW", "IV","MinVar_Unc","MinVar_C","MD","RP","RP_Scaled"), ncol = 2, nrow = 4)

### cumulative return plots
#porfolio return
pr<-edhec1%*%weights
#cumulative returns (1+r1)*(1+r1)....(1+r152)
w_t_EW <- cumprod(na.omit(pr[,1])+1)
plot(w_t_EW, main="Cumulative Return of EW Portfolio")

w_t_IV <- cumprod(na.omit(pr[,2])+1)
plot(w_t_IV, main="Cumulative Return of IV Portfolio")

w_t_MinVar_Uncon <- cumprod(na.omit(pr[,3])+1)
plot(w_t_MinVar_Uncon, main="Cumulative Return of MinVarUnconstrained Portfolio")

w_t_MinVar_c <- cumprod(na.omit(pr[,4])+1)
plot(w_t_MinVar_c, main="Cumulative Return of MinVarConstrained Portfolio")

w_t_MD <- cumprod(na.omit(pr[,5])+1)
plot(w_t_MD, main="Cumulative Return of MD Portfolio")

w_t_RP <- cumprod(na.omit(pr[,6])+1)
plot(w_t_RP, main="Cumulative Return of RP Portfolio")

w_t_RP_s <- cumprod(na.omit(pr[,7])+1)
plot(w_t_s, main="Cumulative Return of RP_Scaled Portfolio")

### STATS

#1 annualised portfolio returns
mu_pr<- t(weights) %*% mu * 12 * 100
mu_pr_alt<-apply(pr,2, mean)*12*100 # another way to calculated mean portfolio return


#2 volatility 
mu_P <- t(weights) %*% mu * 12 * 100
vol_P<- sqrt(diag(t(weights) %*% covMat %*% weights)) * sqrt(12) * 100 
volMax_P <- t(weights) %*% sqrt(diag(covMat)) * sqrt(12) * 100

#3 diversification ratio
dr_P<-volMax_P/vol_P

#4 annulised sharp ratio
sharpeR_P<-mu_P/vol_P

#5 sortino ratio

#Vol_down
# downsidederivation
# add date to montly portfolio return 
date<-index(edhec1)
return<-data.frame(pr)
rownames(return)<-date

vol_down <- DownsideDeviation(return,method='subset')
# alternative definition by dividing with the number of all returns (not just the negative)
vol_down_alt <- DownsideDeviation(return,method='full')

Sortino <- t(mu_P)/vol_down
Sortino_alt <-t(mu_P)/vol_down_alt
Sortino_fct <- SortinoRatio(mu_P, MAR = 0)


#6 Value-at-Risk at 5% with Gaussian distiribution
VaR_5_gaussian <- -qnorm(0.05,mean=mu_P,sd=vol_P)

# 7 Conditional VaR
CVaR_5_gaussian_EW <- -integrate(function(x) x*dnorm(x,mean=mu_P[1],sd=vol_P[1]),lower=-Inf,upper=-VaR_5_gaussian[1])$value

CVaR_5_gaussian_IV <- -integrate(function(x) x*dnorm(x,mean=mu_P[2],sd=vol_P[2]),lower=-Inf,upper=-VaR_5_gaussian[2])$value

CVaR_5_gaussian_MinVar_unc <- -integrate(function(x) x*dnorm(x,mean=mu_P[3],sd=vol_P[3]),lower=-Inf,upper=-VaR_5_gaussian[3])$value

CVaR_5_gaussian_MinVar_c <- -integrate(function(x) x*dnorm(x,mean=mu_P[4],sd=vol_P[4]),lower=-Inf,upper=-VaR_5_gaussian[4])$value

CVaR_5_gaussian_MD <- -integrate(function(x) x*dnorm(x,mean=mu_P[5],sd=vol_P[5]),lower=-Inf,upper=-VaR_5_gaussian[5])$value

CVaR_5_gaussian_RP <- -integrate(function(x) x*dnorm(x,mean=mu_P[6],sd=vol_P[6]),lower=-Inf,upper=-VaR_5_gaussian[6])$value

CVaR_5_gaussian_RP_scaled <- -integrate(function(x) x*dnorm(x,mean=mu_P[7],sd=vol_P[7]),lower=-Inf,upper=-VaR_5_gaussian[7])$value

CVaR_5_gaussian_P<-c(CVaR_5_gaussian_EW, CVaR_5_gaussian_IV,CVaR_5_gaussian_MinVar_unc, CVaR_5_gaussian_MinVar_c, CVaR_5_gaussian_MD, CVaR_5_gaussian_RP, CVaR_5_gaussian_RP_scaled)

CVaR_5_gaussian_P<-matrix(CVaR_5_gaussian_P,1)
colnames(CVaR_5_gaussian_P)<-c("EW", "IV","MinVar_Unc","MinVar_C","MD","RP","RP_Scaled")

#8 skewness
skew <- skewness(return/100)

#9 maximun draw down
DD <- table.Drawdowns(na.omit(return))
MDD <- maxDrawdown(return)
chart.Drawdown(return)


# stats matrix
stats_P <- rbind(t(mu_P),t(vol_P), t(dr_P), t(sharpeR_P), Sortino, t(VaR_5_gaussian), CVaR_5_gaussian_P, skew, MDD)
colnames(stats_P) <- c("Ann. Returns (%)", "Volatility (%)", "Diversification Ratio","Sharpe Ratio","Sortino Ration", "Value_at_Risk(%)","Conditional Value_at_Risk(%)","Skewness", "MDD")

