library(quantmod)
library(fOptions)

#### Question 1 (i) #####
#### Black-scholes equation

#inputs for BS formula

S0 <- 10  #starting equity price
K <-12    # strike price
vola <- 0.3  #implied volatility
r<- 0.01 #the risk free interest rate
T<- 5    #maturity of option in years 

#the BS formula manually = ans = 2.157715

d1 <- 1/(vola*sqrt(T))*(log(S0/K)+(r+0.5*vola^2)*(T))
d2 <- d1 - vola*sqrt(T)

price = pnorm(d1,0,1)*S0 - pnorm(d2)*K*exp(-r*T)

#BS using fOptions
spot<-10
r<-0.01
T<-5
vola<-0.3
strike<-12

# Price with Black Scholes formula (load the fOptions R standard package)
# ans = 2.157716
library(fOptions)
GBSOption(TypeFlag = "c",S = spot,X = strike,Time = T,r = r, b = r,sigma = vola)@price


#### Question 1 (ii) #####
#### Monte Carlo equation

# Simulation period 
T
# Step in 10days with 260 day-count convention
tStep<-10/260  #does he mean 360 days?
# Scenarios
nPaths<-10000
#Evaluation times in units of year  
time<-seq(0,T,tStep)

set.seed(12345)
# Parameters (S0, vola, drift) 
start_value<-S0
vola <- 0.3  #implied volatility
drift <- 0.1


#Definition of the RF box  
#matrix to show path
RF_GBM<-matrix(0,nrow=nPaths,ncol=length(time))

#We generate random numbers N(0,1)
#Always better to generate the innovations in a vectorial fashion and then "distribute" them
#rnorm(n, mean = , sd = ) is used to generate n normal random numbers with arguments mean and sd
#scale, with default settings, will calculate the mean and standard deviation of the entire vector, 
#then "scale" each element by those values by subtracting the mean and dividing by the sd. 
#(If you use  scale(x, scale=FALSE), it will only subtract the mean but not divide by the std deviation.)

pass_rand<-rnorm(nPaths*(length(time)-1)) #generate 10k numbers
pass_rand<-matrix(pass_rand,nrow=nPaths,ncol=(length(time)-1)) #distribute it to matrix
pass_rand<-scale(pass_rand,center = TRUE,scale = TRUE)

#Adding the starting value
pass_rand<-cbind(rep(0,nPaths),pass_rand)
#Accrueting the innovations
pass_rand<-t(apply(pass_rand,1,cumsum))
#drift_matrix: x=x0*exp(mu*t-vola^2*t/2)*exp(vola*N(0,1)*sqrt(t))
pass_drift<-matrix(rep(exp(-0.5*time*vola^2+time*drift),nPaths),nrow=nPaths,ncol=length(time), byrow=TRUE)
#Generation of the paths (scaling the random component by tStep)
RF_GBM<-start_value*exp(pass_rand*vola*sqrt(tStep))*pass_drift
rownames(RF_GBM)<-paste("history",as.character(seq(1,nPaths,1)),sep="")
colnames(RF_GBM)<-paste(as.character(time),"y",sep="")  

# Here we make checks
par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
check_sd<-apply(RF_GBM,2,sd)/sqrt(time)
plot(time,check_sd,ylim = c(0.15,0.40))
abline(h=vola,col=2,lwd=2)
check_mean<-apply(RF_GBM,2,mean)
plot(time,check_mean,type="l",lwd=4)  
abline(a = start_value ,b = drift,col=2,lwd=2)

par(mfrow=c(1,1))
matplot(time,t(RF_GBM)[,1:100],type="l",xlab="t(y)",ylab="RF",main="GBM process")
points(time,apply(t(RF_GBM),1,mean),type="l",col=1,lwd=4)
points(time,apply(t(RF_GBM),1,quantile,0.05),type="l",col=2,lwd=4)
points(time,apply(t(RF_GBM),1,quantile,0.95),type="l",col=2,lwd=4)
grid()

#compute payoff at maturity (slide 21)
#maturity = 5yrs
#S(0) = 10, K=12

#find the colum names where expiry is 5 and glue expiry with y=years.
expiry<-5 
ind<-which(colnames(RF_GBM)==paste0(expiry,"y"))
K<-12
ST<-RF_GBM[,ind] #10thousand values for ending paths
plot(ST,pmax(0,ST-K))

#discount all future payoffs(slide 22)
discount<-exp(-r*T)*pmax(0,ST-K)
#Calculate the mean of discounted payoff distribution
mean(discount)

#ans = 6.121821

##Question 1 (iii)

#Use your Black Scholes function and the paths derived in the previous step to
#derive MtM (Mark- to-Market) distributions for the call option over time.
#Mark to market (MTM) is a measure of the fair value of accounts that can change over time

###############################################################################

rm(list=ls())
gc()

###############################################################################
# source libraries
library(fOptions)

###############################################################################
# define initial/valuation input
S0<-10
sigma<-0.3
mu<-0
T<-5
r<-0.01

#each time units d/dt = #Evaluation times in units of year 
dt<-10/260  #to get nSteps = 10, 130 is chosen. 260 is the number of business days
nPaths<-10000
nSteps<-T/dt  # 5/(130/260) = 10

###############################################################################
# simulate Stock random paths
dw<-matrix(rnorm(nPaths*nSteps,0,1),nPaths,nSteps)
S<-S0*exp(apply((r-sigma^2/2)*dt + sigma*sqrt(dt)*dw,1,cumsum))
S<-rbind(S0,S)
matplot(S[,1:100],type="l")

#Note that W, and consequently its infinitesimal increment dW, 
#represents the only source of uncertainty in the price history of the stock.

###############################################################################
# calculate MtM distributions for a Call option - BS closed formula

strike<-12
time<- T - seq(0,T,dt)   #working backwards as advised in the project notes. i should have a reading every 10 days

# from lecture: the expected payoff in the future is the MtM distribution
C<-S*0
for (i in 1:nrow(S)){
    C[i,]<-GBSOption(TypeFlag="c",S=S[i,],X=strike,Time=time[i],r=r,b=r,sigma=sigma)@price
}

#using my own code of RF_GBM
#need to change RF_GBM, which is discounted steps to expected future payoffs
#RF_GBM is the transpose of S

RF_GBMT <- t(RF_GBM)

C<-S*0
for (i in 1:nrow(RF_GBMT)){
  C[i,]<-GBSOption(TypeFlag="c",S=RF_GBMT[i,],X=strike,Time=time[i],r=r,b=r,sigma=sigma)@price
}

#yess!!!!!! i hope this is right.
#plot to check, should be symmetrical
# plot me results
defaultT<-seq(0,T,dt)
par(mfrow=c(1,2))
matplot(defaultT,C[,1:100],type="l",ylab="C - BS, price [$]",
        xlab="Default Dates [y]",main="Closed Form",ylim=c(min(C),max(C)))

###############################################################################

#(iv)

