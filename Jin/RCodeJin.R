library(quantmod)
library(fOptions)

rm(list=ls()) 

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
tStep<-10/260
# Scenarios
nPaths<-10000
#Evaluation times in units of year  
time<-seq(0,T,tStep) 

set.seed(12345)
# Parameters (S0, vola, drift) 
start_value<-S0
vola <- 0.3  #implied volatility
drift <- 0.01


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
pass_drift<-matrix(rep(exp(-0.5*time*vola^2+time*drift),nPaths),nrow=nPaths,ncol=length(time),
                   byrow=TRUE)
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

#ans = 2.190825


##Question 1 (iii)

#Use your Black Scholes function and the paths derived in the previous step to
#derive MtM (Mark- to-Market) distributions for the call option over time.
#Mark to market (MtM) is a measure of the fair value of accounts that can change over time


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
C0<-exp(-r*T)*mean(pmax(S[nrow(S),] - strike,0))

#MtM shouldn't be discounted to present value - i think, so it should be:
C0<-mean(pmax(S[nrow(S),] - strike,0))

#using my own code of RF_GBM
#need to change RF_GBM, which is discounted steps to expected future payoffs
#RF_GBM is the transpose of S, and also reverse the columns. 
#because MtM looks at time to maturity, and RF_GBM is time away from today. 

#transform RF_GBM into the right form. yes, now it is starting from today to end of path.
RF_GBMT <- t(RF_GBM)

C<-S*0
for (i in 1:nrow(RF_GBMT)){
  C[i,]<-GBSOption(TypeFlag="c",S=RF_GBMT[i,],X=strike,Time=time[i],r=r,b=r,sigma=sigma)@price
}


#price of call option is 3.368507, same as above

#yess!!!!!! i hope this is right.
# plot me results
defaultT<-seq(0,T,dt)
par(mfrow=c(1,2))
matplot(defaultT,C[,1:100],type="l",ylab="Call Option - Black Scholes, price [$]",
        xlab="Default Dates [y]",main="Closed Form",ylim=c(min(C),max(C)))


###############################################################################

#(iv) the MtM distribution is C - code is in CVA Example
#calculate expected exposure
#calculate he Potential Future Exposure (PFE) at 95% confidence level 
# calculate the CVA price using counterparty spread curve:

#########################################################################################################
# Calculate EE, PFE, CVA
#########################################################################################################

#MTM time is time away from today, not time to maturity
tStep<-seq(0,T,dt)

MtM <-C

EE<-apply(pmax(MtM,0),1, mean)
PFE<-apply(pmax(MtM,0),1,quantile,0.95)
par(mfrow=c(1,2))
plot(tStep,EE,type="l")
plot(tStep,PFE,type="l")

# here the CDS curve info
cdsT<-seq(1,5,1)
cds<-c(92,104,112,117,120)
par(mfrow=c(1,1))
plot(cdsT,cds,type="l")

#Spline function
#NB: You need the CDS spread at each repricing date, i.e. from t0 up to T with 10d step. 
#Use the above CDS spread curve and a natural spline [R function: spline() ] 
#to get the CDS spreads for CVA pricing. The Loss Given Default (LGD) equals 40%.

t<-tStep
r<-0.01
DF<-exp(-r*t)
LGD<-0.4
# Basel formula
CDS_curve<-spline(cdsT, cds, xout=t)
temp<-c()
for (i in 2:length(t)){
  temp<-c(temp,max(0,exp(-(CDS_curve$y[i-1]/10000*CDS_curve$x[i-1])/LGD)-
                     exp(-(CDS_curve$y[i]/10000*CDS_curve$x[i])/LGD))*
            (EE[i-1]*DF[i-1]+EE[i]*DF[i])*0.5)
}

CVABasel<-LGD*sum(temp)
print(paste("The CVA (using Basel formula) for the uncollateralized trade equals ",CVABasel))

#answer = 0.12130784488567
#YES!!!

#CVA is a price. It is now a considerable part of the PnL of any financial institution
#CVA is the cost of buying protection on the counterparty that pays the 
#portfolio value in case of default

#########################################################################################################
#(v) new option!

S0<-10
sigma<-0.3
mu<-0
T<-5
r<-0.01

#each time units d/dt = #Evaluation times in units of year 
dt<-10/260  #to get nSteps = 10, 130 is chosen. 260 is the number of business days
nPaths<-10000
nSteps<-T/dt  # 5/(130/260) = 10
#calculate MtM2

#this is the first dw used for the first option.
pass_rand<-rnorm(nPaths*(length(time)-1)) #generate 10k numbers
pass_rand<-matrix(pass_rand,nrow=nPaths,ncol=(length(time)-1)) #distribute it to matrix
pass_rand<-scale(pass_rand,center = TRUE,scale = TRUE)

#to find the second dw2




# generate correlated stock price paths, with correlation of 0.7
# use mtvnorm to generate 10,000 values [dW1dW2 = ρdt]

dt<-10/260
rho = 0.7

require(MASS)
out <- mvrnorm(10000, mu = c(0,0), 
               Sigma = matrix(c(1,rho,rho,1), ncol = 2), 
               empirical = TRUE)
cor(out[,1], out[,2])

dw1<-matrix(out[,1],nrow=nPaths,ncol=(length(time)-1))
dw2<-matrix(out[,2],nrow=nPaths,ncol=(length(time)-1))

S3<-S0*exp(apply((r-sigma^2/2)*dt + sigma*sqrt(dt)*dw1,1,cumsum))
S3<-rbind(S0,S3)

C3<-S3*0
for (i in 1:nrow(S3)){
  C3[i,]<-GBSOption(TypeFlag="c",S=S3[i,],X=12,Time=time[i],r=r,b=r,sigma=sigma)@price
}

S4<-S0*exp(apply((r-sigma^2/2)*dt + sigma*sqrt(dt)*dw2,1,cumsum))
S4<-rbind(S0,S4)

C4<-S4*0
for (i in 1:nrow(S4)){
  C4[i,]<-GBSOption(TypeFlag="c",S=S4[i,],X=14,Time=time[i],r=r,b=r,sigma=sigma)@price
}

MtM3 <-C3+C4

defaultT<-seq(0,T,dt)
par(mfrow=c(1,1))
matplot(defaultT,C3[,1:100],type="l",ylab="Call Option - Black Scholes, price [$]",
        xlab="Default Dates [y]",main="Closed Form",ylim=c(min(C),max(C)))




