```{r monte_carlo_option_pricing_}
#########################################
# Step 0 - Define initial/valuation input
#########################################
#Set seed for random numbers
set.seed(12345)
# Simulation period 
T<-0.5
# Step in y with 360 day-count convention
tStep<-1/(360*2)
# Scenarios
nPaths<-10000
#Evaluation times in units of year  
time<-seq(0,T,tStep)

# Parameters (S0, vola, drift)	
start_value<-100
# yearly vola
vola<-0.25
# Annualized rate of interest
r <-0.01
# Annualized cost-of-carry rate
b <-0.01
# Barrier Value
H <- 75
# yearly drift
drift<-r
# Strike price
strike<-90


#######################################################################################
# Step 1 - Simulate random paths with GBM process for the spot value until maturity T.
#######################################################################################
#Definition of the RF box	
RF_GBM<-matrix(0,nrow=nPaths,ncol=length(time))

#We generate random numbers N(0,1)
#Always better to generate the innovations in a vectorial fashion and then "distribute" them
pass_rand<-rnorm(nPaths*(length(time)-1))
pass_rand<-matrix(pass_rand,nrow=nPaths,ncol=(length(time)-1))
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

par(mfrow=c(1,1))
matplot(time,t(RF_GBM)[,1:100],type="l",xlab="t(y)",ylab="RF",main="GBM process")
points(time,apply(t(RF_GBM),1,mean),type="l",col=1,lwd=4)
points(time,apply(t(RF_GBM),1,quantile,0.05),type="l",col=2,lwd=4)
points(time,apply(t(RF_GBM),1,quantile,0.95),type="l",col=2,lwd=4)
grid()

#########################################################
# Step 2 - Decide whether each path is knocked in or out 
# and thus determine the final payoff of each path.
#########################################################
#Vanilla Opiton
RF_GBM_vanilla <- RF_GBM

#Exotic Option - down-and-in
RF_GBM_di <- RF_GBM[RF_GBM[,ncol(RF_GBM)]<=H,]

#Exotic Option - down-and-out
RF_GBM_do <- RF_GBM[RF_GBM[,ncol(RF_GBM)]>H,]
############################################
# Step 3 - Calculate the payoff at maturity. 
############################################
#Vanilla Opiton
payoff_vanilla <- pmax(RF_GBM_vanilla[,ncol(RF_GBM_vanilla)] - strike,0)


#Exotic Option - down-and-in
payoff_di <- pmax(RF_GBM_di[,ncol(RF_GBM_di)] - strike,0)


#Exotic Option - down-and-out
payoff_do <- pmax(RF_GBM_do[,ncol(RF_GBM_do)] - strike,0)

###################################################################
# Step 4 - Discount all future payoff distribution of values to t0.
###################################################################
#Vanilla Opiton
discounted_payoff_vanilla <- exp(-r*T)*payoff_vanilla

#Exotic Option - down-and-in
discounted_payoff_di <- exp(-r*T)*payoff_di

#Exotic Option - down-and-out
discounted_payoff_do <- exp(-r*T)*payoff_do

################################################################
# Step 5 - Calculate the mean of discounted payoff distribution.
################################################################
#Vanilla Opiton
price_vanilla <- mean(discounted_payoff_vanilla)

#Exotic Option - down-and-in
price_di <- mean(discounted_payoff_di)

#Exotic Option - down-and-out
price_do <- mean(discounted_payoff_do)



#down-in
if(any(RF_GBM < H)){
  payoff_downin <-0
}else{
  payoff_downin <-pmax (RF_GBM_at_maturity - strike, 0)
}

payoff_downin_price <- mean(exp(-r*T)*(payoff_downin))


#downout
if(any(RF_GBM > H)){
  payoff_downout <- pmax (RF_GBM_at_maturity - strike, 0)
}else{
  payoff_downout <- 0
}

payoff_downout_price <- mean(exp(-r*T)*(payoff_downout))
#vanilla
RF_GBM_at_maturity<-RF_GBM[,ncol(RF_GBM)]
payoff_vanilla<-pmax(RF_GBM_at_maturity-strike,0)



payoff_vanilla_price<-mean(exp(-r*T)*(payoff_vanilla))

