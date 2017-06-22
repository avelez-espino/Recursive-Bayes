#########################################################################
#Model to estimate stock recruit relationships using Recursive Bayes Ricker model 

#########################################################################
library(R2WinBUGS)

#Pick the stock for which to run the models
#Bowron (s=4)
#Chilko (s=7)
#Cultus lake (s=11)
#Pitt ES  (s=18)
#Shuswap - ES - Scotch & Seymour (s=8)
#Weaver - Harrison upstream Lake (s=13)
#Nadina-Francois (s=17)
s <- 18

#Load the stock recruit data
SR_data <- read.table("CU Sockdat Jan 13 2014 (final).csv", header=T, sep=",")    #read in the data file
SR_data$RS<- SR_data$rec/SR_data$eff   #Calculate recruits per spawner
SR_data$sex<- SR_data$ETS/SR_data$eff
SR_data <- subset(SR_data, SR_data$eff != "NA")    #remove records with missing data
SR_data <- subset(SR_data, SR_data$rec != "NA")    #remove records with missing data
SR_data_short <- subset(SR_data, SR_data$yr > 1986)    #read in the data file with a shortened time series


#stock number for which the analysis is run
stock <- s

#select the stock recruit data for the stock of interest
R_Obs<-SR_data$rec[SR_data$stk==stock]
S<-SR_data$eff[SR_data$stk==stock] 
N<-length(R_Obs)

#select the stcok recruit data for the stock of interest in case the time series is short
R_Obs_short<-SR_data_short$rec[SR_data_short$stk==stock]
S_short<-SR_data_short$eff[SR_data_short$stk==stock]
N_short<-length(R_Obs_short)

#Create locations to store the model results
#alpha and beta results for the 10,000 MCMC values for the model 
#in case of variable annual alpha, only the average of the last 4 years is stored
alpha.values<-matrix(nrow=10000, ncol=1)# 
beta.values<-matrix(nrow=10000, ncol=1)# 
log.beta.values<-matrix(nrow=10000, ncol=1)
#median alpha and beta results stored for each year
a.values<-matrix(nrow=N, ncol=1)# 
b.values<-matrix(nrow=N, ncol=1)# 

#input for the winbugs model
data<-list("S", "R_Obs", "N")

parameters <- c("alpha", "beta", "tau_R")
inits <- c("inits1.txt", "inits2.txt")
initsKalman <- c("initsKalman1.txt", "initsKalman2.txt")

#Select the number of MCMC samples to take
n.chains=2  #Number of MCMC chains to monitor
n.burnin=5000    #Burn-in to be discarded
n.samples=10000   #Number of MCMC iterations run
n.thin=1      #thinning of the MCMC chain

#Run the WinBUGS model Recursive Bayes Ricker
results <- bugs(data, initsKalman, parameters, "RecursiveBayesUniformative.txt", n.chains, n.burnin, n.iter=n.samples,  n.thin, DIC=T, debug=F, save.history=F)

attach.bugs(results)

for (k in 1:10000){
   alpha.values[k,6]<-mean(alpha[k,(N-3):N])
}
beta.values[,6]<-beta

for (j in 1:N){
    a.values[j,6]<- quantile(alpha[,j], probs=c(0.5))
    b.values[j,6]<- quantile(beta, probs=c(0.5))
}

#input for the winbugs model
data<-list("S", "R_Obs", "N", "mu.beta","tau.beta")
parameters <- c("alpha", "beta", "tau_R")


 detach.bugs()

 detach()

