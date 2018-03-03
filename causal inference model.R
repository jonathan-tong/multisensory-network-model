# clear workspace
rm(list = ls(all = T))


# load packages
library(ez)
library(plyr)
library(ggplot2)

#Auditory position
a = 10
#visual position
v = -90

#Auditory sigma
siga = 8.1
vara = siga^2

#Visual sigma
sigv = 1.7
varv = sigv^2


#causal prior
pcom = 0.5

####calculate likelihoods####
xv = rnorm(10000, v, sigv)
xa = rnorm(10000, a, siga)
s = -90:90
prior = 1/181 #uniform prior over positions


#find estimates for all 10,000 samples
vestimates = rep(NA,10000)
aestimates = rep(NA,10000)

for(i in 1:10000){
  
  #calculate fully integrated estimate
  sac1 = (xa[i]/vara + xv[i]/varv)/(1/vara + 1/varv)
  sac2 = xa[i]
  
  za = (-1*(xa[i]-s)^2)/(2*vara)
  zv = (-1*(xv[i]-s)^2)/(2*varv)
  likea = (1/(siga*sqrt(2*pi)))*exp(za)
  likev = (1/(sigv*sqrt(2*pi)))*exp(zv)
#find single cause likelihood
likec1 = sum(likea*likev*prior)
likec2 = sum(likea*prior)*sum(likev*prior)
 
 postc1 = (likec1*pcom)/((likec1*pcom) + (likec2*(1-pcom)))
 postc2 = 1-postc1
 
 #model averaging
 aestimates[i] = postc1 * sac1 + postc2 * sac2
  vestimates[i] = postc1 * sac1 + postc2 * v
}

sA = mean(aestimates)
sA
sV = mean(vestimates)
sV
