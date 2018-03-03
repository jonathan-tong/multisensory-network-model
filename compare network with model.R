# clear workspace
rm(list = ls(all = T))


# load packages
library(ez)
library(plyr)
library(ggplot2)

#####INPUT PARAMETERS####

#auditory gain
ga = 140
#auditory rise
m = 20
#visual gain
gv = 80
#visual tuning width
sig = 20 

L = 301

####WEIGHT PARAMETERS####
A = 2
V = 5
am = 1
vm = 2

####unit positions####
#input units
alpos = -150:150
arpos = -150:150
vpos = -150:150
#pooling units
alposp = -150:150
arposp = -150:150
vposp = -150:150
mposp = -150:150
#reconstruction units
rposr = -150:150
lposr = -150:150
vposr = -150:150

nl = length(alpos)
nr = length(arpos)
nv = length(vpos)

####CREATE WEIGHT MATRICES####
#initialize the weight matrix of auditory (left and right) and visual units, matrix to be the size 301 by 301
wl = matrix(data = NA, nrow = length(alpos), ncol = length(alposp))
wr = matrix(data = NA, nrow = length(arpos), ncol = length(arposp))
wv = matrix(data = NA, nrow = length(vpos), ncol = length(vposp))

wlm = matrix(data = NA, nrow = length(alpos), ncol = length(mposp))
wrm = matrix(data = NA, nrow = length(arpos), ncol = length(mposp))
wvm = matrix(data = NA, nrow = length(vpos), ncol = length(mposp))

#left auditory weight matrix
n = 1
for(j in alposp){
  p = 1
  for(i in alpos){
    wl[n,p]=A/(nl*(1+exp(-(i-j)/m)))
    
    p = p+1
  }
  n = n+1
}

#right auditory weight matrix
n = 1
for(j in arposp){
  p = 1
  for(i in arpos){
    wr[n,p]=A/(nr*(1+exp((i-j)/m)))
    
    p = p+1
  }
  n = n+1
}

#Visual weight matrix
n = 1
for(j in vposp){
  p = 1
  for(i in vpos){
    delta = abs(i-j)
    delta1 = ifelse(delta<L/2, delta, L-delta)
    wv[n,p]=(V/(sig*sqrt(2*pi)))*exp(-(delta1)^2/(2*sig^2))
    p = p+1
  }
  n = n+1
}

##multisensory
#visual to multisensory weights
n = 1
for(j in mposp){
  p = 1
  for(i in vpos){
    wvm[n,p]=(vm/(sig*sqrt(2*pi)))*exp(-(i-j)^2/(2*sig^2))
    p = p+1
  }
  n = n+1
}

#auditory right to multisensory weights
n = 1
for(j in mposp){
  p = 1
  for(i in arpos){
    wrm[n,p]=am/(nr*(1+exp((i-j)/m)))
    p = p+1
  }
  n = n+1
}

#auditory left to multisensory weights
n = 1
for(j in mposp){
  p = 1
  for(i in alpos){
    wlm[n,p]=am/(nl*(1+exp(-(i-j)/m)))
    p = p+1
  }
  n = n+1
}

####DETERMINE UNIT ACTIVITIES####

#STIMULUS LOCATIONS#
xa =0  #auditory location

aest = rep(NA, 91)
vest = rep(NA, 91)

visualpos = seq(-90,90,2)
count = 1
for (j in visualpos){
xv = j  #visual location

#initial multisensory input#

mult = 10.5

####input activity####
#determine input activity given auditory position

ala = rep(1,301) #adaptation parameter left
ara = rep(1,301) #adaptation parameter right

alinput = (ala*ga)/(1+exp((xa-alpos)/m))
arinput = (ara*ga)/(1+exp(-(xa-arpos)/m))


#determine input activity given visual position
#wrap around
delta = abs(xv-vpos)
delta1 = ifelse(delta<L/2, delta, L-delta)

vinput = gv*exp(-(delta1)^2/(2*sig^2))

####pooling activities####
#auditory pooling units#
alpact = wl %*% alinput
arpact = wr %*% arinput

apact = arpact + alpact

#visual pooling units#
vpact = wv %*% vinput

#multisensory pooling units#

mpact =  (wvm %*% vinput) + (wrm %*% arinput) + (wlm %*% alinput) + mult

#normalization
sumpact = 1+((sum(exp(apact)) + sum(exp(vpact)) + sum(exp(mpact))) / (length(apact) + length(vpact) + length(mpact)))

apact = exp(apact)/sumpact
vpact = exp(vpact)/sumpact
mpact = exp(mpact)/sumpact

####reconstruction activities####

#auditory reconstruction left
lrho = (t(wl) %*% apact) + (t(wlm) %*% mpact)

#auditory reconstruction right
rrho = (t(wr) %*% apact) + (t(wrm) %*% mpact)

arho=lrho*rrho
Sa = rposr[which.max(arho)]
Sa

#visual reconstruction
vrho = (t(wv) %*% vpact) + (t(wvm) %*% mpact)
Sv = vposr[which.max(vrho)]
Sv

aest[count] = Sa
vest[count] = Sv

count = count + 1
print(count)
}

#####CAUSAL INFERENCE####

#Auditory position
a = 0

audests = rep(NA, 91)
visests = rep(NA, 91)

visualpos = seq(-90,90,2)
count = 1
for (j in visualpos){
#visual position
v = j
v
#Auditory sigma 
siga = 8.1
vara = siga^2

#Visual sigma
sigv = 1.8
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
  
  za = (-(xa[i]-s)^2)/(2*vara)
  zv = (-(xv[i]-s)^2)/(2*varv)
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

audests[count] = sA
visests[count] = sV
count = count + 1
print(v)
}

plot(visualpos, aest)
lines(visualpos, audests)

#calculate RMSE for auditory estimates
RMSEa = sqrt((sum((aest - audests)^2))/length(aest))
RMSEa

plot(visualpos, vest)
lines(visualpos,visests)

#calculate RMSE for visual estimates
RMSEv = sqrt((sum((vest - visests)^2))/length(vest))
RMSEv