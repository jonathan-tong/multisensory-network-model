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
sig = 20 #80

L = 301

####WEIGHT PARAMETERS####
A = 2
V = 5 #4.335
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
xv = 20  #visual location

#initial multisensory input, mu#

mult = 10.5

####input activity####
#determine input activity given auditory position

ala = rep(1,301) #initialize adaptation parameter left
ara = rep(1,301) #initialize adaptation parameter right

alinput = (ala*ga)/(1+exp((xa-alpos)/m))
arinput = (ara*ga)/(1+exp(-(xa-arpos)/m))

plot(alpos,alinput)
lines(arpos,arinput)

#determine input activity given visual position
#wrap around
delta = abs(xv-vpos)
delta1 = ifelse(delta<L/2, delta, L-delta)

vinput = gv*exp(-(delta1)^2/(2*sig^2))

plot(vpos,vinput)

####pooling activities####
#auditory pooling units#
alpact = wl %*% alinput
arpact = wr %*% arinput

apact = arpact + alpact

#visual pooling units#
vpact = wv %*% vinput

#multisensory pooling units#

mpact =  (wvm %*% vinput) + (wrm %*% arinput) + (wlm %*% alinput) + mult

#activation non-linear exponential
sumpact = 1+((sum(exp(apact)) + sum(exp(vpact)) + sum(exp(mpact))) / (length(apact) + length(vpact) + length(mpact)))
  
apact = exp(apact)/sumpact
vpact = exp(vpact)/sumpact
mpact = exp(mpact)/sumpact

plot(arposp,apact)
plot(vposp, vpact)
plot(mposp,mpact)


####reconstruction activities####

#auditory reconstruction left
lrho = (t(wl) %*% apact) + (t(wlm) %*% mpact)

#auditory reconstruction right
rrho = (t(wr) %*% apact) + (t(wrm) %*% mpact)

plot(rposr,lrho)
lines(rposr,rrho)

arho=lrho*rrho
plot(rposr,arho)
Sa = rposr[which.max(arho)]
Sa

#visual reconstruction
vrho = (t(wv) %*% vpact) + (t(wvm) %*% mpact)
plot(vposr,vrho)
Sv = vposr[which.max(vrho)]
Sv

# calculate reconstruction error of auditory left and right input units
normarinput = arinput/max(arinput)
normalinput = alinput/max(alinput)
normlrho = lrho/max(lrho)
normrrho = rrho/max(rrho)

lerror = normlrho - normalinput
rerror = normrrho - normarinput

plot(arpos, lerror)
plot(arpos, rerror)

####readjust adaptation weights####
eta = 0.65
tau = 0.009
#if there is a stimulus at current time, the adaptation weights are updated according to the reconstruction error
ala = ala + eta * normalinput * lerror
ara = ara + eta * normarinput * rerror

#if there is no stimulus at current time, the adaption weight decays back to 1.
deltal = 1-ala #determine how far from 1 the adaptation weight currently is, for leftward units
deltar = 1-ara #determine how far from 1 the adaptation weight currently is, for rightward units

ala = ala + sign(deltal)*tau
ara = ara + sign(deltar)*tau