# clear workspace
rm(list = ls(all = T))


# load packages
library(ez)
library(plyr)
library(ggplot2)

##auditory gain
ga = 140
#auditory rise
m = 20
#visual gain
gV = 80
#visual tuning width
sig = 20 #80

L = 301

####adaptation parameters####
eta = 0.65
tau = 0.009

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


#create stimulus array
#use same stimulus sequence as in Bosen et al experiment
stimuli = c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1) #is there a stimulus during this trial?
visual = c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) #is there a visual stimulus with the auditory?
aestimates = rep(NA,length(stimuli)) #create a vector for auditory position estimates
vestimates = rep(NA,length(stimuli)) #create a vector for visual position estimates
ala = rep(1,301) #initialize adaptation parameter left
ara = rep(1,301) #initialize adaptation parameter right

iteration = 1
for (i in stimuli){
  if (i == 1){#stimulus present, first determine if it is a multisensory trial
    gv = ifelse(visual[iteration] == 1, gV, 0)
    #determine unit activities, estimates and reconstruction errors. finally readjust adaptation weights
    ####DETERMINE UNIT ACTIVITIES####
    
    #STIMULUS LOCATIONS#
    xa =0  #auditory location
    xv = 8  #visual location
    
    #initial multisensory input#
    
    mult = 10.7
    
    ####input activity####
    #determine input activity given auditory position
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
    
    #activation non-linear exponential
    sumpact = 1+((sum(exp(apact)) + sum(exp(vpact)) + sum(exp(mpact))) / (length(apact) + length(vpact) + length(mpact)))
    
    apact = exp(apact)/sumpact
    vpact = exp(vpact)/sumpact
    mpact = exp(mpact)/sumpact
    
    ####reconstruction activities####
    
    #auditory reconstruction left
    lrho = (t(wl) %*% apact) + (t(wlm) %*% mpact)
    
    #auditory reconstruction right
    rrho = (t(wr) %*% apact) + (t(wrm) %*% mpact)
    
    #determine auditory estimate
    arho=lrho*rrho
    Sa = rposr[which.max(arho)]
    aestimates[iteration] = Sa
    
    #visual reconstruction
    vrho = (t(wv) %*% vpact) + (t(wvm) %*% mpact)
    
    #determine visual estimate
    Sv = vposr[which.max(vrho)]
    vestimates[iteration] = Sv
    
    # calculate reconstruction error of auditory left and right input units
    normarinput = arinput/max(arinput)
    normalinput = alinput/max(alinput)
    normlrho = lrho/max(lrho)
    normrrho = rrho/max(rrho)
    
    lerror = normlrho - normalinput
    rerror = normrrho - normarinput
    
    #readjust adaptation weights
    ala = ala + eta * normalinput * lerror
    ara = ara + eta * normarinput * rerror
    }
  else {#no stimulus decay adaptation weights towards 1
    aestimates[iteration] = NA
    vestimates[iteration] = NA
    
    deltal = 1-ala #determine how far from 1 the adaptation weight currently is, for leftward units
    deltar = 1-ara #determine how far from 1 the adaptation weight currently is, for rightward units
    
    ala = ala + sign(deltal)*tau
    ara = ara + sign(deltar)*tau
  }
  iteration = iteration + 1
}

trial = 1:length(stimuli)
plot(trial,aestimates)