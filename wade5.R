setwd("/Users/rachel/documents/r/brazilqmra")
library(fitdistrplus)
library(mc2d)
concdatfresh<-read.csv("./concdatfresh.csv",header=TRUE,sep=",")

#Scenario: Swimming in freshwater:
expon<-function(k,d){
  return((1-(exp(1)^(-k*d))))
}

bpois<-function(n50,alpha,d){
  return(1-(1+d*(2^(1/alpha)-1)/n50)^-alpha)
}

pathogen<-data.frame(
  org     =c("ecoli",   "campylobacter", "salmonella", "rotavirus", "cryptosporidium", "ascaris"),
  alpha   =c(NA,        0.145,           0.3126,       0.2531,      NA,                NA),
  n50     =c(NA,        896,             23600,        6.17,        NA,                NA),
  pdi     =c(1,        0.3,             0.3,          0.5,         0.7,               0.39),
  k       =c(0.0000511, NA,              NA,           NA,          0.00419,           1),
  ratio   =c(1,         10^5,            10^5,         10^5,        10^6,              10^6),
  equation=c("expon",   "bpois",         "bpois",      "bpois",     "expon",           "expon")
)


#For each pathogen go through the simulation
simulator<-function(rowpath){
  ratio<-as.numeric(rowpath["ratio"])
  
  #r=ratio,dd=distribution environmental data,c=concentration
  #ask patrick-meanlogs of the original data and of the ratio. if i should run this as an mc stoch or not. 
  #check to see if you should divide or multiply by the ratio. 
  #should i make an mc stoch of "d<-mcstoc(func=rpois,type="V",lambda=dose)"
  #information on r distribution from howard 2007
  #empirical distribution for data. check on this
  
  r<<-mcstoc(func=rlnorm,meanlog=log(ratio),sdlog=1.4)
  
  dd<-mcstoc(func=rlnorm,meanlog= mean(log(concdatfresh$mpn.ml)),sdlog=sd((log(concdatfresh$mpn.ml))))
  c<-(dd*0.843)/r
  
  #i<-3.23-average ingestion per wading event from 3 papers 
  d<-c*3.23
  
  
  #Monte Carlo Simulations
  ndvar(10001)
  
  riski<-if(rowpath["equation"]=="expon"){
    k<-as.numeric(rowpath["k"])
    expon(k,d)
  } else if(rowpath["equation"]=="bpois"){
    n50<-as.numeric(rowpath["n50"])
    alpha<-as.numeric(rowpath["alpha"])
    bpois(n50,alpha,d)
  }
  
  riskd <- riski *as.numeric(rowpath["pdi"])
  
  #risk infection per event
  
  if(rowpath["org"] != "ecoli") {
  print(rowpath["org"])
  wadeinf<-mc(c,d,r,dd,riski)
  print(wadeinf)
}
  #risk disease per event 
  print(rowpath["org"])
  wadedis<-mc(c,d,r,dd,riskd)
  print(wadedis)
  
  
}

(apply(pathogen,1,simulator))



