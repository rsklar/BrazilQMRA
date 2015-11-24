library(fitdistrplus)
library(mc2d)
concdatfresh<-read.csv("./concdatfresh.csv",header=TRUE,sep=",")
#swim 3 has a coefficient, swim 4 file doesn't use a for dose response calc. 

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
  pdi     =c(NA,        0.3,             0.3,          0.5,         0.7,               0.39),
  k       =c(0.0000511, NA,              NA,           NA,          0.00419,           0.0199),
  ratio   =c(1,         10^5,            10^5,         10^5,        10^6,              10^6),
  equation=c("expon",   "bpois",         "bpois",      "bpois",     "expon",           "expon")
)

#Convert FC to pathogen concentration
convertpath_fresh<-function(ratio){
  (concdatfresh$mpn.ml*0.843)/ratio
}

#For each pathogen go through the simulation
simulator<-function(rowpath){
  ratio<-as.numeric(rowpath["ratio"])
  conc<-convertpath_fresh(ratio)
  c<-mcstoc(rlnorm,type="V",meanlog=mean(log(conc)),sdlog=sd(log(conc)))
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
  print(rowpath["org"])
  swiminf<-mc(c,d,riski)
  print(swiminf)
  
  #risk disease per event 
  print(rowpath["org"])
  swimdis<-mc(c,d,riskd)
  print(swimdis)
  
  
}

(apply(pathogen,1,simulator))



