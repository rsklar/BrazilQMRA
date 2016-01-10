setwd("/Users/rachel/documents/r/brazilqmra")
library(fitdistrplus)
library(mc2d)
concdatfresh<-read.csv("./concdatfresh.csv",header=TRUE,sep=",")

ndvar(10001)

#Scenario: Swimming in freshwater:
expon<-function(a,k,d){
  return((a+(1-a))*(1-(exp(1)^(-k*d))))
}
  
bpois<-function(n50,alpha,d){
  return(1-(1+d*(2^(1/alpha)-1)/n50)^-alpha)
}

pathogen<-data.frame(
  org     =c("ecoli",   "campylobacter", "salmonella", "rotavirus", "cryptosporidium", "ascaris"),
  alpha   =c(NA,        0.145,           0.3126,       0.2531,      NA,                NA),
  a       =c(0.000105,  0,               0,            0,           0,                 0),
  n50     =c(NA,        896,             23600,        6.17,        NA,                NA),
  pdi     =c(1,        0.3,             0.3,          0.5,         0.7,               0.39),
  k       =c(0.0000511, NA,              NA,           NA,          0.00419,           1),
  ratio   =c(1,         10^5,            10^5,         10^5,        10^6,              10^6),
  equation=c("expon",   "bpois",         "bpois",      "bpois",     "expon",           "expon")
)

#For each pathogen go through the Monte Carlo Simulations

simulator<-function(rowpath){
  ratio<-as.numeric(rowpath["ratio"])
  
  #We decided concentration will be lognormally distributed
    #i also tried the empirical distribution
  #Schets 2010 reports gamma distribution of ingestion and lorm distribution of duration(time)
  #information on r distribution from howard 2007  
  #?dd<-mcstoc(rempiricalD,values=concdatfresh$mpn.ml)
  
  #r=ratio,dd=distribution environmental data,c=concentration, i=ingestion volume, t=time
  r<<-mcstoc(func=rlnorm,meanlog=log(ratio),sdlog=1.4)
  dd<-mcstoc(func=rlnorm,meanlog= mean(log(concdatfresh$mpn.ml)),sdlog=sd((log(concdatfresh$mpn.ml))))
  c<-(dd*0.843)/r
  
  i<-mcstoc(func=rgamma,rate=0.45,shape=60)
  t<-mcstoc(func=rlnorm,meanlog=3.6,sdlog=0.85)
  
  #dose is a combo of distributions
  #c*i/t?????
  d<-(c*i)/t
  plot(c)
  plot(i)
  plot(t)
  
  riski<-if(rowpath["equation"]=="expon"){
    a<-as.numeric(rowpath["a"])
    k<-as.numeric(rowpath["k"])
    expon(a,k,d)
  } else if(rowpath["equation"]=="bpois"){
    n50<-as.numeric(rowpath["n50"])
    alpha<-as.numeric(rowpath["alpha"])
    bpois(n50,alpha,d)
  }
  plot(riski)
  riskd <- riski*as.numeric(rowpath["pdi"])

  #Risk Infection Per Event
  if(rowpath["org"] != "ecoli"){
    print(rowpath["org"])
    swiminf<<-mc(c,i,t,d,r,dd,riski)
    print(swiminf) 
  }
  
  #risk disease per event 
  print(rowpath["org"])
  swimdis<-mc(c,i,t,d,r,dd,riskd)
  print(swimdis)
  
}

(apply(pathogen,1,simulator))




