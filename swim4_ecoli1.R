library(fitdistrplus)
library(mc2d)
concdatfresh<-read.csv("./concdatfresh.csv",header=TRUE,sep=",")
ndvar(10001)

#swim 4a uses the a neha suggests in her dose response for ecoli.  

#Scenario: Swimming in freshwater:
expon<-function(a,k,d){
  return(a+(1-a)*(1-(exp(1)^(-k*d))))
}
  
bpois<-function(n50,alpha,d){
  return(1-(1+d*(2^(1/alpha)-1)/n50)^-alpha)
}

pathogen<-data.frame(
  org     =c("ecoli",   "campylobacter", "salmonella", "rotavirus", "cryptosporidium", "ascaris"),
  alpha   =c(0.1778,        0.145,           0.3126,       0.2531,      NA,                NA),
  a       =c(1,  1,               1,            1,           1,                 1),
  n50     =c((8.6*10^7),        896,             23600,        6.17,        NA,                NA),
  pdi     =c(1,        0.3,             0.3,          0.5,         0.7,               0.39),
  k       =c(0.0000511, NA,              NA,           NA,          0.00419,           0.0199),
  ratio   =c(1,         10^5,            10^5,         10^5,        10^6,              10^6),
  equation=c("bpois",   "bpois",         "bpois",      "bpois",     "expon",           "expon")
)

#Convert FC to pathogen concentration
convertpath_fresh<-function(ratio){
  (concdatfresh$mpn.ml*0.843)/ratio
}

#For each pathogen go through the Monte Carlo Simulations

simulator<-function(rowpath){
  ratio<-as.numeric(rowpath["ratio"])
  conc<-convertpath_fresh(ratio)
  #We decided concentration will be lognormally distributed
  #Schets 2010 reports gamma distribution of ingestion and lorm distribution of duration(time)
  c<-mcstoc(func=rlnorm,type="V",meanlog=mean(log(conc)),sdlog=sd(log(conc)))
  i<-mcstoc(func=rgamma,type="V",rate=0.45,shape=60)
  t<-mcstoc(func=rlnorm,type="VU",meanlog=3.6,sdlog=0.85)
  d<-c*i*t
  
  riski<<-if(rowpath["equation"]=="expon"){
    a<-as.numeric(rowpath["a"])
    k<-as.numeric(rowpath["k"])
    expon(a,k,d)
  } else if(rowpath["equation"]=="bpois"){
    n50<-as.numeric(rowpath["n50"])
    alpha<-as.numeric(rowpath["alpha"])
    bpois(n50,alpha,d)
  }
  
  riskd <- riski *as.numeric(rowpath["pdi"])

  #risk infection per event
  if(rowpath["org"] != "ecoli"){
    print(rowpath["org"])
    swiminf<<-mc(c,i,t,d,riski)
    #print(swiminf) 
    print(summary(swiminf))

  }

  #risk disease per event 
  print(rowpath["org"])
  swimdis<-mc(c,i,d,riskd)
  #print(swimdis)
  print(summary(swimdis))
  plot(swimdis, na.rm=TRUE)
  
  
}

(apply(pathogen,1,simulator))




