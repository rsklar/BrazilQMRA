setwd("/Users/rachel/documents/r/brazilqmra")
library(fitdistrplus)
library(mc2d)
concdatfresh<-read.csv("./concdatfresh.csv",header=TRUE,sep=",")

ndvar(10001)

#Scenario: Wading in Freshwater
expon<-function(a,k,d){
  return((a+(1-a))*(1-(exp(1)^(-k*d))))
}

bpois<-function(n50,alpha,d){
  return(1-(1+d*(2^(1/alpha)-1)/n50)^-alpha)
}

pathogenwade<-data.frame(
  org     =c("ecoli",   "campylobacter", "salmonella", "rotavirus", "cryptosporidium", "ascaris"),
  alpha   =c(NA,        0.145,           0.3126,       0.2531,      NA,                NA),
  a       =c(0.000105,  0,               0,            0,           0,                 0),
  n50     =c(NA,        896,             23600,        6.17,        NA,                NA),
  pdi     =c(1,        0.3,             0.3,          0.5,         0.7,               0.39),
  k       =c(0.0000511, NA,              NA,           NA,          0.00419,           0.039),
  ratio   =c(1,         10^5,            10^5,         10^5,        10^6,              10^6),
  equation=c("expon",   "bpois",         "bpois",      "bpois",     "expon",           "expon")
)

#For each pathogen go through the Monte Carlo Simulations

simulator<-function(rowpath){
  ratio<-as.numeric(rowpath["ratio"])

  #r=ratio,dd=distribution environmental data,c=concentration, i=ingestion volume, t=time. no separate i or t for wading-using a point value for mean ingestion per event.
  r<<-mcstoc(func=rlnorm,meanlog=log(ratio),sdlog=1.4)
  #dd=pathogens per mL
  dd<<-mcstoc(func=rlnorm,meanlog= mean(log(concdatfresh$mpn.ml)),sdlog=sd((log(concdatfresh$mpn.ml))))
  c<-(dd*0.843)/r

  #i<-3.23-average ingestion per wading event from 3 papers 
  #exposure dose
  d<-c*3.23
  
  
  #Monte Carlo Simulations

  
  riski<-if(rowpath["equation"]=="expon"){
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
  
  if(rowpath["org"] != "ecoli") {
  
  wadeinf<-mc(c,d,r,dd,riski)
}
  return(mc(c,d,r,dd,riskd))
  
}

wdis<-(apply(pathogenwade,1,simulator))
names(wdis)<-pathogenwade$org

risktwade<-function(){
  1-(
    (1-sample(wdis[[1]]$riskd,size=1))*
      (1-sample(wdis[[2]]$riskd,size=1))*
      (1-sample(wdis[[3]]$riskd,size=1))*
      (1-sample(wdis[[4]]$riskd,size=1))*
      (1-sample(wdis[[5]]$riskd,size=1))*
      (1-sample(wdis[[6]]$riskd,size=1))
  )
}


totalriskwade<-replicate(n=10000,expr=risktwade(),simplify = "vector")

boxplot(totalriskwade,log="y")