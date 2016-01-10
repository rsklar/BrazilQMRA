setwd("/Users/rachel/documents/r/brazilqmra")
library(fitdistrplus)
library(mc2d)
concdatsalt<-read.csv("./concdatsalt.csv",header=TRUE,sep=",")

#Scenario:Sailing in saltwater, using dose response for salt water fitted by Alex and Rachel from US data
#using a distribution for i, uniform distribution using UCL as max, 0 as min
expon<-function(a,k,d){
  return((a+(1-a))*(1-(exp(1)^(-k*d))))
}

bpois<-function(n50,alpha,d){
  return(1-(1+d*(2^(1/alpha)-1)/n50)^-alpha)
}

#according to original exponential dose response equation fit at drexel workshop, k=0.001715176

pathogen<-data.frame(
  org     =c("ecoli",   "campylobacter", "salmonella", "rotavirus", "cryptosporidium", "ascaris"),
  alpha   =c(NA,        0.145,           0.3126,       0.2531,      NA,                NA),
  a       =c(0.000105,  0,               0,            0,           0,                 0),
  n50     =c(NA,        896,             23600,        6.17,        NA,                NA),
  pdi     =c(1,        0.3,             0.3,          0.5,         0.7,               0.39),
  k       =c(0.000143199657478504, NA,              NA,           NA,          0.00419,           1),
  ratio   =c(1,         10^5,            10^5,         10^5,        10^6,              10^6),
  equation=c("expon",   "bpois",         "bpois",      "bpois",     "expon",           "expon")
)

#Simulate infection risk and illness risk for each pathogen
simulator<-function(rowpath){
  ratio<-as.numeric(rowpath["ratio"])
  conc<<-concdatsalt$mpn.100.ml
  

  r<-mcstoc(rlnorm,meanlog=log(ratio),sdlog=1.4)
  dd<-mcstoc(rlnorm,meanlog=mean(log(conc)),sdlog=sd(log(conc)))
  c<-(dd*0.01*0.63)/r
  
  i<-mcstoc(runif,min=0,max=11.8)
  #t<-0.52 hr- mean of average target sailing times reported by olympic trial committee
  d<-c*i*0.52
  
  #MC
  ndvar(10001)
  
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
  if(rowpath["org"] != "ecoli"){
  sailinf<-mc(c,dd,d,r,riski)
  print(rowpath["org"])
  print(sailinf)
  summary(sailinf)
}
  
  #risk disease per event 
  saildis<-mc(c,dd,d,r,riskd)
  print(rowpath["org"])
  print(saildis)
  summary(saildis)
  plot(saildis)
}

apply(pathogen,1,simulator)




