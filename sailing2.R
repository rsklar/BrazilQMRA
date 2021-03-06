library(fitdistrplus)
library(mc2d)
concdatsalt<-read.csv("./concdatsalt.csv",header=TRUE,sep=",")

#Scenario:Sailing in saltwater, using dose response for salt water fitted by Alex and Rachel from US data
expon<-function(k,d){
  return((1-(exp(1)^(-k*d))))
}

bpois<-function(n50,alpha,d){
  return(1-(1+d*(2^(1/alpha)-1)/n50)^-alpha)
}

#according to original exponential dose response equation fit at drexel workshop, k=0.001715176

pathogen<-data.frame(
  org=c("ecoli","campylobacter","salmonella","rotavirus","cryptosporidium", "ascaris"),
  alpha=c(NA, 0.145, 0.3126, 0.2531, NA, NA),
  n50=c(NA,896, 23600, 6.17, NA, NA),
  pdi=c(NA,0.3, 0.3, 0.5, 0.7, 0.39),
  k=c(0.000021502,NA, NA, NA, 0.00419, 0.0199),
  ratio=c(1, 10^5, 10^5, 10^5, 10^6, 10^6),
  equation=c("expon", "bpois", "bpois", "bpois", "expon", "expon")
)

#Convert FC to pathogen concentration
  #multiply by .01 to change mpn/100 to mpn/ml
  #0.175 ratio of FC:Ecoli in Freshwater from EPA 

concdatsalt1<-c(concdatsalt$Janeiro,concdatsalt$Fevereiro,concdatsalt$Marco,concdatsalt$Abril,concdatsalt$Maio,concdatsalt$Junho,concdatsalt$Julho, concdatsalt$Agosto,concdatsalt$Setembro, concdatsalt$Outubro)
convertpath_salt<-function(ratio) {
(concdatsalt1*0.01*0.175)/ratio  
}

#Simulate infection risk and illness risk for each pathogen
simulator<-function(rowpath){
  ratio<-as.numeric(rowpath["ratio"])
  conc<-convertpath_salt(ratio)
  c<-mcstoc(rlnorm,type="V",meanlog=mean(log(conc)),sdlog=sd(log(conc)))
  #i<-3.9 ml/hr- kayaking. A conservative estimate. Highest value of all recreational boat activities measured.
  #t<-0.52 hr- mean of average target sailing times reported by olympic trial committee
  d<-c*3.9*0.52
  
  #MC
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
  sailinf<-mc(c,d,riski)
  print(rowpath["org"])
  print(sailinf)
  summary(sailinf)
  
  #risk disease per event 
  saildis<-mc(c,d,riskd)
  print(rowpath["org"])
  print(saildis)
  summary(saildis)
}

apply(pathogen,1,simulator)




