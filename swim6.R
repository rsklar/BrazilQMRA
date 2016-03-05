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
  #create an empty vector to store riskd
  riskvector<-vector(mode="numeric", length=0)

simulator<-function(rowpath){
  ratio<-as.numeric(rowpath["ratio"])
  
  #r=ratio,dd=distribution environmental data,c=concentration, i=ingestion volume, t=time
  r<<-mcstoc(func=rlnorm,meanlog=log(ratio),sdlog=1.4)
  #dd=pathogens per mL
  dd<-mcstoc(func=rlnorm,meanlog= mean(log(concdatfresh$mpn.ml)),sdlog=sd((log(concdatfresh$mpn.ml))))
  c<-(dd*0.843)/r
  
  i<-mcstoc(func=rgamma,rate=0.45,shape=60)

  #dose is a combo of distributions

  d<-(c*i)
  
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
    swiminf<<-mc(c,i,d,r,dd,riski)
  }
  
#swim disease is contained in this return value
  
  return(mc(c,i,d,r,dd,riskd))
}
swimdis<-apply(pathogen,1,simulator)
names(swimdis)<-pathogen$org


#function to generate mean(riskd) for each pathogen

pathriskd<-function(swimmc){
  mean(swimmc$riskd)
}

meanrisk<-lapply(swimdis,pathriskd)
names(meanrisk)<-pathogen$org

#riskdiarrhea<-1-((1-meanrisk$ecoli)*(1-meanrisk$campylobacter)*(1-meanrisk$salmonella)*(1-meanrisk$rotavirus)*(1-meanrisk$cryptosporidium)*(1-meanrisk$ascaris))
#names(riskdiarrhea)<-"Swimming Diarrhea Risk"
#print(riskdiarrhea)

riskdiarrhea<-1-((1-sample(swimdis[[1]]$riskd,size=1))*(1-sample(swimdis[[2]]$riskd,size=1))*(1-sample(swimdis[[3]]$riskd,size=1))*(1-sample(swimdis[[4]]$riskd,size=1))*(1-sample(swimdis[[5]]$riskd,size=1))*(1-sample(swimdis[[6]]$riskd,size=1))
                 
)


