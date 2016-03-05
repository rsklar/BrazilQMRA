setwd("/Users/rachel/documents/r/brazilqmra")
library(fitdistrplus)
library(mc2d)
concdatsalt2<-read.csv("./concdatsalt2.csv",header=TRUE,sep=",")
marinagloria2<-concdatsalt2[concdatsalt2$X=="GN0034",]

#Scenario:Sailing in saltwater, using doses only from marina gloria 
expon<-function(a,k,d){
  return((a+(1-a))*(1-(exp(1)^(-k*d))))
}

bpois<-function(n50,alpha,d){
  return(1-(1+d*(2^(1/alpha)-1)/n50)^-alpha)
}

#according to original exponential dose response equation fit at drexel workshop, k=0.001715176

pathogensail<-data.frame(
  org     =c("ecoli", "enterococci",   "campylobacter", "salmonella", "rotavirus", "cryptosporidium", "ascaris"),
  alpha   =c(NA,          NA,              0.145,           0.3126,       0.2531,        NA,             NA),
  a       =c(0,          0.0015,              0,               0,            0,           0,              0),
  n50     =c(NA,            NA,               896,             23600,        6.17,        NA,            NA),
  pdi     =c(1,              1,               0.3,             0.3,          0.5,         0.7,         0.39),
  k       =c(0.0001432, 0.00018,              NA,              NA,           NA,          0.00419,    0.039),
  ratio   =c(1,             1,                10^5,            10^5,         10^5,        10^6,        10^6),
  equation=c("expon",   "expon",          "bpois",         "bpois",        "bpois",     "expon",    "expon")
)

#Simulate infection risk and illness risk for each pathogen
simulator<-function(rowpath){
  ratio<-as.numeric(rowpath["ratio"])
  conc<-marinagloria2$mpn.100.ml
  
  #r=ratio,dd=distribution environmental data,c=concentration, i=ingestion volume, t=time
  r<-mcstoc(rlnorm,meanlog=log(ratio),sdlog=1.4)
  dd<<-mcstoc(rlnorm,meanlog=mean(log(conc)),sdlog=sd(log(conc)))
  c<-(dd*0.01*0.843)/r
  
  i<<-mcstoc(runif,min=0,max=11.8)
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
}
  
  return(mc(c,i,d,r,dd,riskd))
  
}

saildis<-apply(pathogensail,1,simulator)
names(saildis)<-pathogensail$org

risktsail<-function(){
  1-(
    (1-sample(saildis[[1]]$riskd,size=1))*
    (1-sample(saildis[[2]]$riskd,size=1))*
    (1-sample(saildis[[3]]$riskd,size=1))*
    (1-sample(saildis[[4]]$riskd,size=1))*
    (1-sample(saildis[[5]]$riskd,size=1))*
    (1-sample(saildis[[6]]$riskd,size=1))
  )
}


totalrisksail<-replicate(n=10000,expr=risktsail(),simplify = "vector")

boxplot(totalrisksail,log="y")
