setwd("/Users/rachel/documents/r/brazilqmra")
library(fitdistrplus)
library(mc2d)
concdatsalt1<-read.csv("./concdatsalt1.csv",header=TRUE,sep=",")
marinagloria1<-concdatsalt1[concdatsalt1$X=="GN0034",]

#Scenario:Sailing in saltwater, using doses only from marina gloria 
expon<-function(a,k,d){
  return((a+(1-a))*(1-(exp(1)^(-k*d))))
}

bpois<-function(n50,alpha,d){
  return(1-(1+d*(2^(1/alpha)-1)/n50)^-alpha)
}

#according to original exponential dose response equation fit at drexel workshop, k=0.001715176

pathogensail<-data.frame(
  org     =c("ecoli",   "campylobacter", "salmonella", "rotavirus", "cryptosporidium", "ascaris", "entero"),
  alpha   =c(NA,        0.145,           0.3126,       0.2531,      NA,                NA,           NA),
  a       =c(0,           0,               0,            0,          0,                 0,       0.0015),
  n50     =c(NA,        896,             23600,        6.17,        NA,                NA,           NA),
  pdi     =c(1,        0.3,             0.3,          0.5,         0.7,              0.39,            1),
  k       =c(0.0001432, NA,              NA,           NA,     0.00419,             0.039,      0.00018),
  ratio   =c(1,         10^5,            10^5,         10^5,      10^6,            10^6,      0.00018),
  equation=c("expon",   "bpois",         "bpois",      "bpois",     "expon",      "expon",    "expon")
)

#Simulate infection risk and illness risk for each pathogen
simulator<-function(rowpath){
  ratio<-as.numeric(rowpath["ratio"])
  #when I put conc 1<-marinagloria$mpn.100.ml, turns the column in the dataframe into a factor. 
  conc1<-c(10,133,10,2729,79,1782,86,41,15,202,10)
  
  
  r<-mcstoc(rlnorm,meanlog=log(ratio),sdlog=1.4)
  dd<<-mcstoc(rlnorm,meanlog=mean(log(conc1)),sdlog=sd(log(conc1)))
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
