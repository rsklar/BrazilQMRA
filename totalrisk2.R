#sailing
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
        
        #total risk includes pathogens, not indicator org.
        totalrisksail<-replicate(n=10000,expr=risktsail(),simplify = "vector")

#wading
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
              (1-sample(wdis[[2]]$riskd,size=1))*
              (1-sample(wdis[[3]]$riskd,size=1))*
              (1-sample(wdis[[4]]$riskd,size=1))*
              (1-sample(wdis[[5]]$riskd,size=1))*
              (1-sample(wdis[[6]]$riskd,size=1))
          )
        }
        
        
        #total risk includes pathogens, not indicator org.
        totalriskwade<-replicate(n=10000,expr=risktwade(),simplify = "vector")
        

#Swimming
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
        
        pathogenswim<-data.frame(
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
        swimdis<-apply(pathogenswim,1,simulator)
        names(swimdis)<-pathogenswim$org
        
        risktswim<-function(){
          1-(
              (1-sample(swimdis[[2]]$riskd,size=1))*
              (1-sample(swimdis[[3]]$riskd,size=1))*
              (1-sample(swimdis[[4]]$riskd,size=1))*
              (1-sample(swimdis[[5]]$riskd,size=1))*
              (1-sample(swimdis[[6]]$riskd,size=1))
          )
        }
        
        
        #total risk includes pathogens, not indicator org.
        totalriskswim<-replicate(n=10000,expr=risktswim(),simplify = "vector")
        
#Mean of distributions of each pathogen/exposure scenario (For total pathogen specific risks)
        #Swim
        pathriskswim<-c(
          mean(swimdis$campylobacter[["riskd"]]),
          mean(swimdis$salmonella[["riskd"]]),
          mean(swimdis$rotavirus[["riskd"]]),
          mean(swimdis$cryptosporidium[["riskd"]]),
          mean(swimdis$ascaris[["riskd"]])
        )
        names(pathriskswim)<-c("Campylobacter","Salmonella","Rotavirus","Cryptosporidium","Ascaris")
        #Wade   
        pathriskwade<-c(
          mean(wdis$campylobacter[["riskd"]]),
          mean(wdis$salmonella[["riskd"]]),
          mean(wdis$rotavirus[["riskd"]]),
          mean(wdis$cryptosporidium[["riskd"]]),
          mean(wdis$ascaris[["riskd"]])
        )
        
        names(pathriskwade)<-c("Campylobacter","Salmonella","Rotavirus","Cryptosporidium","Ascaris")
        
        #Sail        
        pathrisksail<-c(
          mean(saildis$campylobacter[["riskd"]]),
          mean(saildis$salmonella[["riskd"]]),
          mean(saildis$rotavirus[["riskd"]]),
          mean(saildis$cryptosporidium[["riskd"]]),
          mean(saildis$ascaris[["riskd"]])
        )
        
        names(pathrisksail)<-c("Campylobacter","Salmonella","Rotavirus","Cryptosporidium","Ascaris")
        
#barplot of total risks in each scenario
        barplot(c(mean(totalriskswim),
                  mean(totalriskwade),
                  mean(totalrisksail)),
                col=c("#0000CD","#0000CD","#00BFFF"),
                xlab="Exposure Scenario",
                names=c("Swimming","Wading","Sailing"),
                ylab="Mean Diarrhea Risk/Exposure Event",
                ylim=c(0,1),main="Diarrhea Risk"
        )
        
        legend("topright",c("Freshwater","Saltwater"),fill=c("#0000CD","#00BFFF"))  
  
#barplot of pathogen specific diarrhea risks
        df<-data.frame(pathriskswim,pathriskwade,pathrisksail)  
        
        barplot(as.matrix(df),xlab="Exposure Scenario",ylab= "Risk/Exposure Scenario (log scale)",
                beside=TRUE,log="y",
                main="Estimated Pathogen Specific Diarrhea Risk",
                font.main=10,
                #col=c("#FFC65D", "#7BC8A4", "#4CC3D9", "#93648D", "#404040","#404040"),
                ylim=c(1e-09,10),
                names=c("Swimming","Wading","Sailing")
        )
        
        #legend("topright",legend=pathogen$org,fill=c("#FFC65D", "#7BC8A4", "#4CC3D9", "#93648D", "#404040","#404040"))
        
#indicator ecoli/ENT versus total risk swim and sail 
        #SWIMMING
        barplot(c(mean(swimdis$ecoli[["riskd"]]),mean(totalriskswim)),
                main="Indicator Versus Pathogen Risk-Swim",
                names=c("E.coli","Criteria Pathogen")
                )
        #WADING
        barplot(c(mean(wdis$ecoli[["riskd"]]),mean(totalriskwade)),
                main="Indicator Versus Pathogen Risk-Wade",
                names=c("E.coli","Criteria Pathogen")
        )
        #Sailing
        barplot(c(mean(saildis$entero[["riskd"]]),mean(totalrisksail)),
                main="Indicator Versus Pathogen Risk-Sail",
                names=c("Enterococci","Criteria Pathogen")
        )
        
#NEED INDICATOR VERSUS PATHOGEN RISK SAIL using ENT 
        #Charts we want:
        
        #1.	estimated pathogen specific diarrhea-File: total  risk1
        #2.	indicator ecoli versus risk predicted from these pathogens
        #3.	indicator enterococci versus risk predicted from these pathogens 
      