#River Vs Sewage 
      #data frame of concentration per mL in river 
      swim_sewage_compare<-data.frame(
        org=c("ecoli","campylobacter","salmonella","rotavirus","cryptosporidium", "ascaris"),
        orgsperml=c(9.50E+09, 1.21E+05, 1.34E+05, 3.62E+04, 1.36E+04, 3.19E+04),
        #sewage min and max per 100 mL for each org.
        conc_min=c(10^6, 10^4, 0.2, 400, 0.1, 0.5),
        conc_max=c(10^7, 10^5, 8000, 85000, 39, 11)
      )
      
      #dataframe of river concentrations(orgs/ml) versus average raw sewage (orgs/100 mL)
      minmax<-cbind(swim_sewage_compare$conc_min,swim_sewage_compare$conc_max)
      meansewageconc<-apply(minmax,1,mean)
      meansewageconcperml<-meansewageconc/100
      rivrVsewage<-cbind(swim_sewage_compare,meansewageconcperml)
      riverVsewagetable<-cbind(riverconc=rivrVsewage$orgsperml,sewageconc=rivrVsewage$meansewageconcperml)
      row.names(riverVsewagetable)<-swim_sewage_compare$org
      print(riverVsewagetable)
      
      #make a bargraph of river concentrations versus raw sewage
        #transpose river vs. sewage table
        riverVsewage<-t(riverVsewagetable)
        print(riverVsewage)
        chartdata<-subset(riverVsewage,select=-c(ecoli))
        
        
        #graph river vs. sewage with ecoli 
          barplot(riverVsewage,xlab="Pathogens",col=c("darkgrey","black"),log="y",ylab="Organism/Ml (log scale)",main="River versus Sewage Concentrations",sub="including E.coli",beside=TRUE)
          legend("topright",c("River Concentration","Sewage Concentration (WHO 2003)"),fill=c("darkgrey","black")) 
          
        #graph river vs sewage without ecoli counts 
          #barplot(as.matrix(chartdata),xlab="Pathogens", ylab= "Organism/ml (log scale)",
          #legend=rownames(riverVsewage),beside=TRUE,log="y",main="River versus Sewage Concentrations",font.main=10)
  
#BAY vs Sewage
        #dataframe of BAY concentrations(orgs/ml) versus average raw sewage (orgs/100 mL)
        bay_sewage_compare<-data.frame(
          org=c("ecoli","campylobacter","salmonella","rotavirus","cryptosporidium", "ascaris"),
          orgsperml=c(2.53E+00, 2.85E-05, 2.84E-05, 2.68E-05, 2.89E-06, 2.58E-06),
          
          #sewage min and max per 100 mL for each org.
          conc_min=c(10^6, 10^4, 0.2, 400, 0.1, 0.5),
          conc_max=c(10^7, 10^5, 8000, 85000, 39, 11)
        )  
        
        #dataframe of bay concentrations(orgs/ml) versus average raw sewage (orgs/100 mL)
        minmax<-cbind(bay_sewage_compare$conc_min,swim_sewage_compare$conc_max)
        meansewageconc<-apply(minmax,1,mean)
        meansewageconcperml<-meansewageconc/100
        bayvsewage<-cbind(bay_sewage_compare,meansewageconcperml)
        bayvsewagetable<-cbind(bayconc=bayvsewage$orgsperml,sewageconc=bayvsewage$meansewageconcperml)
        row.names(bayvsewagetable)<-bay_sewage_compare$org
        print(bayvsewagetable)
        
        #make a bargraph of bay concentrations versus raw sewage
        #transpose bay vs. sewage table
        bayvsewage<-t(bayvsewagetable)
        print(bayvsewage)
        baychartdata<-subset(bayvsewage)
        #baychartdata<-subset(bayvsewage,select=-c(ecoli))
        
        #graph bay vs sewage without ecoli counts 
        
        barplot(as.matrix(baychartdata),xlab="Pathogens", ylab= "Organisms/ml (log scale)",
                beside=TRUE,log="y",main="Bay versus Sewage Concentrations",font.main=10,col=c("darkgrey","black"))
        
        legend("topright",c("Bay Concentration","Sewage Concentration (WHO 2003)"),fill=c("darkgrey","black")) 
        