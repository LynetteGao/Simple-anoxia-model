rm(list=ls())

#example for Ksed input:
fked = c()
for(i in 1:150){
  fked[i] = 0.0002  #per day
}
for(i in 151:300){
  fked[i] = 0.01
}

area  = 3961.2*10000

#example for NEP input:
nep = c()
nep[1:110] = 0
nep[111:290] = 3.5*area/2

#hypo&epil volume over a year in 1980
volume<-read.table(file = "simGLM_ME_volumetric.txt",header = T)
sample_v<-volume[364:729,]
temp<-sample_v[2:3]
vol<-sample_v[6:7]

#only trace epilimninon in those mixing day
# function to determine the day to stratify
#in the example data :day 110
day_stratify<-function(temp_data){
  for(day in 1:365){
    if(temp_data[day,1]-temp_data[day,2]>1){
      break;
    }
  }
  return (day)
}

# function to determine the day to mix
#in the example data : day 290
day_mix<-function(temp_data){
  for(day in 200:365){
    if(temp_data[day,1]>temp_data[day,2]&&temp_data[day,1]-temp_data[day,2]<1){
      break;
    }
  }
  return (day)
}

mix_day<-day_mix(temp)
stratify_day<-day_stratify(temp)


#starting temperature at day one, and O2 saturation at day one
temp = 5

DO_simulate<-function(NEP,Ksed){
  
 
  #parameters for epilimnion
  K600 = k.cole.base(2)
  In = 0
  Out = 0
  #parameters for hypo
  Fhypo=c()
  hypo.02.at.sat.base = o2.at.sat.base(temp = sample_v$avgHypoT,altitude = 300)
  hypo_o2 = c()
  epil.O2 = c()
  kO2 = c()
  
  Fatm=c()
  o2sat=c()
  
  
  o2satat5 = o2.at.sat.base(temp=temp,altitude = 300)
  

  epil.O2[1]=o2satat5*vol[1,1]
  hypo_o2[1] = epil.O2[1]
  
  #sample o2 for > mixing day 
  epi.02.at.sat.base = o2.at.sat.base(temp = sample_v$avgEpiT,altitude = 300)
  
  #simulation for hypolimnion and epilimnion DO
  for(i in 2:365){
    #update the temperature
    if(i<210){
      temp = temp+0.1
    }else{
      temp = temp - 0.13
    }
    #Day: NOT stratified yet(No exchange between epil&hypo)
    if(i<=stratify_day){
     
      #epilimnion:In-out+nep+Fatm
      kO2[i] <- k600.2.kGAS.base(k600=K600,temperature=temp,gas='O2')
      o2sat[i]<-o2.at.sat.base(temp=temp,altitude = 300)
      Fatm[i] <- kO2[i]*(o2sat[i]-epil.O2[i-1]/vol[i-1,1])*area
      epil.O2[i] = epil.O2[i-1]+In-Out+Fatm[i]+NEP[i]
      #Hypo(set it the same as epil)
      hypo_o2[i] = epil.O2[i]
    }
    #Day 111 - 290
    else if(stratify_day<i&i<mix_day){
      #calculate flux(unit?)
      volumechange_hypo = vol[i,2]-vol[i-1,2]  #in m^3?
      volumechange_hypo_proportion =  volumechange_hypo/vol[i-1,2]
      o2delivery = volumechange_hypo_proportion*hypo_o2[i-1]
      #O2 delivery/ Volume of hypo(?) = flux
      Fhypo[i] = o2delivery
      
      #Hypo: Fhypo- Fsed
      k = Ksed[i]
      Fsed = k*hypo_o2[i-1]
      hypo_o2[i] = hypo_o2[i-1]+Fhypo[i]-Fsed
      #epil
      kO2[i] <- k600.2.kGAS.base(k600=K600,temperature=temp,gas='O2')#m/day
      o2sat[i]<-o2.at.sat.base(temp=temp,altitude = 300)
      Fatm[i] <- kO2[i]*(o2sat[i]-epil.O2[i-1]/vol[i-1,1])*area 
      epil.O2[i] = epil.O2[i-1]+Fatm[i]+NEP[i]-Fhypo[i]
    }
    #day = mix day
    else {
      epil.O2[i] = (hypo.02.at.sat.base[i]*sample_v$HypoV[i]+epi.02.at.sat.base[i]*sample_v$EpiV[i])/(sample_v$HypoV[i]+sample_v$EpiV[i])
    }
    # else if(i==mix_day){
    #   epil.O2[i] = (hypo.02.at.sat.base[i]*sample_v$HypoV[i]+epi.02.at.sat.base[i]*sample_v$EpiV[i])/(sample_v$HypoV[i]+sample_v$EpiV[i])
    #   
    #   #epil.O2[i] = (epil.O2[i-1]*vol[i,1]+hypo_o2[i-1]*vol[i,2])/(vol[i,1]+vol[i,2])
    #   hypo_o2[i] = epil.O2[i]
    # }
    # else if(i>mix_day){
    #   #epilimnion:In-out+nep+Fatm
    #   kO2<- k600.2.kGAS.base(k600=K600,temperature=temp,gas='O2')
    #   o2sat<-o2.at.sat.base(temp=temp,altitude = 300)
    #   Fatm <- kO2*(o2sat-epil.O2[i-1]) 
    #   epil.O2[i] = epil.O2[i-1]+In-Out+Fatm+NEP[i]
    #   #Hypo(set it the same as epil)
    #   hypo_o2[i] = epil.O2[i]
    # }
  # }
    
  }  
  return(list(hypo = hypo_o2,epil=epil.O2,ko2 = kO2,fatm = Fatm))
}
result<-DO_simulate(nep,fked)

x <- seq(1,365)#plot every 5 days
#plot(x,y = result$ko2,pch =1,ylim = c(0,50),main="DO simulation",xlab = "day",ylab = "DO")
plot(x,y = result$hypo[x]/vol[x,2],pch =1,ylim = c(0,50),main="DO simulation",xlab = "day",ylab = "DO")
points(x,y=result$epil[x]/vol[x,1],col="blue")



legend("topright",
       c("epilimnion(model)","hypolimnion(model)","epilimnion(1980)","hypolimnion(1980)"),
       fill=c("blue","black","yellow","purple"),cex = 0.75)

#true sampling data
points(sample_v$avgHypoDO*0.032,col="purple")
points(o2.at.sat.base(temp = sample_v$avgEpiT,altitude = 300),col="yellow")


