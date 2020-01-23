rm(list=ls())
######### functions############
# function to determine the day to stratify
day_stratify<-function(temp_data){
  for(day in 1:365){
    if(temp_data[day,1]-temp_data[day,2]>1){
      break;
    }
  }
  return (day)
}

# function to determine the day to mix
day_mix<-function(temp_data){
  for(day in 200:365){
    if(temp_data[day,1]>temp_data[day,2]&&temp_data[day,1]-temp_data[day,2]<1){
      break;
    }
  }
  return (day)
}
###############################

######## data(temp&vol) in 1980 ######
data<-read.table(file = "simGLM_ME_volumetric.txt",header = T)
sample_v<-data[364:729,]
temp_data<-sample_v[2:3]
vol<-sample_v[6:7]
######################################

#######configuration for constant #########
ndays = 365
area  = 3961.2*10000
K600 = k.cole.base(2)
min_z_epil = 1
nep_constant = 1
In = 0
Out = 0
temp = 5 #starting temperature at day one
o2satat5 = o2.at.sat.base(temp=temp,altitude = 300)
mix_day<-day_mix(temp_data)
stratify_day<-day_stratify(temp_data)

##########################################


######### initialize variable for function input #############

nep_input = rep(0,ndays)    
#nep[1:110] = 0
#nep[111:290] = 1*area*2  #g/m^3 * m^2 *m = g/m *m =g

fked = rep(0,ndays)
for(i in 1:stratify_day){
  fked[i] = 0.0002  #per day
}
for(i in (stratify_day+1):ndays){
  fked[i] = 0.01
}

###########################################


##########simulation function #############
DO_simulate<-function(NEP,Ksed){
  
  ###initialize paramters###
  Fhypo=rep(0,ndays)
  Fatm=rep(0,ndays)
  Fsed = rep(0,ndays)
  z_epil = rep(0,365)
  
  kO2 = c()
  kO2[1] = k600.2.kGAS.base(k600=K600,temperature=5,gas='O2')
  
  o2sat=c()
  o2sat[1] = o2.at.sat.base(temp=5,altitude = 300)
  
  epil.O2 = c()
  epil.O2[1]=o2satat5*vol[1,1]
  
  hypo_o2 = c()
  hypo_o2[1] = epil.O2[1]/vol[1,1]*vol[1,2]
  
  ### simulation ###
  for(i in 2:ndays){
    #update the temperature
    if(i<210){
      temp = temp+0.1
    }else{
      temp = temp - 0.13
    }
    
    #(1)epilimnion:In-out+nep+Fatm
    kO2[i-1] <- k600.2.kGAS.base(k600=K600,temperature=temp,gas='O2')
    o2sat[i-1]<-o2.at.sat.base(temp=temp,altitude = 300)
    z_epil[i-1]<-vol[i-1,1]/area
    z_epil[i-1]<-max(min_z_epil,z_epil[i-1])
   
    Fatm[i-1] <- kO2[i-1]*(o2sat[i-1]*vol[i-1,1]-epil.O2[i-1])/z_epil[i-1]
    ###        m/d  * (g   -   g   )  * m^2/m^3 = m/d *g /m =g/d
    
    epil.O2[i] = epil.O2[i-1]+In-Out+Fatm[i-1]+NEP[i-1]

    
    #(2)hypolimnion (different depends on the day)
    if(i<stratify_day){
   
      #Hypo(set it the same as epil)
      hypo_o2[i] = epil.O2[i]/vol[i,1]*vol[i,2]
    }
    #Day 110 - 290
    else if(stratify_day<=i&i<mix_day){
      #calculate flux(unit?)
      volumechange_hypo = vol[i,2]-vol[i-1,2]  #in m^3
      volumechange_hypo_proportion =  volumechange_hypo/vol[i-1,2]
      Fhypo[i-1] =volumechange_hypo_proportion*hypo_o2[i-1]

      Fsed[i-1] = Ksed[i-1] *hypo_o2[i-1]

      #Hypo: Fhypo- Fsed
      hypo_o2[i] = hypo_o2[i-1]+Fhypo[i-1]-Fsed[i-1]
     
      
      NEP[i] = nep_constant * vol[i,1]  ##??
    }
    #day > mix day  ???
    else {
      epil.O2[i] = (epil.O2[i-1]+hypo_o2[i-1])/2
      hypo_o2[i] = epil.O2[i] 
      kO2[i]=0
      o2sat[i]=0
      NEP[i]=0  #??
    }
    
    
  }  
  #return(list(input = In,output = Out, nep = NEP, hypo = hypo_o2,epil=epil.O2,ko2 = kO2,fatm = Fatm,sat=o2sat,epil.con=epil.O2/vol[1:365,1]))
  return(list(k_sed=Ksed,hypo.con=hypo_o2/vol[1:ndays,2],
              Flux_hypo = Fhypo/vol[1:ndays,2], Flux_sed=Fsed/vol[1:ndays,2],
              z = z_epil,nep = NEP/vol[1:ndays,1], hypo = hypo_o2,epil=epil.O2,
              ko2 = kO2,fatm = Fatm/vol[1:365,1],
              sat=o2sat,epil.con=epil.O2/vol[1:365,1])
         )
}

###########################################




############result &graph ########################
result<-DO_simulate(nep_input,fked)
write.csv(result,"test.csv")

x <- seq(1,ndays)
plot(x,y = result$epil[x]/vol[x,1],pch =1,xlim=c(0,365),ylim = c(0,50),
     col="blue",main="DO simulation",xlab = "day",ylab = "DO")

x <- seq(mix_day,stratify_day)
points(x,y=result$hypo[x]/vol[x,2],col="black")

x <- seq(mix_day,365)
points(x,y=result$hypo[x]*2/(vol[x,2]+vol[x,1]),col="black")

legend("topright",
       c("epilimnion(model)","hypolimnion(model)","epilimnion(1980)","hypolimnion(1980)"),
       fill=c("blue","black","yellow","purple"),cex = 0.75)

#true sampling data
points(sample_v$avgHypoDO*0.032,col="purple")
points(o2.at.sat.base(temp = sample_v$avgEpiT,altitude = 300),col="yellow")
