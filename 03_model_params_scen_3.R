#############################################################################
# THIS SCRIPT SETS UP THE MODEL PARAMETERS
# CHANGE PARAMETERS HERE FOR SCENARIO MODELLING
###--------------------------------------------------------------------------
# THIS IS SCENARIO 3: all household interventions, radical cure treatment for everyone
####-------------------------------------------------------------------------
#############################################################################

# # sets the working directory to the current location of this file
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## this script is dependent on set_up_sim_code.R
# source('~/set_up_sim_code.R')
## if only this script is run (instead of all via 00_main): uncomment line above

pacman::p_load(mvtnorm, plotly)

current_scen<-scen_3

#now generate malaria mosquito population
# keep track of hh_id, long, lat, activity (in house, travelling), next location

#general parameters
sigmaH<-1/9.9 #expected exposed period in humans is 9.9 days
gammaH<-1/20 #recovery rate (without hypnozoites)
treat_rate<- 1/14 #estimate: people are treated after 14 infectious days
omegaH<- 0.1 # waning immunity #0.00038
sigmaM<-1/14  #expected exposed period in mozzys is 14 days
Mfull<-1/2    #expected time the mosquito is full is 2 days
flight_rate<-5 #10km/day is roughly the maximum distance travelled by a blood fed mosquito in a study over 22h 

# multi-species params
Lv_prob<-0.75 #probability of developing hypnozoites without treatment
LvT_prob<- 0.68 # probability to develop hypnozoites after standard cure treatment
Ldev_rate<-1/30 # rate of developing hypnozoites when not receiving any treatment
# cross_immu<-0.5 #1 is immunity in one species is always in the other too # chance of cross-species immunity
treat_prob <- 0.1 # probability of getting treated
cure_prob <- 0.99 # probability of receiving radical cure (receiving standard cure is 1-cure_prob)
dur_SC <- 3 # duration of standard cure treatment
dur_RC <- 14 # duration of radical cure treatment
no_hypno <- 1/400 #400d average time until decay of hypnozoites
treat_speed <- 1/2 # time until treatment
relapse_prob <- 0.2 #0.2
relapse_rate <- 1/45 # time until relapse

#bite parameters
transMH_per_bite<-0.05 #0.05
transHM_per_bite<-0.47 #0.47
multiple_meal<-0.18
bite_rate<-1 #1  ################# a guess for now the per mosquito per day bite rate within a household
bite_reduction_IRS<-0.3 #0.3
bite_reduction_LLIN<-0.56 #0.56
bite_from_treatment<-bite_rate*c(1,1-bite_reduction_IRS,(1-bite_reduction_IRS)*(1-bite_reduction_LLIN),1-bite_reduction_LLIN,1)

#death parameters
death_rate<-1/(5*7) #assumes one month baseline lifespan
death_prob_in_hh<-0.5 ############ a guess for now
death_IRS<-0.56 #0.56
death_LLIN<-0.19 #0.19
death_from_treatment<-c(death_prob_in_hh,1-(1-death_IRS)*(1-death_prob_in_hh),1-(1-death_LLIN)*(1-death_IRS)*(1-death_prob_in_hh),1-(1-death_LLIN)*(1-death_prob_in_hh),death_prob_in_hh)

#mozzy influx rate map
Xminmax=c(0,max(ind_level_data$long))
Yminmax=c(0,max(ind_level_data$lat))
#mozzie_gen<-mvrnorm(100,c(mean(latlonggrid_daco$x),mean(latlonggrid_daco$y)),var(long_lat_mat))
nesting_centres=matrix(NA,nrow=3,ncol=2)
nesting_centres[1,]<-c(25,15)
nesting_centres[2,]<-c(15,20)
nesting_centres[3,]<-c(30,30)
centre_magnitude<-c(5,5,5)

#later include wind, humidy and nesting prevention methods
mozzy_map<- function(centres,centre_magnitudes,Xlims,Ylims,gridpoints,var_num){
  Xvec<-seq(from=Xlims[1],to=Xlims[2],by=((Xlims[2]-Xlims[1])/(gridpoints-1)))
  Yvec<-seq(from=Ylims[1],to=Ylims[2],by=((Ylims[2]-Ylims[1])/(gridpoints-1)))
  
  Xmesh<-matrix(Xvec,nrow=gridpoints,ncol=gridpoints)
  Ymesh<-t(matrix(Yvec,nrow=gridpoints,ncol=gridpoints))
  
  Xlong<-as.vector(Xmesh)
  Ylong<-as.vector(Ymesh)
  
  m_map<-matrix(0,nrow=gridpoints^2,ncol=1)
  for(ii in 1:length(centre_magnitudes)){
    ###The mean and variances will change later based on environment
    m_map<-m_map+(centre_magnitudes[ii]/dmvnorm(centres[ii,],centres[ii,],var_num*diag(2)))*dmvnorm(cbind(Xlong,Ylong),centres[ii,],var_num*diag(2))
  }
  #XYZ<-cbind(Xlong,Ylong,m_map)
  XYZ<-cbind(Xvec,Yvec,t(matrix(m_map,nrow=gridpoints,ncol=gridpoints)))
  return(XYZ)
}

temp<-mozzy_map(nesting_centres,centre_magnitude,Xminmax,Yminmax,100,25)
plot_ly(z= ~temp[,3:102],x= ~temp[,1],y= ~temp[,2],type="surface")

mozzy_hh_imports<- function(long_lat,centres,centre_magnitudes,var_num){
  m_map<-matrix(0,nrow=length(long_lat[,1]),ncol=1)
  for(ii in 1:length(centre_magnitudes)){
    #The mean and variances will change later based on environment
    m_map<-m_map+(centre_magnitudes[ii]/dmvnorm(centres[ii,],centres[ii,],var_num*diag(2)))*dmvnorm(long_lat,centres[ii,],var_num*diag(2))
  }
  #XYZ<-cbind(Xlong,Ylong,m_map)
  XYZ<-cbind(long_lat,m_map)
  return(XYZ)
}

# these parameters are specified by some rough inference in relative_suitability_analysed.R

seas_param_a<-bite_rate*transMH_per_bite

#seas_param_a <- 0.01
seas_param_alph <- 0.065

# this is the suitability map scaled based on inteference
mean_seasonality<-rowMeans(seasonal_MAP_hh)
scaled_seasonal_MAP <-seas_param_a*(seasonal_MAP_hh^(seas_param_alph))


hh_import_rate<-mozzy_hh_imports(people_locations,nesting_centres,centre_magnitude,25)
ggplot(people_locations,aes(x,y))+geom_point(aes(fill=hh_import_rate[,3]),colour="black",pch=21,show.legend=FALSE)+ scale_fill_continuous(type = "gradient")

#does bite rate actually work like this? or should  
hh_bite_prop<-((1-bite_reduction_IRS)^(1*(hh_treatment==2 | hh_treatment==3)))*((1-bite_reduction_LLIN)^(1*(hh_treatment==4 | hh_treatment==3)))
hh_deathonbite<-1-((1-death_IRS)^(1*(hh_treatment==2 | hh_treatment==3)))*((1-death_LLIN)^(1*(hh_treatment==4 | hh_treatment==3)))
mozzy_creation<-hh_import_rate[,3]*transHM_per_bite*bite_rate*hh_bite_prop*(1-hh_deathonbite)/N

propto_mozzycreation <-(transHM_per_bite/transMH_per_bite)*hh_bite_prop*(1-hh_deathonbite)/N
seasonal_creation<-matrix(0,nrow=sum(N),ncol=12)

for(ii in 1:length(hh_labels)){
  for(jj in 1:N[ii]){
    #specifies the row of the data set
    if(ii!=1){
      row_num<-sum(N[1:(ii-1)])+jj
    }else{
      row_num<-jj
    }
    hh_row_indicies<-which(hh_labels==ind_level_data$hh_id[row_num])
    
    
    ind_level_data$m_mosquito_rate[row_num]<-mozzy_creation[hh_row_indicies]
    
    seasonal_creation[row_num,] <- scaled_seasonal_MAP[hh_row_indicies,]*propto_mozzycreation[hh_row_indicies]
  }
}