##############################################################################
# this script was written by Camelia Walker, University of Melbourne
# and only minimally adapted by Leandra Bräuninger
# it sets up the mosquito suitability map and individual level malaria dataframe 
##############################################################################

pacman::p_load(raster)


#find all hh id's
hh_labels<-unique(rbind(matrix(hh_full_case_data$id_household,nrow=length(hh_full_case_data$id_household)),matrix(hh_data$id_household,nrow=length(hh_data$id_household))))



people_locations<-matrix(NA,nrow=length(hh_labels),ncol=2)
indicie_MAP<-matrix(NA,nrow=length(hh_labels),ncol=2)
hh_treatment<-matrix(1,nrow=length(hh_labels),ncol=1)
got_malaria<-matrix(NA,nrow=length(hh_labels),ncol=1)
N<-matrix(NA,nrow=length(hh_labels),ncol=1)
N_guess<-matrix(FALSE,nrow=length(hh_labels),ncol=1)


#get the raster
suitability_map<-raster("input_data/malaria suitability rasters/Month1.gri")
#find the x and y grid of the raster
Xvector<-xFromCol(suitability_map)
Yvector<-yFromRow(suitability_map)

Xbounds<-c(107,107.5)
Ybounds<-c(11.85,12.35)

therows<-which(Yvector>=Ybounds[1] & Yvector<=Ybounds[2])
thecols<-which(Xvector>=Xbounds[1] & Xvector<=Xbounds[2])

Xvector<-Xvector[thecols]
Yvector<-Yvector[therows]

#find a relevant vector of locations and treatments
for(ii in 1:length(hh_labels)){
  temp<-hh_data$longitude[which(hh_data$id_household==hh_labels[ii])]
  
  if(length(temp)==0){
    indicie<-which(hh_full_case_data$id_household==hh_labels[ii])
    temp_loc<-cbind(hh_full_case_data$longitude[indicie],hh_full_case_data$latitude[indicie])
    N_guess[ii]<-TRUE
    
    #this is the minimum possible number of individuals in the household
    N[ii]<-length(unique(hh_full_case_data$id_individual[indicie]))

    
  }else{
    indicie<-which(hh_data$id_household==hh_labels[ii])
    temp_loc<-cbind(temp,hh_data$latitude[indicie])
    treatments<-Treatment_type[indicie]
    Ns<-hh_data$hh_number[indicie]
    N[ii]<-Ns[1]
    hh_treatment[ii]<-treatments[1]
    }
  people_locations[ii,]<-temp_loc[1,]
  
  indicie_MAP[ii,]<-c(which.min((temp_loc[1,2]-Yvector)^2),which.min((temp_loc[1,1]-Xvector)^2))
}


#a factor to scale the suitability (just so we can differentiate between the bite rate and travel rate, whereas the creation rate is a product)
MAP_rate_scaling<-1

val_mat<-matrix(NA,nrow=length(therows),ncol=length(thecols))
season_val_mat<-data.frame(matrix(0,nrow=length(therows),ncol=(length(thecols)+1)))
for(ii in 1:length(therows)){
  val_mat[ii,]<-getValues(suitability_map,row=therows[ii])[thecols]
}
norm_const<-((length(Xvector)*length(Yvector))/sum(val_mat))**MAP_rate_scaling
val_mat<-val_mat*norm_const

season_val_mat[1:length(therows),1:length(thecols)]<-val_mat
season_val_mat[1:length(therows),length(thecols)+1]<-1

seasonal_MAP_hh<-matrix(NA,nrow=length(indicie_MAP[,1]),ncol=12)

seasonal_MAP_hh[,1]<-val_mat[indicie_MAP]
for(ii in 2:12){
  suitability_map<-raster(paste("input_data/malaria suitability rasters/Month",as.character(ii),".gri",sep=""))
  for(jj in 1:length(therows)){
    val_mat[jj,]<-getValues(suitability_map,row=therows[jj])[thecols]
  }
  val_mat<-val_mat*norm_const
  season_val_mat[(length(therows)*(ii-1)+1):(length(therows)*ii),1:length(thecols)]<-val_mat
  season_val_mat[(length(therows)*(ii-1)+1):(length(therows)*ii),length(thecols)+1]<-ii
  seasonal_MAP_hh[,ii]<-val_mat[indicie_MAP]
}

people_locations_LATLONG<-data.frame(cbind(people_locations,hh_treatment))
names(people_locations_LATLONG)<-c("long","lat","treatments")
  
people_locations<-latlong2grid(people_locations)  
people_locations[,1]<-people_locations[,1]-min(people_locations[,1])
people_locations[,2]<-people_locations[,2]-min(people_locations[,2])

#now to make everything in terms of a vector of people (as opposed to hhs)
# hh_id | person id | treatment | longatude | latitude | N | sick 
ind_level_data<-data.frame(matrix(NA,nrow=sum(N),ncol=8))
names(ind_level_data)<-c("hh_id","ind_id","treatment","long","lat","N","malaria","m_mosquito_rate")
seasonal_MAP_ind<-matrix(NA,nrow=sum(N),ncol=12)


for(ii in 1:length(hh_labels)){
  for(jj in 1:N[ii]){
    #specifies the row of the data set
    if(ii!=1){
      row_num<-sum(N[1:(ii-1)])+jj
    }else{
      row_num<-jj
    }
    ind_level_data$hh_id[row_num]<-hh_labels[ii]
    ind_level_data$ind_id[row_num]<-jj
    ind_level_data$treatment[row_num]<-hh_treatment[ii]
    ind_level_data$long[row_num]<-people_locations[ii,1]
    ind_level_data$lat[row_num]<-people_locations[ii,2]
    ind_level_data$N[row_num]<-N[ii]
    seasonal_MAP_ind[row_num,]<-seasonal_MAP_hh[ii,]
  }
  
  indicie<-which(hh_full_case_data$id_household==hh_labels[ii])
  if(sum(indicie)!=0){
    #unique list of IDs in the household
    inds_sick<-unique(hh_full_case_data$id_individual[indicie])
    for(jj in 1:length(inds_sick)){
      
      if(ii!=1){
        row_num<-sum(N[1:(ii-1)])+jj
      }else{
        row_num<-jj
      }
      #for each ID in the household check list of malaria test results (only consider last result)
      malaria_ind<-hh_full_case_data$malaria_type[which(hh_full_case_data$id_individual==inds_sick[jj])]
      # print(ii)
      # print(jj)
      # print(row_num)
      if(malaria_ind[length(malaria_ind)]=="P.f"){
        ind_level_data$malaria[row_num]<-1
      }else if(malaria_ind[length(malaria_ind)]=="P.m"){
        ind_level_data$malaria[row_num]<-2
      }else if(malaria_ind[length(malaria_ind)]=="P.v"){
        ind_level_data$malaria[row_num]<-3
      }
    }

  }

  
}
  
temp<-ind_level_data$malaria
temp<-temp[!is.na(temp)]
vill_colours<-map2color(temp,rainbow(200),limits=c(1,10))

