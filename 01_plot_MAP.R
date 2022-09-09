##############################################################################
# this script was written by Camelia Walker, University of Melbourne
# and only minimally adapted by Leandra Bräuninger
# it visualises the input data in space
##############################################################################

pacman::p_load(tidyverse, RColorBrewer, ggpubr, latex2exp, SpatialEpi)


hh_data<-read.csv("input_data/Malaria_cols.csv")

#colour mapping function
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#create colours for accuracy data
vill_colours<-map2color(hh_data$id_household,rainbow(200),limits=c(1,10))


hh_full_case_data<-read.csv("input_data/M_cases.csv")


#create colours for accuracy data

hh_case_data<-hh_full_case_data[hh_full_case_data$malaria_type!="",]
Mtype_to_num<-matrix(data=NA,nrow=length(hh_case_data$malaria_type),ncol=1)
Mtype_to_num[hh_case_data$malaria_type=="P.f",]=1
Mtype_to_num[hh_case_data$malaria_type=="P.m",]=2
Mtype_to_num[hh_case_data$malaria_type=="P.v",]=3
vill_colours<-map2color(Mtype_to_num,rainbow(200),limits=c(1,10))

latlonggrid<-latlong2grid(cbind(hh_data$longitude,hh_data$latitude))

latlonggrid_case<-latlong2grid(cbind(hh_case_data$longitude,hh_case_data$latitude))
latlonggrid_case[,1]<-latlonggrid_case[,1]-9300
latlonggrid[,1]<-latlonggrid[,1]-9300
latlonggrid_case[,2]<-latlonggrid_case[,2]-1320
latlonggrid[,2]<-latlonggrid[,2]-1320

hh_map <- ggplot(latlonggrid_case, aes(x,y))+geom_point(aes(x,y),latlonggrid,pch=24,colour="grey")+geom_point(aes(fill=vill_colours),colour="black",pch=21)+ scale_fill_discrete(name = "Malaria Species", labels = c("P.v.", "P.f.", "P.m."))+labs(x = "km (East from arbitrary point)",y= "km (North from arbitrary point)")
ggsave("output/hh_map.png", hh_map, width = 5, height = 5)

#ggplot(hh_data, aes(longitude,latitude))+geom_point(aes(fill=accuracy),colour="black",pch=21,size=(hh_data$accuracy/4))

Treatment_type<-matrix(data=NA,nrow=length(hh_data$malaria_prevention),ncol=1)
Treatment_type[hh_data$malaria_prevention=="Don't know",]=1
Treatment_type[hh_data$malaria_prevention=="IRS",]=2
Treatment_type[hh_data$malaria_prevention=="IRS + ITN",]=3
Treatment_type[hh_data$malaria_prevention=="ITN",]=4
Treatment_type[hh_data$malaria_prevention=="None",]=5
type_colours<-map2color(Treatment_type,rainbow(200),limits=c(1,5))

treat_map <- ggplot(latlonggrid, aes(x,y))+geom_point(aes(fill=type_colours),colour="black",pch=21)+ scale_fill_discrete(name = "Treatment", labels = c("IRS + ITN", "ITN", "IRS","Unknown","None"))+labs(x = "km (East from arbitrary point)",y= "km (North from arbitrary point)")
ggsave("output/treat_map.png", treat_map, width = 5, height = 5)
