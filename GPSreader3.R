#### a code for the Ibex GPS data ###
#reads a pre-processed xls file, now the one Jonathan sent me 
#filters a few points
#plots the tracks with ggplot
#analyse the data some
#compares male-female
#Orr Spiegel, 2018,2019, 2020

### load packages  ####
require(readxl)
library(ggplot2)
library(ggmap)
require(adehabitatHR);#loads also: require(CircStats); #require(boot) ; require(MASS) ;  
require(geosphere)# for distance to lines (Roads)
require(rgdal)
require(lme4)#for the stat models at the end
require(AICcmodavg)#just for the model comparison table

#require(writexl)#for exporting

### set relevant parameters ####
HomePc=01;#1 for home 0 for lab, for path
ToPlot=0# for the road distance figures
MinDistToRoadThreshhold=50#this is the threshold distance to consider as near by the road
PossibeElev=c(-420,300); #possible elevation range, values more extreme will be filtered
MaxDop=5#anityhinh with DOP above this is deleted
DataSetToUse=3#choosing which dataset to use for the DailyData later 1- just 1h, 3- just h 3\4 long samples,5 -all days
MinFixPerDay=c(15,5,5)# days with less gps fixes will not be considered in the DailyData, for the whole, 1h and 4hour interval datasets
IbexToExcludeAlltogether=c('Menachem_1','NoName','Ezer')#these two have very short duration, so remove from start

Proj="+proj=utm +zone=36 +ellps=WGS84"##NEED TO check cooridnate system
#require(plyr);#require(lme4);require(AICcmodavg); require(fields); 
#require(rgeos);require(sp);#


### set path and read the NEW xls file after some manual processing ####
if (HomePc==1) {
  #home pc-to update
  setwd('D:\\OrrS2\\Box Sync\\R codes\\IbexCodes')#home PC
  Path1="D:\\OrrS2\\Box Sync\\IBEX\\GPSdata\\All_GPS_Combined_ed2020.xlsx"#the file to work on 
    }else{
  #lab pc
  setwd('D:\\OrrS\\Documents\\Box Sync\\R codes\\IbexCodes')#Lab PC
  Path1="D:\\OrrS\\Documents\\Box Sync\\IBEX\\GPSdata\\All_GPS_Combined_ed2020.xlsx"#the file to work on 
}


GPSdatRaw <- read_excel(Path1,col_names=T)#, col_types = c("text","text","text", "text", "text","text","numeric","numeric","numeric","numeric"))

#setwd('C:\\Users\\Orr\\Box Sync\\IBEX\\R_codes')#home
#GPSdatRaw <- read_excel("C:\\Users\\Orr\\Box Sync\\IBEX\\GPSdata\\CombinedGPSdataProccedUpdtApr2018.xlsx",col_names=T)#, col_types = c("text","text","text", "text", "text","text","numeric","numeric","numeric","numeric"))

#View(GPSdatRaw)                      

names(GPSdatRaw)=make.names(names(GPSdatRaw), unique=TRUE);#making sure all names are without spaces
#names(GPSdatRaw) 


### subsetting the columns I need, and correcting formats ####
#GPSdatFltr=GPSdatRaw[,c('AnimalID','CollarID',"UTC_Date","UTC_Time",'FixType',''
#                     "Latitude","Longitude","Height","DOP","Temp","Easting","Northing")]
GPSdatFltr=GPSdatRaw

## getting rid of a few specific ibex with shrot tracking
GPSdatFltr=GPSdatFltr[!GPSdatFltr$AnimalID %in% IbexToExcludeAlltogether,]

##setting colums to be the right datatype, or time format:
GPSdatFltr$AnimalID=as.factor(GPSdatFltr$AnimalID)
GPSdatFltr$Sex=as.factor(GPSdatFltr$Sex)
GPSdatFltr$Season=as.factor(GPSdatFltr$Season)
GPSdatFltr$SmplingIntrvlToday=NA#adding column for storage 

GPSdatFltr$UTC_Time=as.POSIXct(as.numeric(as.POSIXct(GPSdatFltr$UTC_Time)) %% 86400, origin = "2000-01-01",tz ='UTC')
GPSdatFltr$UTC_TimeOnly=strftime(GPSdatFltr$UTC_Time, format="%H:%M:%S",tz ='UTC')
GPSdatFltr$UTC_Date.Time=as.POSIXct(paste(GPSdatFltr$UTC_Date, GPSdatFltr$UTC_TimeOnly), format="%Y-%m-%d %H:%M:%S",tz="UTC")

## getting rid of dupicate lines of the same one
GPSdatFltr=GPSdatFltr %>% dplyr::distinct(CollarID,AnimalID,DateTime,Latitude,Longitude, .keep_all = TRUE)

#GPSdatFltr$TimeDiff=difftime(t1, t2, units = "minutes")
str(GPSdatFltr)


### printing how many animals and datapoints #####
print('using the file All_GPS_Combined_ed2020.xlsx (already partially filtered) ')
print(paste('we have',length(levels(GPSdatFltr$AnimalID)),'individuals in the data', sep=' '))
print(paste('we have',dim(GPSdatFltr)[1],'fixes in the data', sep=' '))

### fitering by elevation and DOP #####
hist(GPSdatFltr$DOP[GPSdatFltr$DOP<MaxDop],0:10,main='DOP values after dropping DOP>=10')
print(paste('dropping out ',sum(GPSdatFltr$DOP>=MaxDop),'observations with DOP>=MaxDop and elevation not between -420 to +300', sep=' '))
GPSdatFltr=GPSdatFltr[which(GPSdatFltr$DOP<MaxDop & GPSdatFltr$FixType=='GPS-3D' & #only low DOP and fixtype 3D!
                              GPSdatFltr$Height>=PossibeElev[1] & GPSdatFltr$Height<=PossibeElev[2]),]
hist(GPSdatFltr$Height, breaks=seq(from=-500, by=25, to=(50+max(GPSdatFltr$Height))),main='Elevation values after dropping extremes -420 +300')
hist(GPSdatFltr$DOP,0:MaxDop,main='DOP values after dropping DOP>=MaxDop')
Indx=seq(from=2,to=dim(GPSdatFltr)[1])

#dropping filtered levels
GPSdatFltr=droplevels(GPSdatFltr)
GPSdatFltr$FixType=NULL#not needed anymore

### sorting by collars, animals, dates& time ####
GPSdatFltr=GPSdatFltr[   order(GPSdatFltr$CollarID,GPSdatFltr$AnimalID,GPSdatFltr$DateTime,decreasing=F),]
#GPSdatFltr[with(GPSdatFltr, order(CollarID, b)), ]
#distinct()


### Estimating Time Difference between fixes (effective sampling intervals) ####
#updating new animal column
GPSdatFltr$IsNewAnimal=c(1,diff(as.numeric(GPSdatFltr$AnimalID)))
#GPSdatFltr$tt=c(1,pracma::strcmpi(as.character(GPSdatFltr$AnimalID)))
GPSdatFltr$IsNewAnimal[GPSdatFltr$IsNewAnimal>0 | GPSdatFltr$IsNewAnimal<0]=1

#GPSdatFltr$IsNewAnimal=c(TRUE, #the first line is TRUE, start comparison from the second line
#                         sapply(Indx,#loop on all values
#                                function(i){if(as.numeric(GPSdatFltr$AnimalID)[i]==as.numeric(GPSdatFltr$AnimalID)[i-1])  #is the ID number the same?
#                                {FALSE}else{TRUE}}))#store FALSE for same animal TRUE if a new one

GPSdatFltr$SameAnimaDay=c(1, #the first line is one, start comparison from the second line
                 sapply(Indx,#loop on all values
                        function(i){if((as.numeric(GPSdatFltr$AnimalID)[i]==as.numeric(GPSdatFltr$AnimalID)[i-1]) & #is the ID number the same?
                                       (as.numeric(GPSdatFltr$UTC_Date)[i]==as.numeric(GPSdatFltr$UTC_Date)[i-1])) #is date the same
                                    {1}else{0}}))#store 1 for same 0 for not

GPSdatFltr$TimeDiff=c(NA,#the first line is missing, start comparison from the second line
                       sapply(Indx,#loop on all values 
                              function(i){if(GPSdatFltr$IsNewAnimal[i]!=TRUE)#calculating only for the same animal  but also across days. change to --> for only within days if(GPSdatFltr$SameAnimaDay[i]==1)#is it the same day and tag?
                                {difftime(GPSdatFltr$UTC_Date.Time[i],GPSdatFltr$UTC_Date.Time[i-1], units = 'secs')}#compare times
                                else{NA}})  ) #dont compare times

if (length(which(GPSdatFltr$TimeDiff<0))>0)print('there is a problem with the data, check negative time differences!!!!!!')
#making sure all complete cases for the more relevant columns
GPSdatFltr=GPSdatFltr[complete.cases(GPSdatFltr[,c("AnimalID","CollarID","UTC_Date.Time","Easting","Northing","Height","Longitude","Latitude")]),]; 
GPSdatFltr$CollarID=as.factor(GPSdatFltr$CollarID)

GPSdatFltr$DayFromBegining=1+ difftime(GPSdatFltr$UTC_Date,min(GPSdatFltr$UTC_Date),units='days')
#
#saving days as burst with a different name for each indiv and days since the dataset starts
GPSdatFltr$DaysAsBurst=as.factor(paste(GPSdatFltr$AnimalID,(1+difftime(GPSdatFltr$UTC_Date,min(GPSdatFltr$UTC_Date),units='days')),sep="_"));
#tt2=tt[GPSdatFltr$SameAnimaDay==0 & GPSdatFltr$IsNewAnimal==0]
#unique(diff(tt2))



### Minimal distance to Paved roads  ####
# list the kmz files in a given folder path
KMZs <- list.files( pattern="*.kmz", full.names=FALSE)#path="Data/kmz-files",

#read the kmz file into a list of different lines
#library(maptools)LonLat=getKMLcoordinates(kmlfile = unzip(KMZs),ignoreAltitude = TRUE)
KMLfile = unzip(KMZs)
LonLat=maptools::getKMLcoordinates(kmlfile = KMLfile,ignoreAltitude = TRUE)

#a small helper function to extract names
lyr <- ogrListLayers(KMLfile)
mykml <- list()
for (i in 1:length(lyr)) {mykml[i] <- readOGR(KMLfile,lyr[i]) }
names(mykml) <- lyr
#RoadSectionNames=mykml$PavedRoadsIbex$Name; rm(KMLfile,mykml,lyr)
names(LonLat)=mykml$PavedRoadsIbex$Name; rm(KMLfile,mykml,lyr)

#to work only on the main road sections
#LonLat[which(names(LonLat) %in% c("Ovnat","KaliaQumeram","Dragot","MitzpeShalem","EinGedy2","EinGedy1" ,"Massada" ))] <- NULL


#Ibex Locations GPS
pnts=as.data.frame(GPSdatFltr[c('Longitude','Latitude')]);print('how many points we have now in filtered dataset?');dim(pnts)[1]
DistanceMatrix=matrix(NA, ncol  = length(LonLat), nrow = dim(pnts)[1])#here i will store the data on min ditance to each polygon

for (lineCnt in 1:length(LonLat) ){#loop on lines in the kml;
  #extacting the currentl line
  line=data.frame(LonLat[[lineCnt]]);colnames(line)=c('Lon','Lat')
  #distances to current line
  distances=dist2Line(pnts, line, distfun=distHaversine)
  DistanceMatrix[,lineCnt]=distances[,1]  #@distance in meteres
  
  
  #plotting
  if (ToPlot==1){#plot the distances for 1000 points for each road section
    #par(pty="s")#square aspect ratio
    plot( makeLine(line), type='l', asp=1,main=names(LonLat)[lineCnt])
    points(line)
    SamplePoints=sample(x=(1:dim(pnts)[1]),size=100)
    points(pnts[SamplePoints,], col='blue', pch=20)
    points(distances[SamplePoints,2], distances[SamplePoints,3], col='red', pch='x')
    #for all points: for (i in 1:nrow(distances)) lines(gcIntermediate(pnts[i,], distances[i,2:3], 10), lwd=2)
    for (i in 1:length(SamplePoints)) lines(gcIntermediate(pnts[SamplePoints[i],], distances[SamplePoints[i],2:3], 10), lwd=2, col='darkgreen')
    }
  
}#loop on lines in the kml
rm(line,lineCnt,ToPlot)

#for each GPS fix find the minimal distance to a paved road and check if below thresh hold
#View(DistanceMatrix)
#Nearerst distance to roads
GPSdatFltr$MinDistToRoad=apply(DistanceMatrix, 1, FUN=min, na.rm=T)
GPSdatFltr$NearbyPavedRoad=GPSdatFltr$MinDistToRoad<=MinDistToRoadThreshhold

#only the main road (90, Arad, Dimona) not including side branches
#DistanceMatrix2=DistanceMatrix[,c(1,6,10,11)]
#GPSdatFltr$MinDistToRoad=apply(DistanceMatrix2, 1, FUN=min, na.rm=T)


rm(DistanceMatrix,SamplePoints,pnts,LonLat,distances)  


### setting ByIbex A dataframe By ibex for sampling interval as a function of day and detect the effective sampling interval####
ByIbex=as.data.frame(levels(GPSdatFltr$AnimalID));names(ByIbex)='AnimalID'
ByIbex$CollarID= NA; ByIbex$Sex=NA;ByIbex$DaysOfTracking= NA;  
ByIbex$FirstDay=NA;  ByIbex$lastDay=NA;
ByIbex$SwitchDayFrom1hInterval= NA; ByIbex$DaysIn1hInterval=NA;  
ByIbex$longerIntervalinHours=NA;ByIbex$DaysInlongerInterval=NA;ByIbex$SamplingIntervals=NA; ByIbex$DatesThisIbex=NA;
ByIbex$NofFixes=NA;ByIbex$FixesNearRoad=NA; ByIbex$MeanDistToRoad=NA;

# manually correcting the smapling protocols for a few ibex
MenachemIndx=which(ByIbex$AnimalID=='Menachem_2')#s8 Mar 2019 to 17 Apr 2019 in 1 h interval

#loop on indivuduals +plot of sampling frequency
for (IbxCntr in 1:length(levels(GPSdatFltr$AnimalID))){
   
  indxthisIbx=which(GPSdatFltr$AnimalID==levels(GPSdatFltr$AnimalID)[IbxCntr])#these are indices of this ibex in the dataset
  DatesThisIbex=unique(GPSdatFltr$UTC_Date[indxthisIbx])
  print(paste(levels(GPSdatFltr$AnimalID)[IbxCntr], 'has data from ',length(DatesThisIbex) ,'days'))
  
  ##loop on days to find sampling interval for each day
  SamplingIntervals=sapply(DatesThisIbex,function(i){mean(na.rm=T,GPSdatFltr$TimeDiff[GPSdatFltr$UTC_Date==i &  GPSdatFltr$AnimalID==levels(GPSdatFltr$AnimalID)[IbxCntr]])})#find the places it the right ibex and day, mean the time differences
  
  #find switch point from 1 hour to longer intervals (3 or 4 hours)
  SwitchIndx=which(SamplingIntervals/3600>=3)[1];if(is.na(SwitchIndx)){SwitchIndx=length(SamplingIntervals)}
  
  #add manual correction for other days here???
  if (IbxCntr==which(ByIbex$AnimalID=='Eshkol')){SwitchIndx=54}
  if (IbxCntr==which(ByIbex$AnimalID=='Eva')){SwitchIndx=40}  
  if (IbxCntr==which(ByIbex$AnimalID=='Yizhar')){SwitchIndx=41}  
  if (IbxCntr==which(ByIbex$AnimalID=='Menachem_2')){SwitchIndx=144}  
  
  
  #setting longer interval
  LongerInterval=round(median(SamplingIntervals[(SwitchIndx+1):length(SamplingIntervals)])/3600,digits=2)
  
  #plotting changes in GPS sampling intervals:
  plot(x=DatesThisIbex,y=SamplingIntervals/3600, type="l",lwd=2, 
       xlab='tracking date',ylab=('sampling interval (h)'),ylim=c(0,8),
       main=paste('mean daily sampling interval for',levels(GPSdatFltr$AnimalID)[IbxCntr])) 
  abline(v=DatesThisIbex[SwitchIndx], col='red',lty=2)
  abline(h=LongerInterval, col='blue',lty=2)
  abline(h=3, col='green',lty=2)
  
  #storing in the dataframe:
  ByIbex$CollarID[IbxCntr]=as.character(GPSdatFltr$CollarID[indxthisIbx[1]])
  ByIbex$DaysOfTracking[IbxCntr]=length(DatesThisIbex)
  ByIbex$FirstDay[IbxCntr]=as.Date(DatesThisIbex[1],origin = "2000-01-01")
  ByIbex$lastDay[IbxCntr]=as.Date(tail(DatesThisIbex,1),origin = "2000-01-01")
  ByIbex$SwitchDayFrom1hInterval[IbxCntr]=as.Date(DatesThisIbex[SwitchIndx],origin = "2000-01-01")
  ByIbex$DaysIn1hInterval[IbxCntr]=DatesThisIbex[SwitchIndx]-DatesThisIbex[1]
  ByIbex$longerIntervalinHours[IbxCntr]=LongerInterval;
  ByIbex$DaysInlongerInterval[IbxCntr]=tail(DatesThisIbex,1)-DatesThisIbex[SwitchIndx]
  ByIbex$DatesThisIbex[IbxCntr]=list(DatesThisIbex)
  ByIbex$SamplingIntervals[IbxCntr]=list(SamplingIntervals)
  ByIbex$NofFixes[IbxCntr]=length(GPSdatFltr$NearbyPavedRoad[indxthisIbx])
  ByIbex$FixesNearRoad[IbxCntr]=sum(GPSdatFltr$NearbyPavedRoad[indxthisIbx])
  ByIbex$MeanDistToRoad[IbxCntr]=mean(GPSdatFltr$MinDistToRoad[indxthisIbx],na.rm = T)
  ByIbex$MeanDistToRoad[IbxCntr]=mean(GPSdatFltr$MinDistToRoad[indxthisIbx],na.rm = T)
  ByIbex$Sex[IbxCntr]=as.character(GPSdatFltr$Sex[indxthisIbx[1]])
  
  #manual correction for Menachem2
  if (IbxCntr==which(ByIbex$AnimalID=='Menachem_2')){ 
    ByIbex$DaysIn1hInterval[IbxCntr]=42#since this is the duration of short interval
    ByIbex$DaysInlongerInterval[IbxCntr]=ByIbex$DaysInlongerInterval[IbxCntr]+101#since he has also 100 days before with long interval
    }  
  
  #adding the sampling interval for each day in the main data.frame
  GPSdatFltr$SmplingIntrvlToday[indxthisIbx]=1;#set all dates of this individual to 1 h
  if (!is.na(LongerInterval)){#check if there are days with longer time interval
    GPSdatFltr$SmplingIntrvlToday[intersect(which(GPSdatFltr$UTC_Date>=DatesThisIbex[SwitchIndx]),indxthisIbx)]=LongerInterval}#update those dates
  
  #manual correction for Menachem2- setting the first long interval period as long
  if (IbxCntr==which(ByIbex$AnimalID=='Menachem_2')){
    #beggining--> the last day of the first long interval period: 
    GPSdatFltr$SmplingIntrvlToday[indxthisIbx[1: tail(which(as.character(GPSdatFltr$UTC_Date[indxthisIbx])=='2019-03-07'),n=1)]] =LongerInterval
  }#manual correction
    #sort(unique(GPSdatFltr$UTC_Date[which(GPSdatFltr$UTC_Date>=DatesThisIbex[SwitchIndx])]))
  #unique(GPSdatFltr$AnimalID[GPSdatFltr$SmplingIntrvlToday==1])
  #unique(GPSdatFltr$UTC_Date[intersect(which(GPSdatFltr$UTC_Date>=DatesThisIbex[SwitchIndx]),indxthisIbx)])
  rm(DatesThisIbex,indxthisIbx,LongerInterval,SwitchIndx)
}

#correcting data.frame back to date format not sure why the 30 years difference needed
ByIbex$FirstDay=as.Date(ByIbex$FirstDay,origin = "1970-01-01",format="%Y-%m-%d")
ByIbex$lastDay=as.Date(ByIbex$lastDay,origin = "1970-01-01",format="%Y-%m-%d")
ByIbex$SwitchDayFrom1hInterval=as.Date(ByIbex$SwitchDayFrom1hInterval,origin = "1970-01-01",format="%Y-%m-%d")
#ByIbex$Sex1m2f=c(1,1,2,2,1,2,1,1,1,2,1,1,1,2,2,1,1)#here we need to update the sex for each individual in the data
GPSdatFltr$SmplingIntrvlToday=floor(GPSdatFltr$SmplingIntrvlToday)#going to the nearest lower full hour sampling
ByIbex$CollarID=as.factor(ByIbex$CollarID)



print('matadata on ibext tracking data')
print(ByIbex[,1:8])
ByIbex$Sex=as.factor(ByIbex$Sex)

## creating a second data frame without all the details 
ByIbex2=ByIbex[,c("AnimalID",'CollarID',"Sex","DaysOfTracking","DaysIn1hInterval",'DaysInlongerInterval',"NofFixes","FixesNearRoad","MeanDistToRoad")]
#ByIbex2=ByIbex[,-c("SwitchDayFrom1hInterval",)]
ByIbex2$PropNearRoad=ByIbex2$FixesNearRoad/ByIbex2$NofFixes
#print(ByIbex2)

### Subsetting only the 1hour and 3/4h intervals dates for each ibex ####
GPSdat1h=subset(GPSdatFltr,SmplingIntrvlToday==1)
GPSdat3or4h=subset(GPSdatFltr,SmplingIntrvlToday>=3)

#dropping filtered levels
GPSdat1h=droplevels(GPSdat1h)
GPSdat3or4h=droplevels(GPSdat3or4h)

### Elevation annotation - no working yet####
#demo in : https://terpconnect.umd.edu/~egurarie/teaching/MovementAtICCB2017/AnnotatingData.html
class(GPSdat1h)

tt=GPSdat1h[1]
#require(elevatr)
#m2 <- dplyr::mutate(GPSdat1h, elevation = get_elev_point(data.frame(x = Longitude, y = Latitude), 
#            prj = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", src = "aws", api_key = "mapzen-91tYY4s"))
#m2$elevation <- as.numeric(m2$elevation@data[,1])

### plotting tracks with ggplot+ggmap all together #####

#ggplot(GPSdatFltr, aes(x,y,color = AnimalID)) + geom_path() + coord_fixed()

#plotting all ibex together
Mapbox1 <- make_bbox(lon=GPSdatFltr$Longitude,lat=GPSdatFltr$Latitude, f=0.1) #defines the borders of the box
SatImagBox1<- get_map(location=Mapbox1, maptype="satellite", source="google",zoom=10);ggmap(SatImagBox1)
#load("./lizmap.Rdata")
#make plot for each variable of interest
mapTry1 <- ggmap(SatImagBox1) +  theme(plot.margin = margin(0,0,0,0, "mm")) +  labs(x="longitude", y="latitude") 
mapTry1 + geom_point(  data=GPSdatFltr,  alpha = 0.5, aes(x=Longitude, y=Latitude, col=AnimalID))
mapTry1 + geom_path (  data=GPSdatFltr,               aes(x=Longitude, y=Latitude, col=AnimalID,group=AnimalID))


### loop to plot with ggplot+ggmap each ibex separatly #####
#levels(GPSdatFltr$AnimalID)[IbxCntr]
Color=c('red','brown','coral4','blue','chocolate','brown2','darkorange1','blueviolet','coral','darkgoldenrod1','cyan2','darksalmon', #the first 12 ibex Fem in Blue, males in red-pink
          'blue','green','rosybrown','royalblue2','darkgreen')#
Color=c(Color,Color,Color)#now becuase there are soo many ibex i double it

for (IbxCntr in 1:length(levels(GPSdatFltr$AnimalID))){
  print(paste(levels(GPSdatFltr$AnimalID)[IbxCntr],'is in color',Color[IbxCntr])) 
  indxthisIbx=which(GPSdatFltr$AnimalID==levels(GPSdatFltr$AnimalID)[IbxCntr])#these are indices of this ibex in the dataset
  ##creating the map for this individal
  Mapbox2 <- make_bbox(lon=GPSdatFltr$Longitude[indxthisIbx],lat=GPSdatFltr$Latitude[indxthisIbx], f=0.4)#defines the borders of the box, now by individual
  SatImagBox2<- get_map(location=Mapbox2, maptype="satellite", source="google",zoom=11)
  mapTry2 <- ggmap(SatImagBox2) +theme(plot.margin = margin(0,0,0,0, "mm")) + labs(x="longitude", y="latitude") +ggtitle(levels(GPSdatFltr$AnimalID)[IbxCntr]) +
    coord_quickmap(xlim=c(Mapbox2[1], Mapbox2[3]),ylim=c(Mapbox2[2], Mapbox2[4]))  
  print(mapTry2 + geom_point(colour = Color[IbxCntr], size = 3, data=GPSdatFltr[indxthisIbx,],  alpha = 0.5, aes(x=Longitude, y=Latitude))) #plotting tracks as points
  print(mapTry2 + geom_path(colour = Color[IbxCntr],            data=GPSdatFltr[indxthisIbx,],               aes(x=Longitude, y=Latitude))) #plotting tracks as a line
}

rm(mapTry2, SatImagBox2,Mapbox2,indxthisIbx)





### converting to Adehabitat object as.ltraj class #####
#Creates class Spatial Points for all locations
#tt=GPSdat1h[which(GPSdat1h$AnimalID=='Alexander') ,];tt=droplevels(tt)

#setting coordinate systems
xysp <- SpatialPoints(GPSdatFltr[,c("Easting","Northing")])
proj4string(xysp) <- CRS(Proj) 
#not sure why this step needed but it worked
GPSdatFltr=data.frame(data.frame(xysp),GPSdatFltr[!names(GPSdatFltr) %in% c("Easting","Northing") ])
coordinates(GPSdatFltr)=data.frame(xysp)
# as.ltraj 
GPSdatFltr_ltrj <- as.ltraj(coordinates(GPSdatFltr),GPSdatFltr$UTC_Date.Time,id=GPSdatFltr$AnimalID,burst = GPSdatFltr$DaysAsBurst,infolocs =GPSdatFltr)
rm(xysp)

xysp <- SpatialPoints(GPSdat3or4h[,c("Easting","Northing")])
proj4string(xysp) <- CRS(Proj) #cooridnate system
#not sure why this step needed but it worked
GPSdat3or4h=data.frame(data.frame(xysp),GPSdat3or4h[!names(GPSdat3or4h) %in% c("Easting","Northing") ])
coordinates(GPSdat3or4h)=data.frame(xysp)
# as.ltraj 
GPSdat3or4h_ltrj <- as.ltraj(coordinates(GPSdat3or4h),GPSdat3or4h$UTC_Date.Time,id=GPSdat3or4h$AnimalID,burst = GPSdat3or4h$DaysAsBurst,infolocs =GPSdat3or4h)
rm(xysp)

xysp <- SpatialPoints(GPSdat1h[,c("Easting","Northing")])
proj4string(xysp) <- CRS(Proj) #cooridnate system
#not sure why this step needed but it worked
GPSdat1h=data.frame(data.frame(xysp),GPSdat1h[!names(GPSdat1h) %in% c("Easting","Northing") ])
coordinates(GPSdat1h)=data.frame(xysp)
# as.ltraj 
GPSdat1h_ltrj <- as.ltraj(coordinates(GPSdat1h),GPSdat1h$UTC_Date.Time,id=GPSdat1h$AnimalID,burst = GPSdat1h$DaysAsBurst,infolocs =GPSdat1h)
rm(xysp)

GPSdatFltrDF=ld(GPSdatFltr_ltrj)
GPSdat3or4hDF=ld(GPSdat3or4h_ltrj)
GPSdat1hDF=ld(GPSdat1h_ltrj)

### DailyData a DataFrame By Day by Ibex #####
#a helper function to put in the right dataset (1n/3h/everything) into DailyData
DailyDataCalculator=function(DataSetToUse){
  #choosing which data to use
  if (DataSetToUse==1){tmpData=GPSdat1hDF;MinFixPerDay=MinFixPerDay[1]} #only th 1 intervals will be used for 1 h data
  if (DataSetToUse==3){tmpData=GPSdat3or4hDF;MinFixPerDay=MinFixPerDay[2]} #only th long smapling intervals will be used for 1 h data
  if (DataSetToUse==5){tmpData=GPSdatFltrDF;MinFixPerDay=MinFixPerDay[3]} #all data
  
  #building it
  DailyData=data.frame(Burst=unique(tmpData$burst))
  DailyData$id=NA;DailyData$date=NA;DailyData$Season=NA; 
  DailyData$Sex=NA;DailyData$DailyN_GPSfix=NA
  DailyData$SumdistM=DailyData$MxDailyDisplcmnt=NA
  DailyData$NetDailyDisplcmnt=NA
  DailyData$TortuosityAllDay=DailyData$TortuosityuntilMax=NA
  DailyData$MeanDayElev=DailyData$DayElevRange=NA
  
  for (burstCnt in 1:length(unique(tmpData$burst))){
    indx=which(as.character(tmpData$burst)==as.character(DailyData$Burst[burstCnt]))#finding the lines of this current burst in the main dataframe
    DailyData$id[burstCnt]=as.character(tmpData$id[indx[1]])
    DailyData$date[burstCnt]=tmpData$UTC_Date[indx[1]]
    DailyData$Season[burstCnt]=as.character(tmpData$Season[indx[1]])
    DailyData$Sex[burstCnt]=as.character(tmpData$Sex[indx[1]])
    DailyData$DailyN_GPSfix[burstCnt]=length(indx)
    DailyData$SumdistM[burstCnt]=sum(tmpData$dist[indx],na.rm=T)
    DailyData$MxDailyDisplcmnt[burstCnt]=sqrt(max(tmpData$R2n[indx],na.rm=T))
    DailyData$NetDailyDisplcmnt[burstCnt]=sqrt(tail(tmpData$R2n[indx],na.rm=T,1))
    if (DailyData$NetDailyDisplcmnt[burstCnt]>0){#making sure that net  displacement is not zero that is undividable
      DailyData$TortuosityAllDay[burstCnt]=DailyData$SumdistM[burstCnt]/DailyData$NetDailyDisplcmnt[burstCnt];
    }else{
      DailyData$TortuosityAllDay[burstCnt]=NA #if net displacement is zero
    }
    
    IndxUntilFurthestPoint=indx[1]:indx[which.max(tmpData$R2n[indx])];#only the locations until the most distant point of the day
    DailyData$TortuosityuntilMax[burstCnt]=sum(tmpData$dist[IndxUntilFurthestPoint],na.rm=T)/#traveled distance until the most distant point
      sqrt(max(tmpData$R2n[IndxUntilFurthestPoint],na.rm=T))#displacment unlit the most furthst point
    
    DailyData$MeanDayElev[burstCnt]=round(mean(tmpData$Height[indx],na.rm=T),digits = 1)#what was the mean elevation
    DailyData$DayElevRange[burstCnt]=range(tmpData$Height[indx],na.rm=T)[2]-range(tmpData$Height[indx],na.rm=T)[1] #what was the elevation range?
    
    rm(IndxUntilFurthestPoint,indx)
    
  }#loop on days by individuals (bursts)
  rm(burstCnt)
  
  ##filtering days with too few fixes
  DailyData=DailyData[DailyData$DailyN_GPSfix>=MinFixPerDay,]
  
  #setting factors again
  DailyData$id=as.factor(DailyData$id)
  DailyData$Season=as.factor(DailyData$Season)
  DailyData$Sex= as.factor(DailyData$Sex)
  DailyData$date=as.POSIXct((DailyData$date), format="%Y-%m-%d ",tz="UTC",origin="1970-01-01")
  DailyData=droplevels(DailyData)
  #DailyData$Sex1m2f=1;#adding sex, all males by defaul
  #DailyData$Sex1m2f[DailyData$id %in% c("Chana","Cheli", "Golda","Nikita", "Pola","Ruth")]=2;#females by names
  #DailyData$Sex1m2f=as.factor(DailyData$Sex1m2f)
  
  
  summary(DailyData)
  return(DailyData)
  rm(tmpData)
}#end of hepler function

DailyData=DailyDataCalculator(DataSetToUse);

#adding a column for a biological guess David wanted.
DailyData$IsMaleRut=0;
DailyData$IsMaleRut[which(DailyData$Sex=='M'&DailyData$Season=='Rut')]=1; #adding a column with 1 or 0 for males during rut
DailyData$IsMaleRut=as.factor(DailyData$IsMaleRut)

## adding the mean daily movement to the byIbex2 dataframe
ByIbex2$NetDailyDisplcmnt=NA;ByIbex2$MxDailyDisplcmnt=NA;ByIbex2$SumdistM=NA;ByIbex2$TortuosityAllDay=NA;ByIbex2$TortuosityuntilMax=NA;

ByIbex2$DailyN_GPSfix     =sapply(as.character(ByIbex2$AnimalID),function(i){mean(na.rm=T,DailyData$DailyN_GPSfix          [as.character(DailyData$id)==i ])})#
ByIbex2$MxDailyDisplcmnt  =sapply(as.character(ByIbex2$AnimalID),function(i){mean(na.rm=T,DailyData$MxDailyDisplcmnt  [as.character(DailyData$id)==i ])})#
ByIbex2$SumdistM          =sapply(as.character(ByIbex2$AnimalID),function(i){mean(na.rm=T,DailyData$SumdistM          [as.character(DailyData$id)==i ])})#
ByIbex2$NetDailyDisplcmnt =sapply(as.character(ByIbex2$AnimalID),function(i){mean(na.rm=T,DailyData$NetDailyDisplcmnt [as.character(DailyData$id)==i ])})#
ByIbex2$TortuosityAllDay  =sapply(as.character(ByIbex2$AnimalID),function(i){mean(na.rm=T,DailyData$TortuosityAllDay  [as.character(DailyData$id)==i ])})#
ByIbex2$TortuosityuntilMax=sapply(as.character(ByIbex2$AnimalID),function(i){mean(na.rm=T,DailyData$TortuosityuntilMax[as.character(DailyData$id)==i ])})#
ByIbex2$MeanDayElev       =sapply(as.character(ByIbex2$AnimalID),function(i){mean(na.rm=T,DailyData$MeanDayElev       [as.character(DailyData$id)==i ])})#
ByIbex2$DayElevRange      =sapply(as.character(ByIbex2$AnimalID),function(i){mean(na.rm=T,DailyData$DayElevRange      [as.character(DailyData$id)==i ])})#
#rm(i)
#tt=ByIbex2[c('MeanDistToRoad','NetDailyDisplcmnt','MxDailyDisplcmnt','SumdistM','MeanDayElev','DayElevRange'),]

### collapse to ByIbexSeason Dataframe ####
library(tidyverse)
#create a data frame for each ibex by season
ByIbexSeason= DailyData %>% group_by(id,Season) %>%
     summarise(Sex= head(Sex,1),
               DaysInSeason=length(DailyN_GPSfix),#how many days in this season
               DailyN_GPSfix=mean(DailyN_GPSfix,na.rm=T),
               MxDailyDisplcmnt=mean(MxDailyDisplcmnt,na.rm=T),
               SumdistM = mean(SumdistM,na.rm=T),
               NetDailyDisplcmnt = mean(NetDailyDisplcmnt,na.rm=T),
               TortuosityuntilMax = mean(TortuosityuntilMax,na.rm=T),
               DayElevRange = mean(DayElevRange,na.rm=T),
               MeanDayElev = mean(MeanDayElev,na.rm=T),
               CV_SumdistM=sd(SumdistM,na.rm=T)/mean(SumdistM,na.rm=T)
               )
ByIbexSeason=as.data.frame(ByIbexSeason)
View(ByIbexSeason)

#mean(df$val[df$sta=="A"])
#just chekcing, looks good!
#length(which(DailyData$id=='Ofra'&DailyData$Season=='Rut'))
#mean(DailyData$MxDailyDisplcmnt[which(DailyData$id=='Ofra'&DailyData$Season=='Rut')])

### export the datafrom to XLSX ####
writexl::write_xlsx(x=ByIbex2,   path='DataExportByIbex2.xlsx') #,sheetName = "ByIbex2"
writexl::write_xlsx(x=DailyData, path='DataExportDailyData.xlsx') #,sheetName = "DailyData"
writexl::write_xlsx(x=ByIbexSeason, path='DataExportByIbexSeason.xlsx') #,sheetName = "DailyData"

### plotting histograms of the Daily movements ####
hist(DailyData$DailyN_GPSfix,breaks=50)
hist(DailyData$SumdistM,breaks=50)
hist(DailyData$MxDailyDisplcmnt,breaks=50)
hist(DailyData$TortuosityAllDay,breaks=50)
hist(DailyData$TortuosityAllDay[DailyData$TortuosityAllDay<50],breaks=50)#not very reliable, to ignore this one
hist(DailyData$TortuosityAllDay[DailyData$TortuosityAllDay<50 & DailyData$SumdistM>1000],breaks=50)
hist(DailyData$TortuosityuntilMax,breaks=50)
hist(DailyData$TortuosityuntilMax[DailyData$SumdistM>1000],breaks=50)
hist(DailyData$MeanDayElev[DailyData$SumdistM>1],breaks=50)
hist(DailyData$DayElevRange[DailyData$SumdistM>1],breaks=50)


### plotting tracks:  #####
#plot(GPSdat3or4h_ltrj);plot all animals together
#plot(GPSdat3or4h_ltrj,id='Stern');plot one  animals together
plot(GPSdatFltr_ltrj[5]);#plot track by burst, not usefull
#hist(GPSdat3or4h_ltrj[1], "UTC_Time", freq = TRUE)
head(GPSdat3or4h_ltrj[[5]])

#there is a problem here with the id() function that is now stopped so need to fix this bit ERROR
for (IbxCntr in 1:length(unique(id(GPSdatFltr_ltrj)))){
  CurrName=unique(id(GPSdatFltr_ltrj))[IbxCntr]
  #indx=which(id(GPSdat3or4h_ltrj)==unique(id(GPSdat3or4h_ltrj))[IbxCntr])
  par(pty="s")
  plot(GPSdatFltr_ltrj,id=CurrName,main=CurrName, xlab= c("Easting"),ylab=("Northing")) 
  
  #infolocs(GPSdat3or4h_ltrj)@AnimalID
} 
rm(CurrName,IbxCntr,Indx,SamplingIntervals)
#t=which.ltraj(GPSdatFltr_ltrj,"id='Alexander'")
#infolocs(GPSdat3or4h_ltrj, id)
#trajdyn(GPSdatFltr_ltrj[4])

## Saving processed file first time, before HR analysis ####
Name='IbexMoveData1.rdata'
save(list=ls(),file=Name)
print('saved  GPS data as ltraj  ')


### t-test and simple box plots ####
#### males or females closer to the road? ##

# Welch t-test- for road distance
t.test(NofFixes ~ Sex, data=ByIbex2)
boxplot(NofFixes ~ Sex,data=ByIbex2, main="GPS fixes", xlab="Sex", ylab="NofFixes")

#filter non relevant ibex
ByIbex2fltrd=ByIbex2[!ByIbex2$AnimalID %in% c('Rifka','Masada','Stern'),]

t.test(FixesNearRoad ~ Sex, data=ByIbex2fltrd)
boxplot(FixesNearRoad ~ Sex,data=ByIbex2fltrd, main="N point near road",xlab="Sex", ylab="FixesNearRoad")

t.test(MeanDistToRoad ~ Sex, data=ByIbex2fltrd)
boxplot(MeanDistToRoad ~ Sex,data=ByIbex2fltrd, main="MeanDistToRoad", xlab="Sex", ylab="MeanDistToRoad")

t.test(PropNearRoad ~ Sex, data=ByIbex2fltrd)
boxplot(PropNearRoad ~ Sex,data=ByIbex2fltrd, main="PropNearRoads",  xlab="Sex", ylab="PropNearRoad")


# testing male-female difference in movement with simple ttest
t.test(SumdistM~Sex,data=ByIbex2)
boxplot(SumdistM ~ Sex,data=ByIbex2, main="Travel dist", xlab="Sex", ylab="SumdistM")

t.test(log(SumdistM)~Sex,data=ByIbex2)### significant differences
boxplot(log(SumdistM) ~ Sex,data=ByIbex2, main="Travel dist", xlab="Sex", ylab="log_SumdistM")

t.test(NetDailyDisplcmnt~Sex,data=ByIbex2)
boxplot(NetDailyDisplcmnt ~ Sex,data=ByIbex2, main="NetDailyDisplcmnt", xlab="Sex", ylab="Net Daily Displcmnt (m)")

t.test(MxDailyDisplcmnt~Sex,data=ByIbex2)
boxplot(MxDailyDisplcmnt ~ Sex,data=ByIbex2, main="MxDailyDisplcmnt", xlab="Sex", ylab="Max Daily Displcmnt (m)")

t.test(TortuosityAllDay~Sex,data=ByIbex2)
boxplot(TortuosityAllDay ~ Sex,data=ByIbex2, main="TortuosityAllDay",  xlab="Sex", ylab="TortuosityAllDay")

t.test(TortuosityuntilMax~Sex,data=ByIbex2)
boxplot(TortuosityuntilMax ~ Sex,data=ByIbex2, main="TortuosityuntilMax",  xlab="Sex", ylab="TortuosityuntilMax")

t.test(MeanDayElev~Sex,data=ByIbex2)### significant differences
boxplot(MeanDayElev ~ Sex,data=ByIbex2, main="MeanDayElev",  xlab="Sex", ylab="MeanDayElev (m)")

t.test(DayElevRange  ~ Sex,data=ByIbex2)### Almost significant differences
boxplot(DayElevRange ~ Sex,data=ByIbex2, main="DayElevRange",  xlab="Sex", ylab="DayElevRange (m)")

### box plots for DailyData, not IbexData #### 
boxplot(SumdistM ~ Sex,data=DailyData, main="Sex effect on movement", xlab="Sex", ylab="SumdistM")
boxplot(NetDailyDisplcmnt ~ Sex,data=DailyData, main="Sex effect on movement", xlab="Sex", ylab="NetDailyDisplcmnt")
boxplot(MxDailyDisplcmnt ~ Sex,data=DailyData, main="Sex effect on movement", xlab="Sex", ylab="MxDailyDisplcmnt")
boxplot(TortuosityuntilMax ~ Sex,data=DailyData, main="Sex effect on movement", xlab="Sex", ylab="TortuosityuntilMax")
boxplot(DayElevRange ~ Sex,data=DailyData, main="Sex effect on daly Elev changes", xlab="Sex", ylab="DayElevRange")
boxplot(MeanDayElev ~ Sex,data=DailyData, main="Sex effect on mean Elev", xlab="Sex", ylab="MeanDayElev")

### boxplot with ggplot for sex and season ####
ggplot(DailyData, aes(x=Season, y=SumdistM,           fill=Sex)) + geom_boxplot(position=position_dodge(1))
ggplot(DailyData, aes(x=Season, y=NetDailyDisplcmnt,  fill=Sex)) + geom_boxplot(position=position_dodge(1))
ggplot(DailyData, aes(x=Season, y=MxDailyDisplcmnt,   fill=Sex)) + geom_boxplot(position=position_dodge(1))
ggplot(DailyData, aes(x=Season, y=TortuosityuntilMax, fill=Sex)) + geom_boxplot(position=position_dodge(1))
ggplot(DailyData, aes(x=Season, y=DayElevRange,       fill=Sex)) + geom_boxplot(position=position_dodge(1))
ggplot(DailyData, aes(x=Season, y=MeanDayElev,        fill=Sex)) + geom_boxplot(position=position_dodge(1))
##same plots with a limit on Y axis can be:
ggplot(DailyData, aes(x=Season, y=SumdistM,           fill=Sex)) + geom_boxplot(position=position_dodge(1))+ylim(0,5000)

### stat models 1 - LMMs for sex only ####
MixModel_SumdistM1           =lmer(SumdistM          ~ 1 + Sex  + (1|id) ,data=DailyData, verbose=F,REML=T)#
MixModel_NetDailyDisplcmnt1  =lmer(NetDailyDisplcmnt ~ 1 + Sex  + (1|id) ,data=DailyData, verbose=F,REML=T)#
MixModel_MxDailyDisplcmnt1   =lmer(MxDailyDisplcmnt  ~ 1 + Sex  + (1|id) ,data=DailyData, verbose=F,REML=T)#
MixModel_TortuosityuntilMax1 =lmer(TortuosityuntilMax~ 1 + Sex  + (1|id) ,data=DailyData, verbose=F,REML=T)#
MixModel_MeanDayElev         =lmer(MeanDayElev       ~ 1 + Sex  + (1|id) ,data=DailyData, verbose=F,REML=T)#
MixModel_DayElevRange        =lmer(DayElevRange      ~ 1 + Sex  + (1|id) ,data=DailyData, verbose=F,REML=T)#

#for exploring model results
summary(MixModel_MeanDayElev)
plot(residuals(MixModel_DayElevRange))
summary(MixModel_MeanDayElev)
require(lmerTest)#https://rdrr.io/cran/lme4/man/pvalues.html
  anova(MixModel_MeanDayElev,ddf = "Kenward-Roger")
  #car::Anova(MixModel_MeanDayElev)
  library(nlme)
  MixModel2_MeanDayElev=lme(MeanDayElev       ~ 1 + Sex, random=~ 1|id ,data=DailyData)#
  anova(MixModel2_MeanDayElev) 
#  Model.pval<-pvals.fnc(Model, nsim = n, withMCMC = TRUE)#needs coda and languageR packages
  
### stat models 2 season*Sex interaction LMM ##### 
MixModelwInrtc_NetDailyDisplcmnt1  =lmer(NetDailyDisplcmnt ~ 1 + Sex*Season + Sex + Season + (1|id) ,data=DailyData, verbose=F,REML=F)#
MixModelwInrtc_SumdistM1           =lmer(SumdistM          ~ 1 + Sex*Season + Sex + Season + (1|id) ,data=DailyData, verbose=F,REML=F)#
MixModelwInrtc_MxDailyDisplcmnt1   =lmer(MxDailyDisplcmnt  ~ 1 + Sex*Season + Sex + Season + (1|id) ,data=DailyData, verbose=F,REML=F)#
MixModelwInrtc_TortuosityuntilMax1 =lmer(TortuosityuntilMax~ 1 + Sex*Season + Sex + Season + (1|id) ,data=DailyData, verbose=F,REML=F)#
MixModelwInrtc_MeanDayElev         =lmer(MeanDayElev       ~ 1 + Sex*Season + Sex + Season + (1|id) ,data=DailyData, verbose=F,REML=F)#
MixModelwInrtc_DayElevRange        =lmer(DayElevRange      ~ 1 + Sex*Season + Sex + Season + (1|id) ,data=DailyData, verbose=F,REML=F)#

#for exploring model results
summary(MixModelwInrtc_NetDailyDisplcmnt1)
plot(residuals(MixModelwInrtc_NetDailyDisplcmnt1))

#plots of model residuals
plot(MixModel_NetDailyDisplcmnt1)
plot(MixModel_SumdistM1)
plot(MixModel_MxDailyDisplcmnt1)
plot(MixModel_TortuosityuntilMax1)


### stat models 3 model comparison code! #####

#this is a helper function that will compare 5 models for a chosen parameter in the DailyData data.frame. its being called later
LMMmodeler1<-function(DataF=DailyData,YvarInd=NA,YvarName=NA){
  require(AICcmodavg)
  ## mixed models that account for the pseudo replocation in plots ###
  MixModelN1 =lmer(DataF[,YvarInd] ~ 1 +                               (1|id) ,data=DataF, verbose=F,REML=F)#
  MixModel2  =lmer(DataF[,YvarInd] ~ 1 +        Season               + (1|id) ,data=DataF, verbose=F,REML=F)#,
  MixModel3  =lmer(DataF[,YvarInd] ~ 1 +  Sex                        + (1|id) ,data=DataF, verbose=F,REML=F)#,
  MixModel4  =lmer(DataF[,YvarInd] ~ 1 +  Sex + Season               + (1|id) ,data=DataF, verbose=F,REML=F)#,
  MixModel5  =lmer(DataF[,YvarInd] ~ 1 +  Sex + Season + Sex*Season  + (1|id) ,data=DataF, verbose=F,REML=F)#,
  
  ## model comparison ###
  #comparing  models -
  ModelNames=c('MixModelN1','MixModel2_Season','MixModel3_Sex','MixModel4_Sex_Season','MixModel5_SexSeasonIntr')
  ModelsList=list(MixModelN1,   MixModel2,      MixModel3,             MixModel4,            MixModel5)
  TableAicMix=aictab(cand.set=ModelsList,modnames=ModelNames,weights = TRUE)
  BestModel=    (ModelsList[[as.numeric(rownames(TableAicMix)[1])]]);
  BestModelName=(ModelNames[[as.numeric(rownames(TableAicMix)[1])]]);
  ScndBestModel=(ModelsList[[as.numeric(rownames(TableAicMix)[2])]]);
  
  #print(YvarName)
  #print(TableAicSimple)
  print(paste('the AIC model comparison table ',YvarName,'is'))
  print(TableAicMix)
  print(paste('and the best model for',YvarName, 'was:',BestModelName,' - here is the summary:'))
  print(summary(BestModel))
  
  return(list(ModelsList,ModelNames,TableAicMix,BestModel,ScndBestModel) )
}

#here I call the helper function, tells her what to look at. change to the specific parameter
names(DailyData)
ParameterToModel='NetDailyDisplcmnt'#Choose the name her for running the models on the chosen parameter
CurrModels=LMMmodeler1(DataF=DailyData,YvarName=ParameterToModel,YvarInd=which(names(DailyData)==ParameterToModel)) #runnig it wiout one outlier: point 43!BiosolYA2[-43,]

ParameterToModel='TortuosityuntilMax'#Choose the name her for running the models on the chosen parameter
CurrModels=LMMmodeler1(DataF=DailyData,YvarName=ParameterToModel,YvarInd=which(names(DailyData)==ParameterToModel)) #runnig it wiout one outlier: point 43!BiosolYA2[-43,]

ParameterToModel='SumdistM'#Choose the name her for running the models on the chosen parameter
CurrModels=LMMmodeler1(DataF=DailyData,YvarName=ParameterToModel,YvarInd=which(names(DailyData)==ParameterToModel)) #runnig it wiout one outlier: point 43!BiosolYA2[-43,]

ParameterToModel='DayElevRange'#Choose the name her for running the models on the chosen parameter
CurrModels=LMMmodeler1(DataF=DailyData,YvarName=ParameterToModel,YvarInd=which(names(DailyData)==ParameterToModel)) #runnig it wiout one outlier: point 43!BiosolYA2[-43,]

ParameterToModel='MeanDayElev'#Choose the name her for running the models on the chosen parameter
CurrModels=LMMmodeler1(DataF=DailyData,YvarName=ParameterToModel,YvarInd=which(names(DailyData)==ParameterToModel)) #runnig it wiout one outlier: point 43!BiosolYA2[-43,]

# #getting the AIC table and Models
# ModelsList=CurrModels[[1]];ModelNames=CurrModels[[2]];TableAicGLM=CurrModels[[3]];rm(CurrModels)
# BestModel       =(ModelsList[[as.numeric(rownames(TableAicGLM)[1])]]);
# BestModelName   =(ModelNames[[as.numeric(rownames(TableAicGLM)[1])]]);
# ScndBestModel=(ModelsList[[as.numeric(rownames(TableAicGLM)[2])]]);
# ScndBestModelName   =(ModelNames[[as.numeric(rownames(TableAicGLM)[2])]]);
# print(ParameterToModel);print(TableAicGLM)
# summary(BestModel)



### stat models 4 model comparison code for no Id, just simple Linear Model #####

#this is a helper function that will compare 5 models for a chosen parameter in the DailyData data.frame. its being called later

#here I call the helper function, tells her what to look at. change to the specific parameter
names(ByIbexSeason)
ParameterToModel='NetDailyDisplcmnt'#Choose the name her for running the models on the chosen parameter
CurrModels=LMMmodeler1(DataF=ByIbexSeason,YvarName=ParameterToModel,YvarInd=which(names(ByIbexSeason)==ParameterToModel)) #runnig it wiout one outlier: point 43!BiosolYA2[-43,]

ParameterToModel='SumdistM'#Choose the name her for running the models on the chosen parameter
CurrModels=LMMmodeler1(DataF=ByIbexSeason,YvarName=ParameterToModel,YvarInd=which(names(ByIbexSeason)==ParameterToModel)) #runnig it wiout one outlier: point 43!BiosolYA2[-43,]

ParameterToModel='TortuosityuntilMax'#Choose the name her for running the models on the chosen parameter
CurrModels=LMMmodeler1(DataF=ByIbexSeason,YvarName=ParameterToModel,YvarInd=which(names(ByIbexSeason)==ParameterToModel)) #runnig it wiout one outlier: point 43!BiosolYA2[-43,]

ParameterToModel='SumdistM'#Choose the name her for running the models on the chosen parameter
CurrModels=LMMmodeler1(DataF=ByIbexSeason,YvarName=ParameterToModel,YvarInd=which(names(ByIbexSeason)==ParameterToModel)) #runnig it wiout one outlier: point 43!BiosolYA2[-43,]

ParameterToModel='DayElevRange'#Choose the name her for running the models on the chosen parameter
CurrModels=LMMmodeler1(DataF=ByIbexSeason,YvarName=ParameterToModel,YvarInd=which(names(ByIbexSeason)==ParameterToModel)) #runnig it wiout one outlier: point 43!BiosolYA2[-43,]

ParameterToModel='MeanDayElev'#Choose the name her for running the models on the chosen parameter
CurrModels=LMMmodeler1(DataF=ByIbexSeason,YvarName=ParameterToModel,YvarInd=which(names(ByIbexSeason)==ParameterToModel)) #runnig it wiout one outlier: point 43!BiosolYA2[-43,]

############################

### DRAFTS #####
## create predictions from model
ForPredict=data.frame(Season=factor(c('Dry', 'Rut', 'Wet','Dry', 'Rut', 'Wet')), 
                      Sex=as.factor(c('M','M','M','F','F','F')),
                      id=as.factor('Ada','Ada','Ada','Ada','Ada','Ada'),
                      row.names =NULL);
PredBstM=predict(BestModel, newdata = ForPredict,type = c( "response"),se.fit=T)


QuantilesToPred=seq(from=0.00,to=1,by=0.1)#c(0.05,0.25,0.5,0.75,0.95);
Quantiles=quantile( x=DataF$VFAcnt,probs=seq(0,1,0.1));
ArrayPred=rbind(as.numeric(rep(Quantiles,9)),rep(levels(DataF$Sample),times=3),rep(levels(DataF$Depth),each=3))
ForPredict=data.frame(t(ArrayPred),row.names =NULL)#,rep(c(0),times=dim(ArrayPred)[2]),rep(c(0),times=dim(ArrayPred)[2]),row.names =NULL)#
names(ForPredict)=c('VFAcnt','Sample','Depth');
ForPredict$Sample=factor(ForPredict$Sample);
ForPredict$Depth=factor(ForPredict$Depth,levels=c('High','Mid','Low'));
ForPredict$VFAcnt=as.numeric(as.character(ForPredict$VFAcnt));
#str(ForPredict)

PredBstM=predict(BstMdl, newdata = ForPredict,type = c( "response"))
PredSecBstM=(predict(ScndBstMdl, newdata = ForPredict,type = c( "response")))
Pred_df=data.frame(ForPredict,PredBstM,PredSecBstM,QuantilesToPred);#names(Pred_df)=c('PredVal','SegCerNum','MovePersonQunatl')
Pred_df$PredC=sapply( 1:length(PredBstM), function(i){weighted.mean(c(Pred_df$PredBstM[i],Pred_df$PredSecBstM[i]),w=TableAic$AICcWt[1:2])})#combine predictions based on 2 top model's weights
#with(Pred_df,plot(x=VFAcnt,y=Pred2))
return(Pred_df)

justLMmodel<-function(DataF=ByIbexSeason,YvarInd=NA,YvarName=NA){
  ## mixed models that account for the pseudo replocation in plots ###
  MixModelN1 =lm(DataF[,YvarInd] ~ 1                                 ,data=DataF)#
  MixModel2  =lm(DataF[,YvarInd] ~ 1 +           Season              ,data=DataF) #,
  MixModel3  =lm(DataF[,YvarInd] ~ 1 +  Sex                          ,data=DataF)#,
  MixModel4  =lm(DataF[,YvarInd] ~ 1 +  Sex + Season                 ,data=DataF)#,
  MixModel5  =lm(DataF[,YvarInd] ~ 1 +  Sex + Season + Sex*Season    ,data=DataF)#,
  
  ## model comparison ###
  #comparing  models - As above. The best model is MixModel2_VFA, followed by MixModel4_VFA_Depth!
  ModelNames=c('MixModelN1','MixModel2_Season','MixModel3_Sex','MixModel4_Sex_Season','MixModel5_SexSeasonIntr')
  ModelsList=list(MixModelN1,   MixModel2,      MixModel3,             MixModel4,            MixModel5)
  TableAicMix=aictab(cand.set=ModelsList,modnames=ModelNames,weights = TRUE)
  BestModel=    (ModelsList[[as.numeric(rownames(TableAicMix)[1])]]);
  BestModelName=(ModelNames[[as.numeric(rownames(TableAicMix)[1])]]);
  ScndBestModel=(ModelsList[[as.numeric(rownames(TableAicMix)[2])]]);
  
  #print(YvarName)
  #print(TableAicSimple)
  print(paste('the AIC model comparison table ',YvarName,'is'))
  print(TableAicMix)
  print(paste('and the best model for',YvarName, 'was:',BestModelName,' - here is the summary:'))
  print(summary(BestModel))
  
  return(list(ModelsList,ModelNames,TableAicMix,BestModel,ScndBestModel) )
}
## trying to calculate HR with Locoh (Getz et al 2007 ) using adehabitat ####
class(GPSdat3or4h)
proj4string(GPSdat3or4h) <- CRS(Proj) #NEED TO check cooridnate system
#head(GPSdat3or4h)
#GPSdat3or4h_sp=ltraj2spdf(GPSdat3or4h_ltrj)
#head(as.data.frame(GPSdat3or4h))

LoCoHkHR=LoCoH.k(xy=GPSdat3or4h[,3], k=5, unin = "m",  unout =  "km2",        duplicates="random", amount = NULL)
LoCoHkHR=LoCoH.k(xy=GPSdat1h[,3], k=5, unin = "m",  unout =  "km2",        duplicates="random", amount = NULL)

plot(clusthr(GPSdat1h))





xysp <- SpatialPoints(GPSdat3or4hDF[,c("Easting","Northing")])
proj4string(xysp) <- CRS("+proj=utm +zone=36 +ellps=WGS84") #NEED TO check cooridnate system
head(as.data.frame(xysp))


#just to work on a small file with bursts (days) from 10 tags
#LizXYdataSample=LizXYdataSubset[id=id(LizXYdataSubset[10])]


as.ltraj(xy = tt[,c("Longitude","Latitude")], 
         date = tt$UTC_Date.Time, 
         id =factor( tt$AnimalID),
         burst = tt$DaysAsBurst,
         typeII=TRUE,proj4string=
           GPSdatFltr);

print('now converting to track ltraj ')
#conversion to class of adehabitatLT
GPSdat1h=as.ltraj(xy = GPSdat1h[,c("Easting","Northing")], 
                       date = GPSdat1h$UTC_Date.Time, 
                       id = GPSdat1h$AnimalID,
                       burst = GPSdat1h$DaysAsBurst,
                       typeII=TRUE,
                       infolocs =GPSdat1h);

GPSdatFltr=as.ltraj(xy = GPSdatFltr[,c("Longitude","Latitude")], 
                      date = GPSdatFltr$UTC_Date.Time, 
                      id =factor( GPSdatFltr$AnimalID),
                      burst = GPSdatFltr$DaysAsBurst,
                      typeII=TRUE,
                      infolocs =GPSdatFltr);



GPSdat3or4h=as.ltraj(xy = GPSdat3or4h[,c("Easting","Northing")], 
                      date = GPSdat3or4h$UTC_Date.Time, 
                      id = GPSdat3or4h$AnimalID,
                      burst = GPSdat3or4h$DaysAsBurst,
                      typeII=TRUE,
                      infolocs =GPSdat3or4h);

tt=subset(GPSdat3or4h, AnimalID=="Pola");tt=droplevels(tt);str(tt)
as.ltraj(xy = tt, 
         date = tt$UTC_Date.Time, 
         id = tt$AnimalID,
         burst = tt$DaysAsBurst,
         typeII=TRUE);

coordinates(tt) = c("Easting", "Northing") # specify column names


#in addiotn the infolocs include the distance and R2N and dx dy and rel.angle turnign angle
head(LizXYdataSubset[[2]])
#just to work on a small file with bursts (days) from 10 lizards
LizXYdataSample=LizXYdataSubset[id=id(LizXYdataSubset[10])]
#LizXYdataSample=LizXYdataFull  [id=id(LizXYdataFull)[1:10]]

GPSdatFltrdf=ld(GPSdatFltr)
for (BurCnt in 1:levels(tt$burst)){
  
}

     
       