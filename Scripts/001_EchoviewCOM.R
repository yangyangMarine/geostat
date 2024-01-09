###############################################################################
############### batch export data Echoview 13.1
############### R 3.5.1 
###############################################################################
  rm(list=ls())
  library(RDCOMClient)
  library(stringr)
  listdir <- list.dirs('D:/MyPassport', full.names = T, recursive =F) 
  EVApp <- COMCreate('EchoviewCom.EvApplication')
  EVFile <- EVApp$OpenFile("D:\\MyPassport\\Shallow water data process\\Vic processing_-70 threshold.EV")

 # loop over 5 round of surveys
 for(j in 1:5){
   dir=(listdir[j])
   setwd(dir)
   listdir_r <- list.dirs(dir, full.names = T, recursive =F) 
   FM = grep('_FM', listdir_r) # FM data only 

   listdir_r <- listdir_r[FM]

      for(i in listdir_r){
        timestart<-Sys.time()
        fileset1<-EVFile[["Filesets"]]$FindByName('Fileset1')
        setwd(i)

        allfiles<-list.files(i, pattern='.raw', recursive=T, full.names=T)
        for(n in allfiles){fileset1[["DataFiles"]]$Add(n)}
        
        objVariable <- EVFile[['Variables']]$FindByName("FMSTD")

        name <- as.data.frame(stringr::str_split(i, "/"))[4,1] 
        objVariable$ExportSingleTargetsByCellsAll(file.path("D:/MyPassport/exported_v3/ST/",paste0('r',7,'_',name,'_FMST.csv')))
        
        PingNo <- objVariable$MeasurementCount()
        #objVariable$DetectFishTracks("fish", 1, PingNo)
        
        #objVariable$ExportFishTracksByRegionsAll(file.path("I:/MyPassport/exported_v3/tracks/",paste0('r',j,'_',name,'_FMtracks.csv')))
        
        fileset1[["DataFiles"]]$RemoveAll()
        timeend<-Sys.time()
        runningtime<-timeend-timestart
        print(runningtime)
      }
  EVApp$Quit()

###############################
########### CW data
###############################
  EVApp <- COMCreate('EchoviewCom.EvApplication')
  EVFile <- EVApp$OpenFile("I:\\MyPassport\\Shallow water data process\\Vic processing_-70 threshold.EV")
  for(j in 1:5){
    dir=(listdir[j])
    setwd(dir)

    listdir_r <- list.dirs(dir, full.names = T, recursive =F) 
    CW = grep('_CW', listdir_r) # only cw folder 
    listdir_r <- listdir_r[CW]
   
    for(i in listdir_r){
      timestart<-Sys.time()
  
      fileset1<-EVFile[["Filesets"]]$FindByName('Fileset1')
    
      allfiles<-list.files(i, pattern='.raw', recursive=T, full.names=T)
      for(n in allfiles){fileset1[["DataFiles"]]$Add(n)}
    
      objVariable <- EVFile[['Variables']]$FindByName("CWSTD")
    
      name <- as.data.frame(stringr::str_split(i, "/"))[4,1] # extract site name from i
      objVariable$ExportSingleTargetsByCellsAll(file.path("D:/MyPassport/exported_v3/ST/",paste0('r',j,'_',name,'_CWST.csv')))
    
      PingNo <- objVariable$MeasurementCount()
      objVariable$DetectFishTracks("fish", 1, PingNo)
      objVariable$ExportFishTracksByRegionsAll(file.path("I:/MyPassport/exported_v3/tracks/",paste0('r',j,'_',name,'_CWtracks.csv')))
    
    nasc <- EVFile[['Variables']]$FindByName("CWNASC")
    nasc$ExportIntegrationByCellsAll(file.path("D:/MyPassport/exported_v3/bathy/",paste0('r',j,'_',name,'_NASC.csv')))
    
    fileset1[["DataFiles"]]$RemoveAll()
    timeend<-Sys.time()
    runningtime<-timeend-timestart
    print(runningtime)
  }
}
  EVApp$Quit()

