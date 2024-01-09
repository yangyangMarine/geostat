# sensitivity analysis for Single target detection algorithm on shallow water broadband survey
# dataset: Igombe 2019
# R version: 4.0.5
# raw data exported from: Echoview 13.1

# loading
  pacman::p_load(ggplot2,dplyr,PerformanceAnalytics, zoo, tidyverse, RColorBrewer)

# automation export STD parameters
  library(RDCOMClient)
  list <- c('Igombe_FM', 'Ihale_FM','Kahumulo_FM','Kayenze_FM','Kome_Nyalubungo_FM','Kome_FM')
  timestart<-Sys.time()
  dir <- paste0("I:/MyPassport/shallow water 2020/",list[i])

  EVApp <- COMCreate('EchoviewCom.EvApplication')
  EVFile <- EVApp$OpenFile("I:\\MyPassport\\shallow water data process\\Vic processing_-70 threshold.EV")
  fileset1<-EVFile[["Filesets"]]$FindByName('Fileset1')
  allfiles<-list.files(dir, pattern='.raw', recursive=T, full.names=T)
  for(n in allfiles){fileset1[["DataFiles"]]$Add(n)}
  
  objVariable <- EVFile[['Variables']]$FindByName("FMSTD")
  EVVar_propSTD<- objVariable[['Properties']][['SingleTargetDetectionWidebandParameters']]
  EVVar_propSTD_CW<- objVariable[['Properties']][['SingleTargetDetectionParameters']]

  #  TS threshold
  EVVar_propSTD_CW[['TsThreshold']] <- -70
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("TSthreshold70",".csv")))
  print(EVApp$GetLastLogMessage())
  
  EVVar_propSTD_CW[['TsThreshold']] <- -65
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("TSthreshold65",".csv")))
  print(EVApp$GetLastLogMessage())
  
  EVVar_propSTD_CW[['TsThreshold']] <- -60
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("TSthreshold60",".csv")))
  print(EVApp$GetLastLogMessage())
  
  # beam compensation
  EVVar_propSTD_CW[['TsThreshold']] <- -65
  EVVar_propSTD[['MaximumBeamCompensation']] <- 4
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("BeamCom4",".csv")))
  print(EVApp$GetLastLogMessage())
  
  EVVar_propSTD[['MaximumBeamCompensation']] <- 8
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("BeamCom8",".csv")))
  print(EVApp$GetLastLogMessage())
  
  EVVar_propSTD[['MaximumBeamCompensation']] <- 12
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("BeamCom12",".csv")))
  print(EVApp$GetLastLogMessage())
  
  # StdDevOfMajorAxisAngles
  EVVar_propSTD[['MaximumBeamCompensation']] <- 6
  EVVar_propSTD[['MaximumStdDevOfMajorAxisAngles']] <- 0.6
  EVVar_propSTD[['MaximumStdDevOfMinorAxisAngles']] <- 0.6
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("AngleStd0.6",".csv")))
  print(EVApp$GetLastLogMessage())
  
  EVVar_propSTD[['MaximumStdDevOfMajorAxisAngles']] <- 1.2
  EVVar_propSTD[['MaximumStdDevOfMinorAxisAngles']] <- 1.2
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("AngleStd1.2",".csv")))
  print(EVApp$GetLastLogMessage())
  
  # MinimumPulseLength
  EVVar_propSTD[['MaximumStdDevOfMajorAxisAngles']] <- 3
  EVVar_propSTD[['MaximumStdDevOfMinorAxisAngles']] <- 3
  EVVar_propSTD_CW[['MinimumPulseLength']] <- 0.8
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("PulseLength0.8_1.5",".csv")))
  print(EVApp$GetLastLogMessage())
  
  EVVar_propSTD_CW[['MaximumPulseLength']] <- 1.3
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("PulseLength0.8_1.3",".csv")))
  print(EVApp$GetLastLogMessage())
  
  EVVar_propSTD_CW[['MinimumPulseLength']] <- 0.5
  EVVar_propSTD_CW[['MaximumPulseLength']] <- 1.3
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("PulseLength0.5_1.3",".csv")))
  print(EVApp$GetLastLogMessage())
  
  # target seperation
  EVVar_propSTD_CW[['MinimumPulseLength']] <- 0.5
  EVVar_propSTD_CW[['MaximumPulseLength']] <- 1.5
  EVVar_propSTD[['MinimumTargetSeparation']] <- 0.1
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("Seperation0.1",".csv")))
  print(EVApp$GetLastLogMessage())

  EVVar_propSTD[['MinimumTargetSeperation']] <- 0.2
  objVariable$ExportSingleTargetsByCellsAll(file.path("I:\\MyPassport\\Shallow water data process\\sensitivity analysis 2023\\FM\\",paste0("Seperation0.2",".csv")))
  print(EVApp$GetLastLogMessage())

  EVApp$Quit()
  timeend<-Sys.time()
  runningtime<-timeend-timestart
  print(runningtime)

  