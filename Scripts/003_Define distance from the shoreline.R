# under R version 4.0.5 
  rm(list=ls())
  library(sf)
  library(raster)
  library(tidyverse)
  library(RColorBrewer)
  library(sfheaders) 
  library(rgdal)
  #library(GDAL) # NEW version of rgdal
  library(sp)
  library(dplyr)
  library(stars)

###########################################################################
########## find regions within certain distance from the shore  ###########
###########################################################################
  
  # import lake.shp and survey.csv  
  lake <- st_read("C:\\Users\\yyy1\\Downloads\\LV bathy\\LakeVictoriaShoreline.shp")
  lake <- st_read("C:\\Users\\yangy\\OneDrive - University of St Andrews\\PhD papers\\chapter1_shallow margins\\LV bathy\\LakeVictoriaShoreline.shp")

  # select only tanzania section
  lake <- st_crop(lake, xmin= 683002.5, xmax= 1021439, ymin= -316155.5, ymax= -104531.4)  # equal to N = -1
  
  # make vacant cells with certain size
  grid <- st_make_grid(lake, cellsize = 5000, what = "polygons")  # polygons, corners, centers
  
  # lake coast
  lakecoast <- st_cast(lake, "MULTILINESTRING")
  
  # crop grid with lake polygons
  grid <- st_intersection(lake, grid)
  paste("the whole area of the lake is", sum(st_area(grid)), "m^2")  # lake[[3]]
  
  # distance from each small grid to coastline
  dist <- st_distance(lakecoast, grid)
  df_dist <- as.numeric(t(dist))
  
  # Area of cells within certain distance from the shore
  Area_within_d <- sum(df_dist <= 2000)/length(t(dist)) * sum(st_area(grid))/10^6 
 
  # rasterise grid and mapview
  grid$dist <- t(dist)
  grid.st = st_rasterize(grid["dist"],)
  plot(grid.st)   
  
  ggplot() +  
    geom_stars(data = grid.st) +
    coord_equal() +
    theme_void() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) 

  mapview::mapview(grid.st) 
  beepr::beep(8)
  
################################################################
#############  get distance from the shore R 3.5.0 #############
################################################################
  
  ST_dir <- ("D:\\MyPassport\\exported_v3\\ST")
  setwd(ST_dir)
  ST_listfile <- list.files(ST_dir, pattern = '_FM')
  ST_listfile = ST_listfile[66:130]
  
  for(i in 1:65){
        df_raw <- read.csv(ST_listfile[i], header = T) %>% 
          filter( Lat_M != '999', Lon_M != '999')
        str(df_raw)
        df_raw <- df_raw %>% dplyr::select("Lon_M","Lat_M","Num_targets","Beam_volume_sum","Height_mean")

        pnts_sf <- st_as_sf(df_raw, crs=4326L, coords = c("Lon_M", "Lat_M"))
        
        lake_pnt <- lake %>%
          st_cast("POLYGON")  %>%
          st_cast("MULTILINESTRING") %>%
          st_cast("LINESTRING") %>%
          st_cast("MULTIPOINT") %>%
          st_cast("POINT")
      
        lake_pnt <- st_transform(lake_pnt, crs = 4326L)  # transform lake shp CRS into 4326
        
        nearest = st_nearest_feature(pnts_sf,lake_pnt)  # only index
        dist = st_distance(pnts_sf, lake_pnt[nearest,], by_element=TRUE) # pairwise = TRUE
        df_raw$dist <- dist
        write.csv(df_raw, file = paste0("Dist_",ST_listfile[i]))
      }
  