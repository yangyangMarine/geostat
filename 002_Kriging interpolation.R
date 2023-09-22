##############################################################################
########   Kriging interpolation for abundance estimate, R 4.0.3 
########
# to performing ordinary kriging serveral steps should be applied before modelling
# 1. Check normality, if not, apply log transform/square root tran
# 2. top cut to remove abnormal values
# 3. anisotropic check, Ominidirectional variogram
# 4. remove the trend of the data, trend = non-stationary 1st order 2nd order 
# 5. cross validation
# 5. calculate the mean and global variance, Paul, F. 1999
#############################################################################

  rm(list=ls())
  # load package
  pacman::p_load(gstat,
               sf,
               sp,
               readxl,
               mapview, # interactive map
               dplyr,
               automap, # auto kriging 
               ggpubr,
               splancs, # hand draw polygons V2.01-43
               raster,# rasterizing 
               RandomFieldsUtils, # dependency of geoR package V0.3.20
               RandomFields, # dependency of geoR package V3.0.62
               geoR  # to remove the 1st/2nd trend of the raw data V1.7-2
) 

# Save results to new
    result <- matrix(NA, 145, 15)
    colnames(result) <- c("mean error", "mean r","zscore","sampledArea","sampledVolume","surface_order","global mean","globalVariance","survey","variomodel","slope start","slope end","length","","")

# import data
    
    ST_dir <- ("D:\\MyPassport\\exported_v3\\ST")
    setwd(ST_dir)
    ST_listfile <- list.files(ST_dir, pattern = '_FM')

    df_raw <- read.csv(ST_listfile[i], header = T) %>% 
      filter( Lat_M != '999', Lon_M != '999')
    df_raw$density <- df_raw$Num_targets/df_raw$Beam_volume_sum*df_raw$Height_mean
    df_raw <- df_raw %>% dplyr::select("Lon_M","Lat_M","density")
    colnames(df_raw) <- c('x','y','density')
    name <- paste(strsplit(ST_listfile[i],split = '_')[[1]][1],strsplit(ST_listfile[i],split = '_')[[1]][2])
    
    # explanatory analysis----
    hist(log(df_raw$density), breaks=20)
    df <- df_raw %>%
      mutate(density = log(df_raw$density)) %>%  # log transform
      filter(density <= 20) 
    colnames(df) <- c('x','y','density')
    
    # check trend and isotropy (omnidirectional) assumptions of stationarity and isotropy----
    df_geoR <- as.geodata(df, coords.col = 1:2, data.col = 3)
    summarydf_geoR <- summary(df_geoR)
    plot(df_geoR) 
    plot(df_geoR, trend='1st') 
    plot(df_geoR, trend='2nd')# check trend and normality-------------trend='1st'

    max.dist = (summarydf_geoR$distances.summary[2] - summarydf_geoR$distances.summary[1])/2
    plot(variog4(df_geoR, trend = '1st', max.dist = max.dist), omnidirectional = T,same.plot=TRUE, legend = F)
    legend(0.006, 1,legend=c("ominid", "0°", "45°", "90°", "135°"),
           col=c("black", "black","black","black","black"), 
           lty = 1:5, cex=0.8) # check anisotropy
    
    # auto model fit----
    coordinates(df) <- ~ x + y 
    crs(df) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    
    trendSurface <- as.formula(density ~ 1); result[i,6] = 0  ## intercept only model
    trendSurface <- as.formula(density ~ x + y);  result[i,6] = 1
    trendSurface <- as.formula(density ~ x + y + I(x^2) + I(y^2) + I(x*y));  result[i,6] = 2
                 
# Variogram fitting
    v0 = variogram(trendSurface,
                   data=df
                   #alpha = c(0, 45, 90, 135)  # Omni direction
                   #cutoff=40, 
                   #width=20,
                   #cloud=T
    )
    plot(v0)
    # assign(paste0('lzn.fit_', i), 
    lzn.fit <- fit.variogram(v0, 
                              model=vgm(c("Exp","Sph","Gau","Lin","Ste")),
                              fit.sills = TRUE, 
                              fit.ranges = TRUE,
                              debug.level = 1,
                              fit.kappa = T)
    plot(v0,lzn.fit)
    
    result[i,10] <- as.character(lzn.fit$model[2])
    assign(paste0('vario_', i),  plot(v0,lzn.fit, main=name))
    assign(paste0('fitted.model_', i),lzn.fit$model[2])

    # make a grid for interpolation and set grid resolution----
    x <- seq(min(df$x)-0.002, max(df$x)+0.002, by = (max(df$x)-min(df$x))/nrow(df)*6) 
    y <- seq(max(df$y)+0.002, min(df$y)-0.002, by = -abs((max(df$y)-min(df$y))/nrow(df)*6))
    grd <- expand.grid(x=x,y=y)
    mapview(df,cex='density')

    plot(df_raw[,c(1,2)], xlim=c(min(grd$x),max(grd$x)), ylim=c(min(grd$y), max(grd$y)), col='red')
    bound <- getpoly()  # manually draw the target regions by click
    
    grd.irregular <- grd[inout(grd,bound),]  # get selected poly
    colnames(grd.irregular) <- c("x", "y")
    coordinates(grd.irregular) <- ~ x + y
    crs(grd.irregular) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    gridded(grd.irregular)  <- TRUE 
    plot(grd.irregular)
 
    assign(paste0('grd.irregular_', i), grd.irregular)

    # kriging
      if(result[i,6] == 0) {
        lzn.kriged <- krige(density ~ 1, df, grd.irregular, lzn.fit)
      } else if (result[i,6] == 1) {
        lzn.kriged <- krige(density ~ x + y, df, grd.irregular, lzn.fit)
      } else {
        lzn.kriged <- krige(density ~ x + y + I(x^2) + I(y^2) + I(x * y), df, grd.irregular, lzn.fit)
      }
    # backtransform data----
    # https://gis.stackexchange.com/questions/237574/backtransformation-of-kriging-predictions-and-variances
    bt <- exp(lzn.kriged$var1.pred + (lzn.kriged$var1.var/2))
    mu_bt <- mean(bt)
    mu_original <- mean(exp(df$density))
    
    # (mu_original-mu_bt)/max.dist
    btt <- bt * (mu_original/mu_bt)  # add a correction factor    
    lzn.kriged$var1.pred <- btt  # rewrite output result
    
    #see kriging result----
    coord <- as.data.frame(lzn.kriged@coords)
    value <- as.data.frame(lzn.kriged$var1.pred)
    krige_result <- cbind(coord,value)
    colnames(krige_result) <- c('x','y','Density')
    ggplot(aes(x=x, y=y),data = krige_result) + 
      geom_tile(aes(fill=Density)) + coord_equal() +
      scale_fill_gradient(low = "yellow", high="red") +
      theme_bw()
    
    result[i,7] <- mean(lzn.kriged$var1.pred)  # artimathic mean 
    result[i,8] <- sum((lzn.kriged$var1.pred - mean(lzn.kriged$var1.pred))^2)/(length(krige_result$Density) - 1)  # global variance/CVgeo
    
    assign(paste0('dfr_density_', i), rasterFromXYZ(krige_result, crs = "+proj=longlat +datum=WGS84 +no_defs"))
    dfr_density <- rasterFromXYZ(krige_result, crs = "+proj=longlat +datum=WGS84 +no_defs")
    
    mapview(dfr_density,layer.name ='Density', col.regions = blues9, na.color=NA)
    
    # cross validation----  
    x <- krige.cv(trendSurface, df, lzn.fit, nmax=40, nfold=5)
     #bubble(x, "residual", main = "5-fold CV residuals")

    # mean error, ideally 0:
    result[i,1] <- mean(x$residual)

    # MSPE, ideally small
    result[i,2] <- mean(x$residual^2)

    # Mean square normalized error, ideally close to 1
    result[i,3] <- mean(x$zscore^2)

######################################################
# bathy for sampled volume and area in the region
######################################################
    bathy_dir <- ("D:\\MyPassport\\exported_v3\\bathy")
    setwd(bathy_dir)
    bathy_listfile <- list.files(bathy_dir, pattern = '_CW')

    #import data
    df_depth <- read.csv(bathy_listfile[i],header = T) %>% 
      filter( Lat_M != '999', Lon_M != '999', Exclude_below_line_depth_mean >= 0) 
    df_depth <- df_depth %>% dplyr::select("Lon_M","Lat_M","Exclude_below_line_depth_mean")
    colnames(df_depth) <- c('x','y','depth')
    coordinates(df_depth) <- ~ x + y 
    crs(df_depth) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    # IDW interpolation depth
    
    IDW <- gstat::idw(depth ~ 1, df_depth, grd.irregular, idp=3.5)
    coord <- as.data.frame(IDW@coords)
    value <- as.numeric(unlist(as.data.frame(IDW$var1.pred)))
    
    IDW_result <- cbind(coord,value)
    colnames(IDW_result) <- c('x','y','depth')
    
    # mapping bathymetry
    assign(paste0('dfr_depth_', i), rasterFromXYZ(IDW_result, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    dfr_depth <- rasterFromXYZ(IDW_result, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
   
    # get area and sampled volume
    area <- raster::area(dfr_depth)
    cellArea <- mean(area@data@values)  # km2
    cellCount <- nrow(IDW_result)  # or = area@ncols*area@nrows - length(dfr[is.na(dfr)])
    result[i,4]<- cellArea * cellCount  # sum of area in km2
    result[i,5] <- sum(IDW_result$depth * cellArea * 10^(-3))  # km3  
    
    assign(paste0('dfr_', i), rasterFromXYZ(krige_result,crs = "+proj=longlat +datum=WGS84 +no_defs"))
    
    result[i,9] <- name
    mapview(dfr_depth)

