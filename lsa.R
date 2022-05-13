for (jjj in rev(4*c(30,35,40,45,50,55,60,65,70,75,80,85,90,95))) {
  time1 <- timestamp()
  setwd("D:/YTC/")
  folder.name <- "200000_last"
  dir.create(file.path (folder.name), showWarnings = F)
  wind.size = jjj; mov.step = wind.size / 4; DEM.path = paste0(file.path (folder.name),"/DEM.tif"); border.path = paste0(file.path (folder.name), "/border.shp")
  
  #####
  ### LOAD LIBRARIES
  library(raster)
  library(parallel)
  library(rgeos)
  library(rgdal)
  
  #####
  ### READ RASTER DEM
  gt <- raster(DEM.path, band = 1) # считываем исходный растр
  par(mfrow = c(1,2))
  plot(gt)
  area <- readOGR(border.path)
  area <- spTransform(area, crs(gt))
  area <- gBuffer(spgeom = area, width = wind.size*res(gt)[1])
  plot(area, add = T)
  
  i.tiles <- ceiling((gt@ncols - floor(wind.size/2)) / mov.step) # столбцов в нарезке
  j.tiles <- ceiling((gt@nrows - floor(wind.size/2)) / mov.step) # строк в нарезке
  ld_corner <- xyFromCell(gt, cell = ncell(gt) - ncol(gt) + 1)
  ld_center <- ld_corner + res(gt)*wind.size/2
  sg1 <- SpatialGrid(GridTopology(cellcentre.offset = as.vector(ld_center), cellsize = mov.step*res(gt), cells.dim = c(i.tiles,j.tiles)))
  sp1 <- SpatialPoints(sg1)
  
  sp1 <- gIntersection(sp1, area)
  points(sp1[sample(1:length(sp1), 1000),])
  
  cl <- makeCluster(getOption("cl.cores", 62))
  clusterCall(cl = cl, fun = function() library(raster))
  clusterExport(cl = cl, varlist = c("sp1", "gt", "wind.size"))
  system.time(tt <- parallel::parSapply(cl = cl, X = 1:length(sp1), FUN = function(x) extent(sp1[x,]) + res(gt)*wind.size))
  clusterExport(cl = cl, varlist = c("tt"))
  
  dir.mat <- matrix(data = NA, nrow = wind.size, ncol = wind.size)
  for (i in 1:wind.size) {
    for (j in 1:wind.size) {
      if ((j - ((wind.size/2) + 1)) < 0) {
        dir.mat[i, j] <- atan((i - ((wind.size/2) + 1)) / (j - ((wind.size/2) + 1))) + pi
      } else {
        dir.mat[i, j] <- atan((i - ((wind.size/2) + 1)) / (j - ((wind.size/2) + 1)))
      }
    }
  }
  dir.mat <- dir.mat + pi/2
  part <- ceiling(sqrt(wind.size))
  
  foo1 <- function(gt = gt, tt = tt, i) {
    (ras_crop <- crop(gt, tt[[i]]))
    # plot(ras_crop)
    trend1 <- geodiv::fitplane(ras_crop, 1)
    names(ras_crop) <- NULL
    ras_crop <- as.matrix(ras_crop) - trend1
    m1 <- matrix(0, nrow = wind.size, ncol = wind.size)
    m1[(wind.size - nrow(ras_crop) + 1):wind.size,1:ncol(ras_crop)] <- ras_crop
    ras_crop <- m1
    four.im <- fft(ras_crop)
    abs_fi  <- abs(four.im) / prod(dim(ras_crop))
    
    val <- unique(round(sort(abs_fi, decreasing = T), di = 6))
    
    indc <- (which(abs_fi > val[part+1]))
    four.im2 <- four.im
    four.im2[-indc] <- 0
    
    plt <- Re(fft(four.im2, inverse = T) / prod(dim(ras_crop)))
    
    (max.imp <- 1 - (sd(ras_crop - plt) / sd(ras_crop)))
    (max.mag <- 4*max(abs_fi))
    
    dir <- dir.mat[indc]
    weit <- abs_fi[indc]
    xyt <- cbind(weit*sin(dir), weit*cos(dir))
    t <- prcomp(xyt)$rotation[1,1] / prcomp(xyt)$rotation[2,1]
    (pr.dir1 <- sin(pi*(90 - 180*atan(t)/pi))/180)
    (pr.dir2 <- cos(pi*(90 - 180*atan(t)/pi))/180)
    (pr.deg <- cor(xyt[,2] , xyt[,1])^2)
    
    if (!is.na(max(abs_fi))) {
      (max.freq.lat <- (which.max(abs_fi) %% wind.size) - 1)
      (max.freq.long <- floor(which.max(abs_fi) / wind.size))
      (wavelen <- (res(gt)[1]*wind.size) / (max.freq.lat^2 + max.freq.long^2)^0.5)
    }
    
    fr <- res(gt)[1]*wind.size / (ifelse(col(ras_crop) - 1 < ncol(ras_crop) + 1 - col(ras_crop), col(ras_crop) - 1, ncol(ras_crop) + 1 - col(ras_crop))^2 + 
                                    ifelse(row(ras_crop) - 1 < nrow(ras_crop) + 1 - row(ras_crop), row(ras_crop) - 1, nrow(ras_crop) + 1 - row(ras_crop))^2)^0.5
    fr <- round(fr, 3)
    dt <- aggregate(as.vector(abs_fi), by = list(as.vector(fr)), mean)
    nlsm <- minpack.lm::nlsLM(formula = x ~ a*exp(1)^(-l*Group.1), data = as.data.frame(dt), start = list(a = cellStats(gt, max), l = 1))
    (a_exp <- coefficients(nlsm)[1])
    (l_exp <- coefficients(nlsm)[2])
    (dev_exp <- nlsm$m$deviance())
    c(max.imp, max.mag, pr.dir1, pr.dir2, pr.deg, max.freq.lat, max.freq.long, wavelen,a_exp, l_exp, dev_exp )
  }
  
  (clusterCall(cl = cl, fun = function() {library(minpack.lm) ; library(geodiv)}))
  clusterExport(cl = cl, varlist = c("foo1", "part", "dir.mat"))
  timestamp()
  system.time(gt.stack <- parallel::parLapply(cl = cl, X = seq_along(tt), fun = function(x) try(foo1(gt = gt, tt = tt, i = x))))
  
  index1 <- unlist(lapply(X = gt.stack, FUN = function(x) ifelse(length(grep(x, pattern = "Error")) != 0, F, T)))
  
  
  test_m <- matrix(NA, nrow = length(tt), ncol = 11)
  test_m[index1,] <- do.call(rbind, gt.stack[index1])
  
  plot(rasterFromXYZ(cbind(coordinates(sp1), test_m[,2])))
  time2 <- timestamp()
  time2 ; time1
  
  vars <- c("max.imp", "max.mag", "pr.dir1", "pr.dir2", "pr.deg", "max.freq.lat", "max.freq.long", "wavelen", "a_exp", "l_exp", "dev_exp" )
  for ( i in 1:length(vars)) {
    writeRaster(rasterFromXYZ(cbind(coordinates(sp1), test_m[,i])), paste0(folder.name,"/",wind.size,"_",vars[i],".tif"), "GTiff", overwrite=TRUE)
  }
  stopCluster(cl)
}