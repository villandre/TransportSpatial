library(sp)
library(raster)

identifyClosestPoints <- function(coorVector, coorMatrix, numNeighbours)
{
  distFromPoint <- apply(coorMatrix, MARGIN = 1, FUN = function(x)
  {
    sum((coorVector - x)^2)
  })
  which(rank(distFromPoint) < numNeighbours+2)
}

makeCircularPolygon <- function(numPointsCircle, centerCoor, radius)
{
  posDomain <- seq(from = 0, to = 2*pi, length.out = numPointsCircle)
  xCoor <- radius*cos(posDomain) + centerCoor[[1]]
  yCoor <- radius*sin(posDomain) + centerCoor[[2]]
  simplePolygon <- Polygon(cbind(xCoor, yCoor), hole = FALSE)
  simplePolygon
}

plotCircle <- function(LonDec, LatDec, Km, plotIt = TRUE) {#Corrected function
    #LatDec = latitude in decimal degrees of the center of the circle
    #LonDec = longitude in decimal degrees
    #Km = radius of the circle in kilometers
    ER <- 6371 #Mean Earth radius in kilometers. Change this to 3959 and you will have your function working in miles.
    AngDeg <- seq(from = 3, to = 360, by = 3) #angles in degrees
    Lat1Rad <- LatDec*(pi/180)#Latitude of the center of the circle in radians
    Lon1Rad <- LonDec*(pi/180)#Longitude of the center of the circle in radians
    AngRad <- AngDeg*(pi/180)#angles in radians
    Lat2Rad <-asin(sin(Lat1Rad)*cos(Km/ER)+cos(Lat1Rad)*sin(Km/ER)*cos(AngRad)) #Latitude of each point of the circle rearding to angle in radians
    Lon2Rad <- Lon1Rad+atan2(sin(AngRad)*sin(Km/ER)*cos(Lat1Rad),cos(Km/ER)-sin(Lat1Rad)*sin(Lat2Rad))#Longitude of each point of the circle rearding to angle in radians
    Lat2Deg <- Lat2Rad*(180/pi)#Latitude of each point of the circle rearding to angle in degrees (conversion of radians to degrees deg = rad*(180/pi) )
    Lon2Deg <- Lon2Rad*(180/pi)#Longitude of each point of the circle rearding to angle in degrees (conversion of radians to degrees deg = rad*(180/pi) )
    if (plotIt)
    {
      polygon(Lon2Deg,Lat2Deg,lty=2)
    }
    coordTable <- cbind(Lon2Deg,Lat2Deg)
    coordTable <- rbind(coordTable, head(coordTable, 1)) ## To close the circle...
    Polygon(coords = cbind(Lon2Deg,Lat2Deg), hole = FALSE)
}

## From Hadley Wickham's Advanced R. The result of this infix function is returning a default value when a function evaluates to NULL.
`%||%` <- function(a, b) if (!is.null(a)) a else b
## Replacement functions can also be handy, e.g. colnames(aaa) <- c("foo", "bar") is an example of this. Note however that a copy of colnames() is created in the process.

# This function transforms a matrix of coordinates for a polygon into a SpatialPolygon objet. Very useful to check if points are contained within an arbitrary zone of the map.

coo2sp <- function(coo) {
  n <- nrow(coo)
  if (any(coo[1,]!=coo[n,]))
  coo <- coo[c(1:n,1),]
  SpatialPolygons(list(Polygons(list(Polygon(coo)), '0')))
}

getPointsInPolygon <- function(pointCoordinatesMatrix, polygonMatrix) {
  pointsInSpatial <- SpatialPoints(pointCoordinatesMatrix)
  polygonInSpatial <- coo2sp(polygonMatrix)
  pointInIntersection <- intersect(pointsInSpatial, polygonInSpatial)
  intersectAsMatrix <- pointInIntersection@coords
  colnames(intersectAsMatrix) <- colnames(pointCoordinatesMatrix)
  intersectAsMatrix
}

populateRaster <- function(rasterObject, coordsMat, valuesVec, radiusInKm, categFun = min, continuousFun = mean) {
  circlesPolygonsList <- lapply(1:nrow(coordsMat), function(lineNum) {
    aPolygon <- plotCircle(LonDec = coordsMat[lineNum, "lon"], LatDec = coordsMat[lineNum, "lat"], Km = radiusInKm, plotIt = FALSE)
    Polygons(srl = list(aPolygon), ID = as.character(lineNum))
  })
  circlesPolygons <- SpatialPolygons(Srl = circlesPolygonsList)
  funToUse <- continuousFun
  if ("factor" %in% class(valuesVec)) {
    funToUse <- categFun
  }
  newRaster <- rasterize(x = circlesPolygons, y = rasterObject, field = as.numeric(valuesVec), fun = funToUse, update = TRUE, updateValue = "NA")
  newRaster
}

covariateImputeKriging <- function(pointCoordinatesMatrix, covValues, pointsToImpute, maxdist = Inf, nmin = 0, nmax = Inf, useIDW = FALSE, idp = 2, returnCovar = FALSE) {
  require(gstat)
  colnames(pointCoordinatesMatrix) <- c("x","y")
  covFrameWithCoor <- as.data.frame(cbind(covariate = covValues, pointCoordinatesMatrix))
  coordinates(covFrameWithCoor) <- ~x+y

  covFrameWithCoor@data$covariate <- covValues - mean(covValues) # Do we really need to centralize values?

  pointWithProj <- SpatialPoints(pointsToImpute)
  if (useIDW) {
    idwEst <- gstat::idw(covariate~1, locations = covFrameWithCoor, newdata = pointWithProj, maxdist = maxdist, nmin = nmin, nmax = nmax, idp = idp)
    return(idwEst@data$var1.pred + mean(covValues))
  }
  lzn.vgm <- variogram(covariate~1, data = covFrameWithCoor) # calculates sample variogram values
  nuggetVal <- 0.95*min(lzn.vgm$gamma)
  psillVal <- tail(lzn.vgm$gamma, n = 1) - nuggetVal
  lzn.fit <- fit.variogram(lzn.vgm, model=vgm(psill = psillVal, model = "Exp", nugget = nuggetVal), fit.sills = FALSE, fit.ranges = FALSE) # fit model
  if (returnCovar) {
    lzn.kriged.withCovar <- krige0(formula = covariate~1, data = covFrameWithCoor, newdata = pointWithProj, model=lzn.fit, beta = 0, computeVar = TRUE, fullCovariance = TRUE)
    return(lzn.kriged.withCovar)
  }
  require(geoR)
  ordiKrig <- krige.conv(geodata = list(coords = covFrameWithCoor@coords, data = covFrameWithCoor@data),locations = pointWithProj@coords, krige = krige.control(cov.pars = c(psillVal, median(lzn.vgm$dist))))
  # lzn.kriged <- krige(covariate~1, locations = covFrameWithCoor, newdata = SpatialPoints(pointsToImpute), model=lzn.fit, maxdist = maxdist, nmin = nmin, nmax = nmax, beta = 0)
  ordiKrig@data$var1.pred + mean(covValues)
}

coordToLineID <- function(minLons, maxLons, minLats, maxLats, lonCoord, latCoord) { # minLons must have named elements.
  # minLons <- sapply(linesList, function(x) min(x$lon))
  # maxLons <- sapply(linesList, function(x) max(x$lon))
  # minLats <- sapply(linesList, function(x) min(x$lat))
  # maxLats <- sapply(linesList, function(x) max(x$lat))

  # We include only segments that are in range.
  compatibleLines <- which((minLons <= lonCoord) & (maxLons >= lonCoord) & (minLats <= latCoord) & (maxLats >= latCoord))
  if (length(compatibleLines) == 0) return(NA)

  names(minLons)[compatibleLines]
}

distPointToLine <- function(point, point1onLine, point2onLine) {
  abs((point2onLine[[2]] - point1onLine[[2]])*point[[1]] - (point2onLine[[1]] - point1onLine[[1]])*point[[2]] + point2onLine[[1]]*point1onLine[[2]] - point2onLine[[2]]*point1onLine[[1]])/sqrt((point2onLine[[2]] - point1onLine[[2]])^2 + (point2onLine[[1]] - point1onLine[[1]])^2)
}

getLatexTableFromLPPM <- function(lppmOutput, digits = 3) {
  require(xtable)
  filename <- "/tmp/lppmOut.info"
  previousWidth <- options()$width
  options(width = 200) # We chage the print width to prevent the print output to be split on multiple lines, which screws up the formatting in the table.
  capture.output(lppmOutput, file = filename, append = FALSE)
  options(width = previousWidth) # It must be reset afterwards, since changing system variables produces a side-effect.
  argsForCommandBeg <- paste("-n Estimate", filename)
  outputForBeginning <- system2(command = "grep", args = argsForCommandBeg, stdout = TRUE, stderr = TRUE)
  splitNumbers <- as.numeric(strsplit(outputForBeginning, split = ":")[[1]][1])
  numLinesToSkip <- splitNumbers
  argsForCommandEnd <- paste("-n \"Original data\"", filename)
  outputForEnd <- system2(command = "grep", args = argsForCommandEnd, stdout = TRUE, stderr = TRUE)
  endLineNum <- as.numeric(strsplit(outputForEnd, split = ":")[[1]][1])
  tableAsText <- readLines(filename)[(numLinesToSkip+1):(endLineNum-1)]
  system2(command = "rm", args = filename)
  tableAsTextSplitAtStar <- strsplit(tableAsText, split = "\\*")
  linesWithoutStars <- lapply(tableAsTextSplitAtStar, FUN = function(x) {
    lineWithoutStars <- x
    if (length(x) > 1) {
      lineWithoutStars <- x[[1]]
    }
    stringsWithEmpty <- strsplit(lineWithoutStars, split = " ")[[1]]
    nonEmptyStrings <- stringsWithEmpty[stringsWithEmpty != ""]
    if (length(x) == 1) {
      nonEmptyStrings <- head(nonEmptyStrings, n = -1)
    }
    nonEmptyStrings
  })
  tableRemade <- as.data.frame(do.call("rbind", linesWithoutStars), stringsAsFactors = FALSE)
  theRownames <- tableRemade[ , 1]
  tableRemade <- tableRemade[ , -1]
  tableRemade <- sapply(tableRemade, function(x) as.numeric(x))
  colnames(tableRemade) <- c("Estimate", "SE", "CI-low", "CI-high")
  rownames(tableRemade) <- theRownames
  xtable(tableRemade, digits = digits)
}

# The function will return the fitted intensity for a given set of covariates if length is 0. If it is greater than 0, then it will return the expected number of crashes, i.e. exp(fitted intensity * segment length)

fittedValueFromCoefVec <- function(coefVec, newValues, segmentLength = 0) {
  fittedValues <- sapply(setdiff(names(coefVec), "(Intercept)"), FUN = function(varName) {
    if (grepl(pattern = ":", x = varName)) {
      splitName <- strsplit(varName, ":")[[1]]
      return(coefVec[[varName]]*newValues[[splitName[[1]]]]*newValues[[splitName[[2]]]])
    }
    coefVec[[varName]]*newValues[[varName]]
  })
  fittedIntensity <- sum(fittedValues) + coefVec[["(Intercept)"]]
  if (segmentLength > 0) {
    return(exp(fittedIntensity*segmentLength))
  }
  fittedIntensity
}
