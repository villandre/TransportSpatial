# --------------------------------------------------------------------------->
# header
# ---------------------------------------------------------------------------<
library(sp)
library(rgdal)
library(viridis)
library(INLA)
library(igraph)
library(raster)
library(deldir)
library(rgeos)
library(spatstat)
library(osmar) ## This will let us retrieve information from OpenStreetMap

source("~/Projet_Transport/code/R/helperFunctions.R")
source("~/Projet_Transport/code/R/spde-tutorial-functions.R") # Obtained from https://www.math.ntnu.no/inla/r-inla.org/tutorials/spde/R/spde-tutorial-functions.R

#-------------------------------------------------------------------------------------------------- >
# Data import section (always run this section entirely)
# ------------------------------------------------------------------------------------------------- <

quebecCity <- rgdal::readOGR(dsn = "/home/campus/luc.villandre/Projet_Transport/data", layer = "NEW_ROAD_SAFETY_DATA")

crashDataRaw <- read.csv("~/Projet_Transport/data/rawCrashes10years.csv", header = TRUE, stringsAsFactors = FALSE)

colnames(crashDataRaw) <- replace(colnames(crashDataRaw), match("y", colnames(crashDataRaw)), "lat")
colnames(crashDataRaw) <- replace(colnames(crashDataRaw), match("x", colnames(crashDataRaw)), "lon")

# We load both the intersection and the link data with coordinates. For a link to belong to a segment in the grid or the mesh, its midpoint must be within it.

load("~/Projet_Transport/data/intersectionDataWithCoor.Rdata") # The intersection data file I created with OSM data. Creates intersectCovarWithCoor.
intersectionData <- intersectCovarWithCoor
rm(intersectCovarWithCoor)

# The coordinates in intersectionData were inputted as lists, for some reason.

intersectionData$lon <- do.call("c", intersectionData$lon)
intersectionData$lat <- do.call("c", intersectionData$lat)
intersectionData$osm_id <- as.character(intersectionData$osm_id)
rownames(intersectionData) <- intersectionData$osm_id

# I add missing coordinates (approximated by looking at the history in OSM)

missingCoorInfo <- read.table("~/Projet_Transport/data/missingNodeInfo.txt", sep = ",", skip = 2, header = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "numeric", "numeric"))
load("~/Projet_Transport/data/missingNodesCoorQuebecMap.Rdata") # Creates missingNodesCoor

rownames(missingCoorInfo) <- missingCoorInfo$osm_id
missingNodes <- intersectionData$osm_id[is.na(intersectionData$lon)]

intersectionData[missingNodes, c("lon", "lat")] <- missingCoorInfo[missingNodes, c("lon","lat")]

# Note: There are two intersections with the exact same coordinates: "2029" and "2030"

load("~/Projet_Transport/data/linkDataWithCoordinatesChecked.Rdata") # Creates linkDataWithCoorClean.
linkData <- linkDataWithCoorClean
colnames(linkData) <- replace(colnames(linkData), match(c("coor1", "coor2"), colnames(linkData)), c("lon", "lat"))
rm(linkDataWithCoorClean)

maxLon <- -71.129772 # Easier to exclude points east of Quebec. Besides, Josh does not seem to have mapped those accidents to intersections or links.

accidentsForPlot <- subset(crashDataRaw, subset = (lon < maxLon) & ((lon > quebecCity@bbox[1,"min"]) & (lat < quebecCity@bbox[2,"max"]) & (lat > quebecCity@bbox[2,"min"])))

# I extracted coordinates for a polygon that defines the boundaries of our map.

coorFrame <- read.table(file = "~/Projet_Transport/data/boundaryPolygonCoor.info", sep = ",")
colnames(coorFrame) <- c("lat", "lon")

# We exclude all accidents outside the boundaries.

boundaryMatrix <- cbind(lon = coorFrame$lon, lat = coorFrame$lat)

pl.study <- coo2sp(boundaryMatrix)

crashesWithinBoundary <- gIntersection(SpatialPoints(crashDataRaw[ , c("lon", "lat")]), pl.study, byid = TRUE)
rownames(crashesWithinBoundary@coords) <- substr(rownames(crashesWithinBoundary@coords), start = 1, stop = nchar(rownames(crashesWithinBoundary@coords))-2)
# Test
# table(signif(crashDataRaw[rownames(crashesWithinBoundary@coords), "lon"],8) == signif(crashesWithinBoundary@coords[ , "x"],8))
# This returns all TRUE, which is a good sign.
#
# We create a bounding box around a region of interest.

bbHalfWidth <- 5000
bbHalfHeight <- 5000
bbCenterCoord <- c(lon = -71.277977, lat = 46.772500) ## This is a point on Laurier blvd.
load("~/Projet_Transport/data/localOSMdataLaurier.Rdata") # Creates object ua.

# The following commands load the maps for the region centered at Centre Videotron.

# bbHalfWidth <- 4000
# bbHalfHeight <- 4000
# bbCenterCoord <- c(lon = -71.247618, lat = 46.828724) # This is the Videotron Center.
# load("~/Projet_Transport/data/localOSMdataVideotron.Rdata") # Creates object ua.

# bbCenterCoord <- c(lon = -71.396149, lat = 46.844955) # This is a point in Val-Belair

bb <- center_bbox(bbCenterCoord[["lon"]], bbCenterCoord[["lat"]] , bbHalfWidth, bbHalfHeight)

# ------------------------------------------------------------------------------- >
# Preliminary analyses (No need to run this unless you really want to)
# ------------------------------------------------------------------------------- <

# # Do the intersection coordinates make sense? Let's plot them...
#
# plot(x = intersectCovarWithCoor$lon, y = intersectCovarWithCoor$lat, pch = 18, cex = 0.3, col = "blue") # Seems ok to me...
# # We plot the map with the boundary I created.
#
# pdf("~/Projet_Transport/QuebecMapWithBoundary.pdf", paper = "special", width = 40, height = 40)
# plot(quebecCity, xlab = "Longitude", ylab = "Latitude")
# lines(coo2sp(boundaryMatrix), lwd = 3, col = "red")
# dev.off()
#
# # We first plot the raw accident data on the map.
#
# pdf("~/Projet_Transport/mapWithRawAccidents.pdf", paper = "special", width = 40, height = 40)
# plot(quebecCity, xlab = "Longitude", ylab = "Latitude")
# points(x = accidentsForPlot$lon, y = accidentsForPlot$lat, cex = 0.5, pch = 13, col = "red")
# dev.off()
#
#
# # It'd be nice to overlay the map on top of a raster indicating crash frequency!
#
# crashCounts <- bin2(x = cbind(accidentsForPlot$lon, accidentsForPlot$lat), ab = quebecCity@bbox, nbin = c(200,200))
# rasterObject <- raster(ncol = nrow(crashCounts$nc), nrow = ncol(crashCounts$nc), xmx = quebecCity@bbox["x","max"], xmn = quebecCity@bbox["x","min"], ymn = quebecCity@bbox["y","min"], ymx = quebecCity@bbox["y","max"])
# proj4string(rasterObject) <- proj4string(quebecCity)
# values(rasterObject) <- as.vector(crashCounts$nc[ ,ncol(crashCounts$nc):1]) ## Still don't get why this transformation is needed though...
#
# ## Can we overlay the map on the raster?
#
# pdf("~/Projet_Transport/mapWithRawAccidentsRaster.pdf", paper = "special", width = 30, height = 30)
# plot(rasterObject)
# plot(quebecCity, xlim = bbox(rasterObject)["s1", ], ylim = bbox(rasterObject)["s2", ], add = TRUE)
#  ## elements in the vector containing the values for the raster object are drawn row wise, rather than columnwise.
# dev.off()

# ---------------------------------------------------------------------------- >
# We now fit the log-Cox process to our data, using the SPDE approach to approximate the random field.
# Code is based on tutorial in Lund et al. 2016
# ---------------------------------------------------------------------------- <

pointsForMesh <- cbind(intersectionData$lon, intersectionData$lat)
# pointsForMeshWithNA <- rbind(cbind(lon = linkData$lon, lat = linkData$lat),cbind(lon = intersectionData$lon, lat = intersectionData$lat))
# pointsForMesh <- subset(pointsForMeshWithNA, subset = !is.na(pointsForMeshWithNA[ , "lon"]))

mapLimits <- t(bbox(quebecCity))
mapLimits <- rbind(mapLimits["min", ], replace(mapLimits["min", ], 2, mapLimits["max", 2]), mapLimits["max", ], replace(mapLimits["max", ], 2, mapLimits["min", 2]), mapLimits["min", ])

segm.bnd <- inla.mesh.segment(head(boundaryMatrix, n = -1))

#(nv <- (mesh <- inla.mesh.2d(loc.d = mapLimits, off = 0.02, max.e = 0.01))$n)

# OR

(nv <- (mesh <- inla.mesh.2d(loc=pointsForMesh, max.edge=0.1, offset = 0.005, interior = segm.bnd))$n)

# We define the dual mesh.

dmesh <- inla.mesh.dual(mesh)

sum(w <- sapply(1:length(dmesh), function(i) {
    if (gIntersects(dmesh[i,], pl.study))
        return(gArea(gIntersection(dmesh[i,], pl.study)))
    else return(0)

})) ## 48 polygons are totally outside the boundaries. They are therefore given a weight of 0.

# We define Voronoi polygons

# dd <- deldir(mesh$loc[,1], mesh$loc[,2])
# tiles <- tile.list(dd)
#
# w <- sapply(tiles, function(p) {
#   pl <- coo2sp(cbind(p$x, p$y))
#   intersectTest <- tryCatch(gIntersects(pl, pl.study), error = function(e) e)
#   if (!("error" %in% class(intersectTest))) {
#     if (intersectTest) {
#       return(tryCatch(gArea(gIntersection(pl, pl.study)), error = function(e) e))
#     }
#     return(0)
#   }
#   return(intersectTest)
# })

# To understand this section, see Eq. 3-4 in Simpson et al., Biometrika (2016)
y.pp <- rep(0:1, c(nv, nrow(crashesWithinBoundary@coords))) # The observed data. We have no events at the nodes used to construct the mesh, and, trivially, 1 event at each coordinate where an accident was logged within the boundaries.
e.pp <- c(w, rep(0, nrow(crashesWithinBoundary@coords))) # The expected number of events. The expected number of events in the polygons delimited by the dual mesh is given by vector w. The expected number of events at any given point is 0, since the process is continuous over the space.

lmat <- inla.spde.make.A(mesh=mesh, loc=cbind(crashesWithinBoundary@coords[, "x"], crashesWithinBoundary@coords[ , "y"]))
imat <- Diagonal(nv, rep(1, nv))
A.pp <- rBind(imat, lmat)

stk.pp <- inla.stack(data=list(y=y.pp, e=e.pp),
A=list(1,A.pp), tag='pp',
effects=list(list(b0=rep(1,nv+nrow(crashesWithinBoundary@coords))), list(i=1:nv)))

# We define the SPDE model.

sigma0 <- 1
size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
range0 <- size / 5
#kappa0 <- sqrt(8) / range0
kappa0 <- 1
tau0 <- 1 / (sqrt(4 * pi) * kappa0 * sigma0)
spde <- inla.spde2.matern(mesh,
  B.tau = cbind(log(tau0), -1, +1),
  B.kappa = cbind(log(kappa0), 0, -1),
  theta.prior.mean = c(0, 1),
  theta.prior.prec = c(1, 10) # Change to 100
  # prior.kappa = 1)
  )

# For now, I ignore covariates...

pp.res <- inla(y ~ 0 + b0 + f(i, model=spde),
  family='poisson', data=inla.stack.data(stk.pp),
  control.predictor=list(A=inla.stack.A(stk.pp)),
  E=inla.stack.data(stk.pp)$e)

save(pp.res, file = "~/Projet_Transport/inlaFitPointProcess.Rdata", compress = TRUE)

# ----------------------------------------------------------------------------- >
# This section includes code to analyse the input and output of the previously-fitted model.
# -------------------------------------------------------------------------------<

# What does the mesh look like?

# par(mar=c(0,0,0,0))
#
# pdf("~/Projet_Transport/MeshPoint.pdf")
# plot(mesh, asp=1, main='')
# points(pointsForMesh, col=4, pch=19, cex = 0.1); lines(mapLimits, col=3)
# dev.off()

# Plotting the Voronoi polygons

# pdf("~/Projet_Transport/Voronoi.pdf")
# par(mar=c(2,2,1,1), mgp=2:0)
# plot(mesh$loc, asp=1, pch=19, xlab='', ylab='', cex = 0.1)
# for (i in 1:length(tiles)) lines(c(tiles[[i]]$x, tiles[[i]]$x[1]), c(tiles[[i]]$y, tiles[[i]]$y[1]))
# lines(mapLimits, col=3)
# dev.off()

# Plotting the map in the Voronoi polygons.

# pdf("~/Projet_Transport/VoronoiAndMap.pdf", paper = "special", width = 30, height = 30)
# #plot(mesh$loc, asp=1, pch=19, xlab='', ylab='', cex = 0.1)
# plot(quebecCity, xlab = "Longitude", ylab = "Latitude")
# for (i in 1:length(tiles)) lines(c(tiles[[i]]$x, tiles[[i]]$x[1]), c(tiles[[i]]$y, tiles[[i]]$y[1]), col = "grey75", lwd = 0.5)
# lines(mapLimits, col=3)
# points(x = intersectionData$lon, y = intersectionData$lat, col = "green", pch = 19, cex = 0.1)
# points(x = crashesWithinBoundary[, "lon"], y = crashesWithinBoundary[ ,"lat"], col = "red", pch = 17, cex = 0.1)
# dev.off()

# Plotting the map in the dual mesh.

# pdf("~/Projet_Transport/dualMeshAndMap.pdf", paper = "special", width = 30, height = 30)
# plot(quebecCity, xlab = "Longitude", ylab = "Latitude")
# lines(dmesh, col = "grey75")
# points(x = crashesWithinBoundary[, "lon"], y = crashesWithinBoundary[ ,"lat"], col = "red", pch = 17, cex = 0.1)
# dev.off()

# Can we plot the field?

# pgrid0 <- inla.mesh.projector(mesh, dims=c(751,751))
# prd0.m <- inla.mesh.project(pgrid0, pp.res$summary.ran$i$mean)

# pdf("~/Projet_Transport/LatentFieldMean.pdf", paper = "special", width = 30, height = 30)
# lattice::levelplot(prd0.m, col.regions=topo.colors(99), main='latent field mean', xlab='', ylab='', scales=list(draw=FALSE))
# dev.off()

# The model involves a fair bit of smoothing. Maybe it's a problem, maybe it's not. After all, covariates may end up explaining a lot of the intensity.

# We plot the parameters post.

# pp.rf <- inla.spde2.result(pp.res, 'i', spde)
#
# pdf("~/Projet_Transport/spatialParaPost.pdf")
# par(mfrow=c(2,2), mar=c(3,3,1,0.3), mgp=c(2,1,0))
# plot(pp.res$marginals.fix[[1]], type='l', xlab=expression(beta[0]), ylab='Density')
# plot(pp.rf$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma^2), ylab='Density')
# plot(pp.rf$marginals.kappa[[1]], type='l', xlab=expression(kappa), ylab='Density')
# plot(pp.rf$marginals.range.nominal[[1]], type='l', xlab='Nominal range', ylab='Density')
# dev.off()

# ---------------------------------------------------------------------------------->
# The spatstat package lets us fit a non-homogeneous point process on a linear network.
# We'll use the lppm function to do so.
# We identify an arbitrary region of interest and define a bounding box around it.
# Then, we extract the data from OSM and build the required objects.
# ----------------------------------------------------------------------------------
# Phase 1: Preparing data for lppm (importing map data from OSM)
# ----------------------------------------------------------------------------------<

# We identify crashes in the bounding box.

bboxCoor <- rbind(c(bb[["left"]], bb[["bottom"]]), c(bb[["right"]], bb[["bottom"]]), c(bb[["right"]], bb[["top"]]), c(bb[["left"]], bb[["top"]]))
boundBoxPoly <- coo2sp(bboxCoor)
localCrashes <- gIntersection(crashesWithinBoundary, boundBoxPoly, byid = TRUE)
rownames(localCrashes@coords) <- substr(rownames(localCrashes@coords), start = 1, stop = nchar(rownames(localCrashes@coords)) - 2)

# We create an empty raster map from the covariate information

squareCellEdgeLengthInMeters <- 5
numRasterCols <- ceiling(2*bbHalfWidth/squareCellEdgeLengthInMeters)
numRasterRows <- ceiling(2*bbHalfHeight/squareCellEdgeLengthInMeters)

quebecRaster <- raster(ncol = numRasterCols, nrow = numRasterRows, xmx = bb["right"], xmn = bb["left"], ymn = bb["bottom"], ymx = bb["top"], crs = proj4string(quebecCity))

intersectionAndLinkCovInfo <- rbind(subset(intersectionData, select = -c(minor, major, fatal, crash, model_id)), subset(linkData, select = -c(minor, major, fatal, crash, model_id, recorded, name)))
intersectionAndLinkCovInfo$logTrips <- log(intersectionAndLinkCovInfo$trips)
intersectionAndLinkCovInfo$highway <- factor(x = intersectionAndLinkCovInfo$highway, levels = c("primary", "motorway", "residential"), ordered = TRUE)
covariateNames <- c("highway", "logTrips", "drate", "ci", "v", "cvs")
localIntersAndLinks <- gIntersection(SpatialPoints(intersectionAndLinkCovInfo[ , c("lon", "lat")]), boundBoxPoly, byid = TRUE)
rownames(localIntersAndLinks@coords) <- substr(rownames(localIntersAndLinks@coords), start = 1, stop = nchar(rownames(localIntersAndLinks@coords)) - 2)
localCoordsMat <- localIntersAndLinks@coords
colnames(localCoordsMat) <- c("lon", "lat")

# We now retrieve OSM data directly (loaded in the data input section).

# ua <- get_osm(bb)
# save(ua, file = "~/Projet_Transport/data/localOSMdataValBelair.Rdata", compress = TRUE) # Better to save it in memory, to spare us the trouble of having to query OSM every time.
# We keep only features of interests, namely public roads...

codesForRoadsInMap <- find(ua, way(tags((k == "highway") & (v %in% c("primary", "residential", "motorway", "secondary", "trunk", "tertiary", "motorway_link", "trunk", "primary_link", "secondary_link", "tertiary_link")))))
nodesFromOSM <- as.matrix(subset(ua$nodes$attrs, select = c("lon", "lat")))
rownames(nodesFromOSM) <- ua$nodes$attrs$id

linesAndNodesInMap <- find_down(ua, way(codesForRoadsInMap))
mapFeaturesOfInterest <- subset(ua, ids = linesAndNodesInMap)

# SpatialLinesDataFrame object
mapFeaturesOfInterestAsSp <- as_sp(mapFeaturesOfInterest, "lines")

compLocalMapAsGraph <- as_igraph(mapFeaturesOfInterest)

nodesCoordinates <- nodesFromOSM[names(V(compLocalMapAsGraph)),] # No NA's!

# We create a raster for intersections...
# Intersections are nodes shared by two or more ways. I wonder if finding intersection between node sets in all identified ways could help identify the different intersections?

radiusIntersectionsInKm <- 0.01
nodesInEachWay <- lapply(codesForRoadsInMap, function(wayCode) {
  find_down(ua, way(wayCode))$node_ids
})
names(nodesInEachWay) <- codesForRoadsInMap

allNodes <- do.call("c", nodesInEachWay)
nodeFrequency <- table(allNodes)
intersectionNodes <- names(nodeFrequency)[nodeFrequency > 1] # This will probably catch arbitrary subdivisions in a road...
intersectionCoords <- nodesCoordinates[intersectionNodes,]

# We also have a bit of information about intersections in the crash data...

localCrashesWithInfo <- crashDataRaw[rownames(localCrashes@coords), ]
localIntersectionsFromRaw <- subset(localCrashesWithInfo, match_type == "inter")
localIntersectionsFromRawCoords <- localIntersectionsFromRaw[ , c("lon", "lat")]

intersectionRaster <- populateRaster(rasterObject = quebecRaster, coordsMat = rbind(intersectionCoords, localIntersectionsFromRawCoords), valuesVec = rep(1, nrow(intersectionCoords) + nrow(localIntersectionsFromRawCoords)), radiusInKm = radiusIntersectionsInKm)
values(intersectionRaster) <- replace(values(intersectionRaster), which(is.na(values(intersectionRaster))), 0)

# ----------------------------------------------------------------------------- >
# Phase 2: Preparing linnet object and fitting the model.
# ----------------------------------------------------------------------------- <
# We need to produce a linnet object from the osm data we extracted.
# To produce this object, we need to have the location of all vertices and an adjacency matrix / edge matrix
# This can be obtained by producing an igraph object and extracting the necessary components.

# Some nodes fall outside the bounding box. We need to remove them.

netPointsInsideBbox <- gIntersection(SpatialPoints(nodesCoordinates), boundBoxPoly, byid = TRUE)

## Add a mark indicating whether each point is for a motorway, primary, or residential link. ##

nodesCoordinatesPPP <- ppp(as.vector(netPointsInsideBbox@coords[, "x"]), as.vector(netPointsInsideBbox@coords[, "y"]), window = owin(x = c(bb[["left"]], bb[["right"]]), y = c(bb[["bottom"]], bb[["top"]]))) ## vertices in linnet.
nodesInsideBbox <- substr(rownames(netPointsInsideBbox@coords), start = 1, stop = nchar(rownames(netPointsInsideBbox@coords)) - 2)

netAdjacencyMatrix <- as_adjacency_matrix(as.undirected(compLocalMapAsGraph, mode = "collapse"), type = "both")
subNetAdjacencyMatrix <- netAdjacencyMatrix[nodesInsideBbox, nodesInsideBbox]

localLinNet <- linnet(vertices = nodesCoordinatesPPP, m = as.matrix(subNetAdjacencyMatrix != 0), sparse = TRUE)
componentIndices <- connected(localLinNet)
componentFreqTable <- table(componentIndices)
thinnedNet <- thinNetwork(localLinNet, retainvertices = which(componentIndices == as.numeric(names(componentFreqTable)[which.max(componentFreqTable)]))) # It is connected!

# TEST ##########################33
# pdf("~/Projet_Transport/localLinNet.pdf", paper = "special", width = 30, height = 30)
# plot(thinnedNet)
# dev.off() # It seems to work: disconnected components seem to have disappeared.

localLPP <- lpp(localCrashes@coords, L = thinnedNet)

lapply(setdiff(covariateNames, "highway") , FUN = function(covariateName) {
  myFun <- function(x, y, seg, tp) { # The formals are imposed by the linfun function.
    coordMatrixForExtract <- as.matrix(cbind(x,y))
    invisible(covariateImputeKriging(pointCoordinatesMatrix = localCoordsMat, covValues = intersectionAndLinkCovInfo[rownames(localCoordsMat), covariateName, drop = TRUE], pointsToImpute = coordMatrixForExtract, maxdist = Inf, nmin = 0, nmax = Inf, useIDW = TRUE, idp = 2, returnCovar = FALSE))
  }
  funName <- paste(covariateName, "Fun", sep = "")
  linfunToAssign <- linfun(myFun, L = thinnedNet)
  assign(funName, value = linfunToAssign, envir = .GlobalEnv)
  NULL
}) # Function op
# Test: vFun(x = localCrashes@coords[1:1000 , "x"], y = localCrashes@coords[1:1000 , "y"])

intersectionFun <- function(x, y, seg, tp) { # The formals are imposed by the linfun function.
  coordMatrixForExtract <- as.matrix(cbind(x,y))
  invisible(extract(intersectionRaster, coordMatrixForExtract))
}

# The function associating street type to points is slightly more complicated. Each segment in the vector map has a type, e.g. motorway, primary, trunk, which I extract. I then re-classify those values, grouping certain values into a unique factor level, e.g. primary_link, secondary_link, "tertiary_link", "motorway_link" goes into the "ramp" level.
# I then assign a priority in assignment, e.g. a point on an intersection of a primary and a residential road will be classified as falling on a primary road.

rightTable <- subset(mapFeaturesOfInterest$ways$tags, k == "highway", select = c(id,v))
mergedTable <- merge(mapFeaturesOfInterestAsSp@data, rightTable, by = "id", all.x = TRUE, all.y = FALSE, sort = FALSE) # Same number of rows as the table on the left, good sign! The ordering of id hasn't changed either.
mergedTable$v <- factor(as.character(mergedTable$v), ordered = TRUE)
levels(mergedTable$v) <- list(ramp = c("primary_link", "secondary_link", "tertiary_link", "motorway_link"), primary = "primary", secondary = "secondary", motorway = c("trunk", "motorway"), residential = "residential", tertiary = "tertiary") # We reformat the street type column. "trunk" is a major highway, and is recoded "motorway". "*_link" are exit ramps or access points to a highway. For now, we group them in a  "ramp" category. Note that the factor levels are ordered: ramp < primary < secondary < motorway < residential < tertiary. The ordering may reflect the priority of assignment at intersections.

rownames(mergedTable) <- as.character(mergedTable$id) # I made sure the id's were unique.

linesInMapReformat <- mapply(id = rownames(mergedTable), LineObject = mapFeaturesOfInterestAsSp@lines, function(id, LineObject) {
  as.data.frame(LineObject@Lines[[1]]@coords)
}, SIMPLIFY = FALSE)
names(linesInMapReformat) <- rownames(mergedTable)

minLons <- sapply(linesInMapReformat, function(x) min(x$lon))
maxLons <- sapply(linesInMapReformat, function(x) max(x$lon))
minLats <- sapply(linesInMapReformat, function(x) min(x$lat))
maxLats <- sapply(linesInMapReformat, function(x) max(x$lat))

lapply(levels(mergedTable$v) , FUN = function(covariateName) {
  myFun <- function(x, y, seg, tp) { # The formals are imposed by the linfun function.
    coordMatrixForExtract <- as.matrix(cbind(x,y))
    colnames(coordMatrixForExtract) <- c("lon", "lat")
    funForApply <- function(lonLat) {
      compatibleIDs <- coordToLineID(minLons = minLons, maxLons = maxLons, minLats = minLats, maxLats = maxLats, lonCoord = lonLat["lon"], latCoord = lonLat["lat"])
      if (length(compatibleIDs) > 1) {
        distsToLine <- sapply(compatibleIDs, FUN = function(id) {
          numRows <- nrow(linesInMapReformat[[id]])
          distancesForRoadSegment <- sapply(2:numRows, FUN = function(rowIndex) {
            distPointToLine(point = lonLat, point1onLine = linesInMapReformat[[id]][rowIndex, ], point2onLine = linesInMapReformat[[id]][rowIndex - 1, ])
          })
          min(distancesForRoadSegment)
        })
        # If a point is on an intersection, it should return distance of 0 for two line segments, which is why I don't use which.min, which only returns one value.
        compatibleIDs <- compatibleIDs[which(distsToLine == min(distsToLine))]
      }
      if (is.na(compatibleIDs)[[1]]) return(NA)

      as.numeric(as.character(min(mergedTable[compatibleIDs,]$v)) == covariateName)
    }
    apply(coordMatrixForExtract, MARGIN = 1, FUN = funForApply)
  }
  funName <- paste(covariateName, "Fun", sep = "")
  linfunToAssign <- linfun(myFun, L = thinnedNet)
  assign(funName, value = linfunToAssign, envir = .GlobalEnv)
  NULL
})

# Test: ################
# crashesToConsider <- 1:nrow(localCrashes@coords)
# primaryIndices <- primaryFun(x = localCrashes@coords[crashesToConsider, "x"], y = localCrashes@coords[crashesToConsider, "y"])
# secondaryIndices <- secondaryFun(x = localCrashes@coords[crashesToConsider, "x"], y = localCrashes@coords[crashesToConsider, "y"])
# residentialIndices <- residentialFun(x = localCrashes@coords[crashesToConsider, "x"], y = localCrashes@coords[crashesToConsider, "y"])
# motorwayIndices <- motorwayFun(x = localCrashes@coords[crashesToConsider, "x"], y = localCrashes@coords[crashesToConsider, "y"])
# rampIndices <- rampFun(x = localCrashes@coords[crashesToConsider, "x"], y = localCrashes@coords[crashesToConsider, "y"])
# tertiaryIndices <- tertiaryFun(x = localCrashes@coords[crashesToConsider, "x"], y = localCrashes@coords[crashesToConsider, "y"])
#
# pdf("~/Projet_Transport/streetTypeFlagTestLaurierSector.pdf", paper = "special", width = 30, height = 30)
# plot(localLPP, cols = "black")
# plot(localCrashes[which(residentialIndices == 1)], col = "green", add = TRUE, pch = 8, cex =  7)
# plot(localCrashes[which(primaryIndices == 1)], col = "purple", add = TRUE, pch = 9, cex =  7)
# plot(localCrashes[which(secondaryIndices == 1)], col = "red", add = TRUE, pch = 10, cex =  7)
# plot(localCrashes[which(motorwayIndices == 1)], col = "blue", add = TRUE, pch = 7, cex =  7)
# plot(localCrashes[which(rampIndices == 1)], col = "orange", add = TRUE, pch = 6, cex =  7)
# plot(localCrashes[which(tertiaryIndices == 1)], col = "gold", add = TRUE, pch = 5, cex =  7)
# dev.off() # Does seem ok. Map makes a lot of sense. I am surprised to see certain crashes close to secondary roads assigned to residential roads.
# sum(primaryIndices)+sum(residentialIndices)+sum(motorwayIndices)+sum(secondaryIndices)+sum(rampIndices)+sum(tertiaryIndices) ## Sum is 2059, which matches the number of points in localCrashes. We're ok.
######################

fittedLocalLPPM <- lppm(localLPP ~ intersectionFun + logTripsFun + rampFun + primaryFun + secondaryFun + tertiaryFun + residentialFun + vFun + drateFun + ciFun + cvsFun + vFun:(primaryFun + secondaryFun + tertiaryFun + residentialFun) + drateFun:(primaryFun + secondaryFun + tertiaryFun + residentialFun) + ciFun:(primaryFun + secondaryFun + tertiaryFun + residentialFun) + cvsFun:(primaryFun + secondaryFun + tertiaryFun + residentialFun) + logTripsFun:(primaryFun + secondaryFun + tertiaryFun + residentialFun)) # Base level is motorway.

# getLatexTableFromLPPM(fittedLocalLPPM, digits = 3)
# refPoint <- c(-73.645086, 45.544921)
# fittedValueFromCoefVec(coefVec = coef(fittedLocalLPPM), newValues = c(vFun = 35*1000/3600, logTripsFun = log(7000), ciFun = mean(intersectionAndLinkCovInfo$ci), cvsFun = mean(intersectionAndLinkCovInfo$cvs), drateFun = mean(intersectionAndLinkCovInfo$drate), primaryFun = 0, secondaryFun = 1, rampFun = 0, tertiaryFun = 0, residentialFun = 0, intersectionFun = 0), segmentLength = 0.009)
# spDistsN1(pts = matrix(refPoint, 1, 2), pt = refPoint+c(0,0.009), longlat = TRUE) equals approx. 1km.
# What would be the intensity on an arbitrary primary segment
# ----------------------------------------------------------------------------------------------- >
# A few plots to analyse the LPPM output and inputs...

pdf("~/Projet_Transport/LaurierSector.pdf", paper = "special", width = 30, height = 30)
plot(localLPP, cols = "red")
dev.off()

# Plotting intensity (Laurier sector)

pdf("~/Projet_Transport/fittedIntensityLaurier.pdf", paper = "special", width = 30, height = 30)
plot(fittedLocalLPPM, type = "cif", main = NULL, ribwid = 0.025, useRaster = TRUE)
plot(localLPP, add = TRUE, cols = "green", col = "cadetblue1")
dev.off()

# Plotting values obtained through kriging...
# mapPointsInRaster <- rasterize(x = mapFeaturesOfInterestAsSp, y = quebecRaster, field = 1, update = TRUE, updateValue = "NA")
# intersectionPoints <- which(!is.na(values(mapPointsInRaster)))
pointsForImputation <- rasterToPoints(quebecRaster, spatial = TRUE)
imputedValues <- vFun(pointsForImputation@coords[ , "x"], pointsForImputation@coords[ , "y"])
newRaster <- quebecRaster
values(newRaster) <- imputedValues

pdf("~/Projet_Transport/speedIDWrasterLaurier.pdf", paper = "special", width = 30, height = 30)
plot(localLPP, main = NULL)
plot(newRaster, add = TRUE, alpha = 0.7)
dev.off()
