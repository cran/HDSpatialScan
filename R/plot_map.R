########################################################################################################################################################
##' @title Schema of the clusters
##'
##' @description This function plots a schema of the sites and the clusters
##'
##' @param output_clusters list. List of the sites in the clusters: it is the sites_clusters of the output of NPFSS, PFSS, DFFSS, URBFSS, MDFFSS, MRBFSS, MG, MNP, UG or UNP, or the sites_clusters_LH/sites_clusters_W/sites_clusters_P/sites_clusters_R of the MPFSS.
##' @param sites_coord numeric matrix. Coordinates of the sites, in the same order that the data for the cluster detection.
##' @param system character. System in which the coordinates are expressed: "Euclidean" or "WGS84".
##' @param system_conv character. System to convert the coordinates for the plot. Only considered if system is "WGS84". Must be entered as in the PROJ.4 documentation
##' @param colors character. Colors of the clusters. If length(colors)=1 all the clusters will be in this color. Else it should be a vector of length the number of clusters to plot.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' coords <- coordinates(map_sites)
##'
##' res_npfss <- NPFSS(data = funi_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)
##'
##' plot_schema(output_clusters = res_npfss$sites_clusters,
##' sites_coord = coords, system = "WGS84", system_conv = "+init=epsg:2154")}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_npfss <- NPFSS(data = funi_data[indices,], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 99)
##' if(length(res_npfss$sites_clusters)>0){
##' plot_schema(output_clusters = res_npfss$sites_clusters,
##' sites_coord = coords, system = "WGS84", system_conv = "+init=epsg:2154")
##' }
##'
##' }
##'
##' @return No value returned, plots a schema of the sites and the clusters.
##'
##' @export
##'
plot_schema <- function(output_clusters, sites_coord, system, system_conv = NULL, colors = "red"){
  if(length(system)!=1){
    stop("Only one system must be specified")
  }
  if(is.null(system)){
    stop("Specify a correct system: Euclidean or WGS84")
  }
  if(system != "Euclidean" & system != "WGS84"){
    stop("Specify a correct system: Euclidean or WGS84")
  }
  if(ncol(sites_coord)!=2){
    stop("sites_coord must be a matrix with two columns")
  }
  if(system == "WGS84"){
    if(!class(system_conv)[1] == "character"){
      stop("system_conv must be a character")
    }else{
      new_coords <- suppressWarnings(rgdal::project(sites_coord, system_conv))
    }
  }else{
    new_coords <- sites_coord
  }


  indices <- c(1:length(output_clusters))

  if(length(colors)==1){
    if("black" %in% colors){
      stop("Black is not available for clusters since it is the color of the other sites")
    }
    plot(x = new_coords[,1], y = new_coords[,2], col = "black", xlim = range(new_coords[,1]), ylim = range(new_coords[,2]), asp = 1, pch = 19, xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', bty = 'n')
    for(indice in indices){
      points(x = new_coords[output_clusters[[indice]],1], y = new_coords[output_clusters[[indice]],2], col = colors, pch = 19)

      barycenter <- colMeans(matrix(unique(new_coords[output_clusters[[indice]],], MARGIN = 1), ncol = 2))
      TeachingDemos::shadowtext(x = barycenter[1], y = barycenter[2], labels = indice, bg = "white", col = "black")

    }
  }else{
    if(length(colors)!=length(indices)){
      stop("There must be one color or the same number of colors than the desired clusters to be plotted")
    }
    if("black" %in% colors){
      stop("Black is not available for clusters since it is the color of the other sites")
    }
    plot(x = new_coords[,1], y = new_coords[,2], col = "black", xlim = range(new_coords[,1]), ylim = range(new_coords[,2]), asp = 1, pch = 19, xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', bty = 'n')
    for(i in 1:length(indices)){
      points(x = new_coords[output_clusters[[indices[i]]],1], y = new_coords[output_clusters[[indices[i]]],2], col = colors[i], pch = 19)

      barycenter <- colMeans(matrix(unique(new_coords[output_clusters[[indices[i]]],], MARGIN = 1), ncol = 2))
      TeachingDemos::shadowtext(x = barycenter[1], y = barycenter[2], labels = indices[i], bg = "white", col = "black")

    }
  }
}

########################################################################################################################################################
##' @title Map of circular clusters
##'
##' @description This function plots a map of the sites and the circular clusters.
##'
##' @param spobject SpObject. SpatialObject with the same coordinates system that centres (the same that sites_coord in the scan functions)
##' @param centres numeric matrix or vector if only one cluster was detected. Coordinates of the centres of each cluster.
##' @param radius numeric vector. Radius of each cluster in the user units if system = "Euclidean", or in km if system = "WGS84" (in the output of the scan functions)
##' @param system character. System in which the coordinates are expressed: "Euclidean" or "WGS84".
##' @param colors character. Colors of the clusters. If length(colors)=1 all the clusters will be in this color. Else it should be a vector of length the number of clusters to plot.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' coords <- coordinates(map_sites)
##'
##' res_npfss <- NPFSS(data = funi_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)
##'
##'
##' plot_map(spobject = map_sites, centres = res_npfss$centres_clusters,
##' radius = res_npfss$radius_clusters, system = "WGS84")}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_npfss <- NPFSS(data = funi_data[indices,], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 99)
##' if(length(res_npfss$sites_clusters)>0){
##' plot_map(spobject = map_sites[indices,], centres = res_npfss$centres_clusters,
##' radius = res_npfss$radius_clusters, system = "WGS84")
##' }
##'
##' }
##'
##' @return No value returned, plots a map of the sites and the circular clusters.
##'
##' @export
##'
plot_map <- function(spobject, centres, radius, system, colors = "red"){

  if(length(system)!=1){
    stop("Only one system must be specified")
  }
  if(is.null(system)){
    stop("Specify a correct system: Euclidean or WGS84")
  }
  if(system != "Euclidean" & system != "WGS84"){
    stop("Specify a correct system: Euclidean or WGS84")
  }

  centres <- matrix(centres, ncol = 2)
  if(system == "Euclidean"){
    plot(spobject, asp = 1)
  }else{
    plot(spobject)
  }

  points(coordinates(spobject), pch = 16, cex = 0.5)
  if(length(radius)!=nrow(centres)){
    stop("radius must have the same number of elements that the number of rows in centres")
  }

  if(length(colors)==1){
    for(i in 1:length(radius)){
      if(radius[i]>0){
        if(system == "Euclidean"){
          plotrix::draw.circle(x = centres[i,1], y = centres[i,2], radius = radius[i], border = colors)
        }else{
          temp_circle <- swfscMisc::circle.polygon(x = centres[i,1], y = centres[i,2], radius = radius[i], units = "km", ellipsoid = datum(model = c("wgs84")), destination.type = "ellipsoid", poly.type = "gc.earth")
          lines(temp_circle, col = colors, lwd = 2)
        }
      }else{
        points(x = centres[i,1], y = centres[i,2], col = colors, pch = 16, cex = 0.7)
      }

      TeachingDemos::shadowtext(x = centres[i,1], y = centres[i,2], labels = i, bg = "white", col = "black")
    }
  }else{
    if(length(colors)!=length(radius)){
      stop("There must be one color or the same number of colors than the number of clusters to be plotted")
    }
    for(i in 1:length(radius)){
      if(radius[i]>0){
        if(system == "Euclidean"){
          plotrix::draw.circle(x = centres[i,1], y = centres[i,2], radius = radius[i], border = colors[i])
        }else{
          temp_circle <- swfscMisc::circle.polygon(x = centres[i,1], y = centres[i,2], radius = radius[i], units = "km", ellipsoid = datum(model = c("wgs84")), destination.type = "ellipsoid", poly.type = "gc.earth")
          lines(temp_circle, col = colors[i], lwd = 2)
        }
        }else{
          points(x = centres[i,1], y = centres[i,2], col = colors[i], pch = 16, cex = 0.7)
        }
      TeachingDemos::shadowtext(x = centres[i,1], y = centres[i,2], labels = i, bg = "white", col = "black")
    }
  }

}

########################################################################################################################################################
##' @title Map of the clusters
##'
##' @description This function plots a map of the sites and the clusters
##'
##' @param spobject SpObject. SpatialObject corresponding the sites.
##' @param sites_coord numeric matrix. Coordinates of the sites or the individuals, in the same order that the data for the cluster detection.
##' @param output_clusters list. List of the sites in the clusters: it is the sites_clusters of the output of NPFSS, PFSS, DFFSS, URBFSS, MDFFSS, MRBFSS, MG, MNP, UG or UNP, or the sites_clusters_LH/sites_clusters_W/sites_clusters_P/sites_clusters_R of the MPFSS.
##' @param system character. System in which the coordinates are expressed: "Euclidean" or "WGS84".
##' @param colors character. Colors of the clusters. If length(colors)=1 all the clusters will be in this color. Else it should be a vector of length the number of clusters to plot.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' coords <- coordinates(map_sites)
##'
##' res_npfss <- NPFSS(data = funi_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)
##'
##' plot_map2(spobject = map_sites, sites_coord = coords,
##' output_clusters = res_npfss$sites_clusters, system = "WGS84")}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_npfss <- NPFSS(data = funi_data[indices,], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 99)
##' if(length(res_npfss$sites_clusters)>0){
##' plot_map2(spobject = map_sites[indices,], sites_coord = coords,
##' output_clusters = res_npfss$sites_clusters, system = "WGS84")
##' }
##'
##' }
##'
##' @return No value returned, plots a map of the sites and the clusters.
##'
##' @export
##'
plot_map2 <- function(spobject, sites_coord, output_clusters, system, colors = "red"){

  if(length(system)!=1){
    stop("Only one system must be specified")
  }
  if(is.null(system)){
    stop("Specify a correct system: Euclidean or WGS84")
  }
  if(system != "Euclidean" & system != "WGS84"){
    stop("Specify a correct system: Euclidean or WGS84")
  }

  if(system == "Euclidean"){
    plot(spobject, asp = 1)
  }else{
    plot(spobject)
  }

  for(cl in 1:length(output_clusters)){
    # if the sites_coord are in fact individuals we have to get the real sites coordinates
    coords <- unique(sites_coord[output_clusters[[cl]],], MARGIN = 1)
    coords <- matrix(coords, ncol = 2)
    # now we have to identify the sites
    coords_all <- coordinates(spobject)
    indices_sites <- numeric(nrow(coords))

    for(site in 1:nrow(coords)){
      indices_sites[site] <- which(coords_all[,1] == coords[site,1] & coords_all[,2] == coords[site,2])
    }

    if(length(colors)==1){
      c_union <- suppressWarnings(gUnaryUnion(spobject[indices_sites,]))
      adj_color <- adjustcolor(colors, alpha.f = 0.6)
      plot(c_union, add = TRUE, col = adj_color)

      barycenter <- colMeans(matrix(unique(sites_coord[output_clusters[[cl]],], MARGIN = 1), ncol = 2))
      TeachingDemos::shadowtext(x = barycenter[1], y = barycenter[2], labels = cl, bg = "white", col = "black")


    }else{
      if(length(colors)!=length(output_clusters)){
        stop("There must be one color or the same number of colors than the number of clusters to be plotted")
      }

      c_union <- suppressWarnings(gUnaryUnion(spobject[indices_sites,]))
      adj_color <- adjustcolor(colors[cl], alpha.f = 0.6)

      plot(c_union, add = TRUE, col = adj_color)

      barycenter <- colMeans(matrix(unique(sites_coord[output_clusters[[cl]],], MARGIN = 1), ncol = 2))
      TeachingDemos::shadowtext(x = barycenter[1], y = barycenter[2], labels = cl, bg = "white", col = "black")

    }

  }

}
