################################################################
##' @title UNP scan procedure
##'
##' @description This function returns the significant clusters and their associated p-value for the UNP
##'
##' @param data vector. Vector of the data, each element corresponds to a site (or an individual if the observations are by individuals and not by sites).
##' @param sites_coord numeric matrix. Coordinates of the sites (or the individuals, in that case there can be many individuals with the same coordinates). If system = "WGS84" the first column corresponds to the longitude and the second to the latitude. If system = "Euclidean", the first column corresponds to the x and the second one to the y
##' @param system character. System in which the coordinates are expressed: "Euclidean" or "WGS84".
##' @param mini numeric. A minimum for the clusters (see type_minimaxi). Changing the default value may bias the inference.
##' @param maxi numeric. A Maximum for the clusters (see type_minimaxi). Changing the default value may bias the inference.
##' @param type_minimaxi character. Type of minimum and maximum: by default "sites/indiv": the mini and maxi are on the number of sites or individuals in the potential clusters. Other possible values are "area": the minimum and maximum area of the clusters, or "radius": the minimum and maximum radius.
##' @param mini_post numeric. A minimum to filter the significant clusters a posteriori (see type_minimaxi_post). The default NULL is for no filtering with a a posteriori minimum.
##' @param maxi_post numeric. A maximum to filter the significant clusters a posteriori (see type_minimaxi_post). The default NULL is for no filtering with a a posteriori maximum.
##' @param type_minimaxi_post character. Type of minimum and maximum a posteriori: by default "sites/indiv": the mini_post and maxi_post are on the number of sites or individuals in the significant clusters. Other possible values are "area": the minimum and maximum area of the clusters, or "radius": the minimum and maximum radius.
##' @param sites_areas numeric vector. Areas of the sites. It must contain the same number of elements than the rows of sites_coord. If the data is on individuals and not on sites, there can be duplicated values. By default: NULL
##' @param MC numeric. Number of Monte-Carlo permutations to evaluate the statistical significance of the clusters. By default: 999.
##' @param typeI numeric. The desired type I error. A cluster will be evaluated as significant if its associated p-value is less than typeI. By default 0.05.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("multi_data")
##' uni_data <- multi_data[,1]
##' coords <- coordinates(map_sites)
##' res_unp <- UNP(data=uni_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("multi_data")
##' uni_data <- multi_data[,1]
##' coords <- coordinates(map_sites)
##' res_unp <- UNP(data=uni_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2, MC = 9)
##' }
##'
##' @return The list of the following elements:
##' \itemize{
##' \item sites_clusters: the index of the sites (or individuals) included in the significant clusters. The clusters are listed in their order of detection. The secondary clusters are defined according to Kulldorff.
##' \item pval_clusters: the associated p-values
##' \item centres_clusters: the coordinates of the centres of each cluster
##' \item radius_clusters: the radius of the clusters in km if system = "WGS84" or in the coordinates unit otherwise
##' \item areas_clusters: the areas of the clusters (in the same units as sites_areas). Only if sites_areas is not NULL.
##' \item system: the system of coordinates
##' }
##' @export
##'
##'
UNP <- function(data, sites_coord = NULL, system = NULL, mini = 1, maxi = nrow(sites_coord)/2, type_minimaxi = "sites/indiv", mini_post = NULL, maxi_post = NULL, type_minimaxi_post = "sites/indiv", sites_areas = NULL, MC = 999, typeI = 0.05){
  matrix_clusters_provided <- NULL
  area_clusters_provided <- NULL

  if(is.null(mini_post) == FALSE | is.null(maxi_post) == FALSE){
    if(!(type_minimaxi_post %in% c("sites/indiv", "area", "radius"))){
      stop("The value of type_minimaxi_post must be sites/indiv or area or radius")
    }

    if(type_minimaxi_post == "area" & is.null(matrix_clusters_provided) == FALSE){
      if(is.null(area_clusters_provided)){
        stop("You must provide area_clusters_provided for the a posteriori filtering")
      }
    }
    if(type_minimaxi_post == "radius" & is.null(matrix_clusters_provided) == FALSE){
      stop("radius a posteriori filtering is only available when matrix_clusters_provided = NULL")
    }
    if(type_minimaxi_post == "area" & is.null(matrix_clusters_provided) == TRUE){
      if(is.null(sites_areas)){
        stop("You must provide sites_areas for the a posteriori filtering")
      }
    }

    # a posteriori filtering:
    filtering_post <- TRUE
  }else{
    filtering_post <- FALSE
  }

  if(is.null(matrix_clusters_provided)){
    if(is.null(sites_coord) | is.null(system) | is.null(mini) | is.null(maxi)){
      stop("You must specify the coordinates of the sites, the system of coordinates used, and the minimum and maximum number of sites in a cluster")
    }else{
      if(! (type_minimaxi %in% c("sites/indiv", "area", "radius")) ){
        stop("The value of type_minimaxi must be sites/indiv or area or radius")
      }
      if(type_minimaxi == "area" & is.null(sites_areas)){
        stop("If type_minimaxi = area you must specify the areas of the sites or individuals")
      }
      if(is.null(sites_areas)==FALSE & length(sites_areas)!=nrow(sites_coord)){
        stop("sites_areas must contain the same number of elements than the number of rows in sites_coord")
      }
      details_clusters <- clusters(sites_coord, system, mini, maxi, type_minimaxi, sites_areas)
      matrix_clusters <- details_clusters$matrix_clusters
      centres <- details_clusters$centres
      radius <- details_clusters$radius
      system <- details_clusters$system
      if(is.null(sites_areas)==FALSE){
        areas <- details_clusters$areas
      }
      details_clusters <- NULL
    }
  }else{
    matrix_clusters <- matrix_clusters_provided
    if(sum(!(matrix_clusters %in% c(0,1)))>0){
      stop("The matrix of potential clusters must contain only 0 and 1 values")
    }
    if(is.null(area_clusters_provided) == FALSE & length(area_clusters_provided) != nrow(matrix_clusters_provided)){
      stop("Length of area_clusters_provided must be the same as the number of columns of matrix_clusters_provided")
    }
    areas <- area_clusters_provided
  }

  if(is.vector(data)){
    if(nrow(matrix_clusters)!=length(data)){
      stop("The number of sites must be the same in the matrix of clusters and in the data")
    }

    rank_data <- rank(data, ties.method = "average")
    index <- wmw_uni(matrix(rank_data, ncol = 1), matrix_clusters)

    nb_clusters <- ncol(matrix_clusters)

    nb_sites <- nrow(matrix_clusters)
    permutations <- permutate(1:nb_sites, MC)
    generation_signif <- sapply(1:MC, function(g) rank_data[permutations[g,]])

    results <- wmw_uni(matrix(generation_signif, ncol = MC), matrix_clusters)

    stat_MC <- rowMaxs(results)
    pvals <- sapply(1:nb_clusters, function(j) (length(which(stat_MC >= index[1,j]))+1)/(MC+1))

    index_clusters_temp <- which(pvals <= typeI)

    if(length(index_clusters_temp)>0){
      ordre <- order(index[index_clusters_temp], decreasing = TRUE)
      index_clusters_temp <- index_clusters_temp[ordre]

      # a posteriori filtering
      if(filtering_post == TRUE){
        if(type_minimaxi_post == "sites/indiv"){
          index_clusters <- post_filt_nb_sites(mini_post, maxi_post, nb_sites, index_clusters_temp, matrix_clusters)
        }
        if(type_minimaxi_post == "radius"){
          index_clusters <- post_filt_radius(mini_post, maxi_post, radius, index_clusters_temp)
        }
        if(type_minimaxi_post == "area"){

          if(is.null(matrix_clusters_provided)){
            areas_clusters <- areas
          }else{
            areas_clusters <- area_clusters_provided
          }

          index_clusters <- post_filt_area(mini_post, maxi_post, areas_clusters, index_clusters_temp)

        }
      }else{
        index_clusters <- index_clusters_temp
      }
    }else{
      index_clusters <- index_clusters_temp
    }


    # non overlapping clusters:
    final_index <- non_overlap(index_clusters, matrix_clusters)

    pval_clusters <- pvals[final_index]
    sites_clusters <- lapply(final_index, function(j) which(matrix_clusters[,j]==1))
    if(is.null(matrix_clusters_provided)){
      centres_clusters <- matrix(centres[final_index,], ncol = 2)
      radius_clusters <- radius[final_index]
      if(is.null(sites_areas) == FALSE){
        areas_clusters <- areas[final_index]
      }
    }else{
      if(is.null(area_clusters_provided) == FALSE){
        areas_clusters <- areas[final_index]
      }
    }

  }else{
    stop("The data must be a vector")
    }
  if(is.null(matrix_clusters_provided)){
    if(is.null(sites_areas)){
      return(list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, centres_clusters = centres_clusters, radius_clusters = radius_clusters, system = system))
    }else{
      return(list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, centres_clusters = centres_clusters, radius_clusters = radius_clusters, areas_clusters = areas_clusters, system = system))
    }
  }else{
    if(is.null(area_clusters_provided)){
      return(list(sites_clusters = sites_clusters, pval_clusters = pval_clusters))
    }else{
      return(list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, areas_clusters = areas_clusters))
    }
  }
}


################################################################
##' @title UG scan procedure
##'
##' @description This function returns the significant clusters and their associated p-value for the UG.
##'
##' @param data vector. Vector of the data, each element corresponds to a site (or an individual if the observations are by individuals and not by sites).
##' @param sites_coord numeric matrix. Coordinates of the sites (or the individuals, in that case there can be many individuals with the same coordinates).
##' @param system character. System in which the coordinates are expressed: "Euclidean" or "WGS84".
##' @param mini numeric. A minimum for the clusters (see type_minimaxi). Changing the default value may bias the inference.
##' @param maxi numeric. A Maximum for the clusters (see type_minimaxi). Changing the default value may bias the inference.
##' @param type_minimaxi character. Type of minimum and maximum: by default "sites/indiv": the mini and maxi are on the number of sites or individuals in the potential clusters. Other possible values are "area": the minimum and maximum area of the clusters, or "radius": the minimum and maximum radius.
##' @param mini_post numeric. A minimum to filter the significant clusters a posteriori (see type_minimaxi_post). The default NULL is for no filtering with a a posteriori minimum.
##' @param maxi_post numeric. A maximum to filter the significant clusters a posteriori (see type_minimaxi_post). The default NULL is for no filtering with a a posteriori maximum.
##' @param type_minimaxi_post character. Type of minimum and maximum a posteriori: by default "sites/indiv": the mini_post and maxi_post are on the number of sites or individuals in the significant clusters. Other possible values are "area": the minimum and maximum area of the clusters, or "radius": the minimum and maximum radius.
##' @param sites_areas numeric vector. Areas of the sites. It must contain the same number of elements than the rows of sites_coord. If the data is on individuals and not on sites, there can be duplicated values. By default: NULL
##' @param MC numeric. Number of Monte-Carlo permutations to evaluate the statistical significance of the clusters. By default: 999.
##' @param typeI numeric. The desired type I error. A cluster will be evaluated as significant if its associated p-value is less than typeI. By default 0.05.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("multi_data")
##' uni_data <- multi_data[,1]
##' coords <- coordinates(map_sites)
##' res_ug <- UG(data=uni_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("multi_data")
##' uni_data <- multi_data[,1]
##' coords <- coordinates(map_sites)
##' res_ug <- UG(data=uni_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2, MC = 9)}
##'
##'
##'
##' @return The list of the following elements:
##' \itemize{
##' \item sites_clusters: the index of the sites (or the individuals) included in the significant clusters. The clusters are listed in their order of detection. The secondary clusters are defined according to Kulldorff.
##' \item pval_clusters: the associated p-values
##' \item centres_clusters: the coordinates of the centres of each cluster
##' \item radius_clusters: the radius of the clusters in km if system = WGS84 or in the user units
##' \item areas_clusters: the areas of the clusters (in the same units as sites_areas). Only if sites_areas is not NULL.
##' \item system: the system of coordinates
##' }
##' @export
##'
##'
UG <- function(data, sites_coord = NULL, system = NULL, mini = 1, maxi = nrow(sites_coord)/2, type_minimaxi = "sites/indiv", mini_post = NULL, maxi_post = NULL, type_minimaxi_post = "sites/indiv", sites_areas = NULL, MC = 999, typeI = 0.05){

  matrix_clusters_provided <- NULL
  area_clusters_provided <- NULL

  if(is.null(mini_post) == FALSE | is.null(maxi_post) == FALSE){
    if(!(type_minimaxi_post %in% c("sites/indiv", "area", "radius"))){
      stop("The value of type_minimaxi_post must be sites/indiv or area or radius")
    }

    if(type_minimaxi_post == "area" & is.null(matrix_clusters_provided) == FALSE){
      if(is.null(area_clusters_provided)){
        stop("You must provide area_clusters_provided for the a posteriori filtering")
      }
    }
    if(type_minimaxi_post == "radius" & is.null(matrix_clusters_provided) == FALSE){
      stop("radius a posteriori filtering is only available when matrix_clusters_provided = NULL")
    }
    if(type_minimaxi_post == "area" & is.null(matrix_clusters_provided) == TRUE){
      if(is.null(sites_areas)){
        stop("You must provide sites_areas for the a posteriori filtering")
      }
    }

    # a posteriori filtering:
    filtering_post <- TRUE
  }else{
    filtering_post <- FALSE
  }

  if(is.null(matrix_clusters_provided)){
    if(is.null(sites_coord) | is.null(system) | is.null(mini) | is.null(maxi)){
      stop("You must specify the coordinates of the sites, the system of coordinates used, and the minimum and maximum number of sites in a cluster")
    }else{
      if(! (type_minimaxi %in% c("sites/indiv", "area", "radius")) ){
        stop("The value of type_minimaxi must be sites/indiv or area or radius")
      }
      if(type_minimaxi == "area" & is.null(sites_areas)){
        stop("If type_minimaxi = area you must specify the areas of the sites or individuals")
      }
      if(is.null(sites_areas)==FALSE & length(sites_areas)!=nrow(sites_coord)){
        stop("sites_areas must contain the same number of elements than the number of rows in sites_coord")
      }
      details_clusters <- clusters(sites_coord, system, mini, maxi, type_minimaxi, sites_areas)
      matrix_clusters <- details_clusters$matrix_clusters
      centres <- details_clusters$centres
      radius <- details_clusters$radius
      system <- details_clusters$system
      if(is.null(sites_areas)==FALSE){
        areas <- details_clusters$areas
      }
      details_clusters <- NULL
    }
  }else{
    matrix_clusters <- matrix_clusters_provided
    if(sum(!(matrix_clusters %in% c(0,1)))>0){
      stop("The matrix of potential clusters must contain only 0 and 1 values")
    }
    if(is.null(area_clusters_provided) == FALSE & length(area_clusters_provided) != nrow(matrix_clusters_provided)){
      stop("Length of area_clusters_provided must be the same as the number of columns of matrix_clusters_provided")
    }
    areas <- area_clusters_provided
  }

  if(is.vector(data)){
    if(nrow(matrix_clusters)!=length(data)){
      stop("The number of sites must be the same in the matrix of clusters and in the data")
    }

    index <- dfree(matrix(data, ncol = 1), matrix_clusters)

    nb_clusters <- ncol(matrix_clusters)
    nb_sites <- nrow(matrix_clusters)
    permutations <- permutate(1:nb_sites, MC)
    generation_signif <- sapply(1:MC, function(g) data[permutations[g,]])

    results <- dfree(matrix(generation_signif, ncol = MC), matrix_clusters)

    stat_MC <- rowMaxs(results)
    pvals <- sapply(1:nb_clusters, function(j) (length(which(stat_MC >= index[1,j]))+1)/(MC+1))

    index_clusters_temp <- which(pvals <= typeI)

    if(length(index_clusters_temp) > 0){
      ordre <- order(index[index_clusters_temp], decreasing = TRUE)
      index_clusters_temp <- index_clusters_temp[ordre]

      # a posteriori filtering
      if(filtering_post == TRUE){
        if(type_minimaxi_post == "sites/indiv"){
          index_clusters <- post_filt_nb_sites(mini_post, maxi_post, nb_sites, index_clusters_temp, matrix_clusters)
        }
        if(type_minimaxi_post == "radius"){
          index_clusters <- post_filt_radius(mini_post, maxi_post, radius, index_clusters_temp)
        }
        if(type_minimaxi_post == "area"){
          if(is.null(matrix_clusters_provided)){
            areas_clusters <- areas
          }else{
            areas_clusters <- area_clusters_provided
          }

          index_clusters <- post_filt_area(mini_post, maxi_post, areas_clusters, index_clusters_temp)

        }
      }else{
        index_clusters <- index_clusters_temp
      }

    }else{
      index_clusters <- index_clusters_temp
    }



    # non overlapping clusters:
    final_index <- non_overlap(index_clusters, matrix_clusters)

    pval_clusters <- pvals[final_index]
    sites_clusters <- lapply(final_index, function(j) which(matrix_clusters[,j]==1))
    if(is.null(matrix_clusters_provided)){
      centres_clusters <- matrix(centres[final_index,], ncol = 2)
      radius_clusters <- radius[final_index]
      if(is.null(sites_areas) == FALSE){
        areas_clusters <- areas[final_index]
      }
    }else{
      if(is.null(area_clusters_provided) == FALSE){
        areas_clusters <- areas[final_index]
      }
    }
  }else{
    stop("The data must be a vector")
  }

  if(is.null(matrix_clusters_provided)){
    if(is.null(sites_areas)){
      return(list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, centres_clusters = centres_clusters, radius_clusters = radius_clusters, system = system))
    }else{
      return(list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, centres_clusters = centres_clusters, radius_clusters = radius_clusters, areas_clusters = areas_clusters, system = system))
    }
  }else{
    if(is.null(area_clusters_provided)){
      return(list(sites_clusters = sites_clusters, pval_clusters = pval_clusters))
    }else{
      return(list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, areas_clusters = areas_clusters))
    }
  }

}
