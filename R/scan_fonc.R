################################################################
##' @title NPFSS scan procedure (univariate functional or multivariate functional)
##'
##' @description This function returns the significant clusters and their associated p-value for the NPFSS (multivariate or univariate functional data).
##'
##' @param data list of numeric matrices or a matrix. List of nb_sites (or nb_individuals if the observations are by individuals and not by site) matrices of the data, the rows correspond to the variables and each column represents an observation time (multivariate case) ; or Matrix of the data, the rows correspond to the sites (or to the individuals) and each column represents an observation time (univariate case). The times must be equally spaced and the same for each site/individual.
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
##' @param nbCPU numeric. Number of CPU. If nbCPU > 1 parallelization is done. By default: 1.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' coords <- coordinates(map_sites)
##'
##' res_npfss <- NPFSS(data = funi_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_npfss <- NPFSS(data = funi_data[indices,], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = 1, MC = 99)
##' }
##'
##' @return The list of the following elements:
##' \itemize{
##' \item sites_clusters: the index of the sites (or individuals) included in the significant clusters. The clusters are listed in their order of detection. The secondary clusters are defined according to Kulldorff.
##' \item pval_clusters: the associated p-values
##' \item centres_clusters: the coordinates of the centres of each cluster
##' \item radius_clusters: the radius of the clusters in km if system = WGS84 or in the user units
##' \item areas_clusters: the areas of the clusters (in the same units as sites_areas). Only if sites_areas is not NULL.
##' \item system: the system of coordinates
##' }
##' @export
##'
##'
NPFSS <- function(data, sites_coord = NULL, system = NULL, mini = 1, maxi = nrow(sites_coord)/2, type_minimaxi = "sites/indiv", mini_post = NULL, maxi_post = NULL, type_minimaxi_post = "sites/indiv", sites_areas = NULL, MC = 999, typeI = 0.05, nbCPU = 1){

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

  # univariate or multivariate ?
  if(class(data)[1]=="list"){
    # multivariate
    if(length(data)!=nrow(matrix_clusters)){
      stop("The data must contain the same number of sites than the matrix of clusters")
    }
    nb_rows <- sapply(1:length(data), function(r) nrow(data[[r]]))
    nb_cols <- sapply(1:length(data), function(r) ncol(data[[r]]))

    if(sum((nb_rows!=nb_rows[1])) + sum((nb_cols!=nb_cols[1])) > 0){
      stop("All matrices of the list must be of the same dimensions")
    }
    if(nb_rows[1]==1){
      stop("There must be at least two variables (rows)")
    }
    if(nb_cols[1]==1){
      stop("There must be at least two observation times (columns)")
    }

    signs <- multi_signs_matrix(data)
    index <- multi_fWMW(signs, matrix_clusters)

    nb_clusters <- ncol(matrix_clusters)
    nb_sites <- nrow(matrix_clusters)
    permutations <- permutate(1:nb_sites, MC)
    generation_signif <- list()
    for(g in 1:MC){
      generation_signif[[g]] <- list()
      for(s in 1:nb_sites){
        generation_signif[[g]][[s]] <- signs[[permutations[g,s]]]
      }
    }

    num_cores <- detectCores()

    if(nbCPU <=1 | num_cores <2){
      results <- lapply(1:MC, function(i) multi_fWMW(generation_signif[[i]], matrix_clusters))
    }else{

      if(num_cores < nbCPU){
        nbCPU <- num_cores
      }

      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(nbCPU)
        results <- parLapply(cl, 1:MC, function(i) multi_fWMW(generation_signif[[i]], matrix_clusters))
        stopCluster(cl)

      }else{
        results <- mclapply(1:MC, function(i) multi_fWMW(generation_signif[[i]], matrix_clusters), mc.cores = nbCPU)
      }

    }
    results <- matrix(unlist(results), ncol = nb_clusters, byrow = TRUE)
    stat_MC <- rowMaxs(results)
    pvals <- sapply(1:nb_clusters, function(j) (length(which(stat_MC >= index[j]))+1)/(MC+1))

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
    if(class(data)[1]=="matrix"){
      # univariate
      if(nrow(matrix_clusters)!=nrow(data)){
        stop("The number of sites must be the same in the matrix of clusters and in the data")
      }
      if(ncol(data)==1){
        stop("There must be at least two observation times")
      }

      signs <- uni_signs_matrix(data)
      index <- uni_fWMW(signs, matrix_clusters)

      nb_clusters <- ncol(matrix_clusters)
      nb_sites <- nrow(matrix_clusters)
      permutations <- permutate(1:nb_sites, MC)

      generation_signif <- lapply(1:MC, function(g) signs[permutations[g,],])

      num_cores <- detectCores()

      if(nbCPU <=1 | num_cores <2){
        results <- lapply(1:MC, function(i) uni_fWMW(generation_signif[[i]], matrix_clusters))
      }else{

        if(num_cores < nbCPU){
          nbCPU <- num_cores
        }

        if(Sys.info()["sysname"] == "Windows"){
          cl <- makeCluster(nbCPU)
          results <- parLapply(cl, 1:MC, function(i) uni_fWMW(generation_signif[[i]], matrix_clusters))
          stopCluster(cl)

        }else{
          results <- mclapply(1:MC, function(i) uni_fWMW(generation_signif[[i]], matrix_clusters), mc.cores = nbCPU)
        }

      }
      results <- matrix(unlist(results), ncol = nb_clusters, byrow = TRUE)
      stat_MC <- rowMaxs(results)
      pvals <- sapply(1:nb_clusters, function(j) (length(which(stat_MC >= index[j]))+1)/(MC+1))

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
      stop("The data must be a list (multivariate case) or a matrix (univariate case)")
    }
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
##' @title PFSS scan procedure
##'
##' @description This function returns the significant clusters and their associated p-value for the PFSS.
##'
##' @param data matrix. Matrix of the data, the rows correspond to the sites (or to the individuals if the observations are by individuals and not by sites) and each column represents an observation time. The times must be equally spaced and the same for each site/individual.
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
##' @param nbCPU numeric. Number of CPU. If nbCPU > 1 parallelization is done. By default: 1.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' coords <- coordinates(map_sites)
##'
##' res_pfss <- PFSS(data = funi_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_pfss <- PFSS(data = funi_data[indices,], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 9)
##' }
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
PFSS <- function(data, sites_coord = NULL, system = NULL, mini = 1, maxi = nrow(sites_coord)/2, type_minimaxi = "sites/indiv", mini_post = NULL, maxi_post = NULL, type_minimaxi_post = "sites/indiv", sites_areas = NULL, MC = 999, typeI = 0.05, nbCPU = 1){

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

  if(class(data)[1]=="matrix"){
    if(nrow(matrix_clusters)!=nrow(data)){
      stop("The number of sites must be the same in the matrix of clusters and in the data")
    }
    if(ncol(data)==1){
      stop("There must be at least two observation times")
    }

    index <- fanova_cpp(data, matrix_clusters)

    nb_clusters <- ncol(matrix_clusters)
    nb_sites <- nrow(matrix_clusters)
    permutations <- permutate(1:nb_sites, MC)

    generation_signif <- lapply(1:MC, function(g) data[permutations[g,],])

    num_cores <- detectCores()

    if(nbCPU <=1 | num_cores <2){
      results <- lapply(1:MC, function(i) fanova_cpp(generation_signif[[i]], matrix_clusters))
    }else{
      if(num_cores < nbCPU){
        nbCPU <- num_cores
      }

      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(nbCPU)
        results <- parLapply(cl, 1:MC, function(i) fanova_cpp(generation_signif[[i]], matrix_clusters))
        stopCluster(cl)

      }else{
        results <- mclapply(1:MC, function(i) fanova_cpp(generation_signif[[i]], matrix_clusters), mc.cores = nbCPU)
      }
    }
    results <- matrix(unlist(results), ncol = nb_clusters, byrow = TRUE)
    stat_MC <- rowMaxs(results)
    pvals <- sapply(1:nb_clusters, function(j) (length(which(stat_MC >= index[j]))+1)/(MC+1))

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
    stop("The data must be a matrix")
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
##' @title MPFSS scan procedure
##'
##' @description This function returns the significant clusters and their associated p-value for the MPFSS.
##'
##' @param data list of numeric matrices. List of nb_sites (or nb_individuals if the observations are by individuals and not by sites) matrices of the data, the rows correspond to the variables and each column represents an observation time. The times must be equally spaced and the same for each site/individual.
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
##' @param method character vector. The methods to compute the significant clusters. Options: "LH", "W", "P", "R" for respectively the Lawley-Hotelling trace test statistic, The Wilks lambda test statistic, the Pillai trace test statistic and the Roy's maximum root test statistic. By default all are computed.
##' @param nbCPU numeric. Number of CPU. If nbCPU > 1 parallelization is done. By default: 1.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("fmulti_data")
##' coords <- coordinates(map_sites)
##'
##' res_mpfss <- MPFSS(data = fmulti_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("fmulti_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_mpfss <- MPFSS(data = fmulti_data[indices], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 9)
##' }
##'
##'
##' @return The list of the following elements according to the methods wanted:
##' \itemize{
##' \item sites_clusters_LH: the index of the sites (or individuals) included in the significant clusters for the LH method. The clusters are listed in their order of detection. The secondary clusters are defined according to Kulldorff.
##' \item pval_clusters_LH: the associated p-values
##' \item centres_clusters_LH: the coordinates of the centres of each cluster
##' \item radius_clusters_LH: the radius of the clusters in km if system = WGS84 or in the user units
##' \item areas_clusters_LH: the areas of the clusters (in the same units as sites_areas). Only if sites_areas is not NULL.
##' \item sites_clusters_W: the index of the sites (or individuals) included in the significant clusters for the W method. The clusters are listed in their order of detection. The secondary clusters are defined according to Kulldorff.
##' \item pval_clusters_W: the associated p-values
##' \item centres_clusters_W: the coordinates of the centres of each cluster
##' \item radius_clusters_W: the radius of the clusters in km if system = WGS84 or in the user units
##' \item areas_clusters_W: the areas of the clusters (in the same units as sites_areas). Only if sites_areas is not NULL.
##' \item sites_clusters_P: the index of the sites (or individuals) included in the significant clusters for the P method. The clusters are listed in their order of detection. The secondary clusters are defined according to Kulldorff.
##' \item pval_clusters_P: the associated p-values
##' \item centres_clusters_P: the coordinates of the centres of each cluster
##' \item radius_clusters_P: the radius of the clusters in km if system = WGS84 or in the user units
##' \item areas_clusters_P: the areas of the clusters (in the same units as sites_areas). Only if sites_areas is not NULL.
##' \item sites_clusters_R: the index of the sites (or individuals) included in the significant clusters for the R method. The clusters are listed in their order of detection. The secondary clusters are defined according to Kulldorff.
##' \item pval_clusters_R: the associated p-values
##' \item centres_clusters_R: the coordinates of the centres of each cluster
##' \item radius_clusters_R: the radius of the clusters in km if system = WGS84 or in the user units
##' \item areas_clusters_R: the areas of the clusters (in the same units as sites_areas). Only if sites_areas is not NULL.
##' \item system: the system of coordinates
##' }
##' @export
##'
##'
MPFSS <- function(data, sites_coord = NULL, system = NULL, mini = 1, maxi = nrow(sites_coord)/2, type_minimaxi = "sites/indiv", mini_post = NULL, maxi_post = NULL, type_minimaxi_post = "sites/indiv", sites_areas = NULL, MC = 999, typeI = 0.05, method = c("LH","W","P","R"), nbCPU = 1){

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

  if(class(data)[1]=="list"){
    if(length(data)!=nrow(matrix_clusters)){
      stop("The data must contain the same number of sites than the matrix of clusters")
    }
    nb_rows <- sapply(1:length(data), function(r) nrow(data[[r]]))
    nb_cols <- sapply(1:length(data), function(r) ncol(data[[r]]))

    if(sum((nb_rows!=nb_rows[1])) + sum((nb_cols!=nb_cols[1])) > 0){
      stop("All matrices of the list must be of the same dimensions")
    }
    if(nb_rows[1]==1){
      stop("There must be at least two variables (rows)")
    }
    if(nb_cols[1]==1){
      stop("There must be at least two observation times (columns)")
    }
    if(sum(!(method %in% c("LH","W","P","R")))>0){
      stop("The methods must be among LH, W, P and R")
    }

    nb_clusters <- ncol(matrix_clusters)
    nb_sites <- nrow(matrix_clusters)
    nb_times <- ncol(data[[1]])
    nb_var <- nrow(data[[1]])

    mean_tot <- matrix(0, nrow = nb_var, ncol = nb_times)
    for(i in 1:nb_sites){
      mean_tot <- mean_tot + data[[i]]
    }
    mean_tot <- mean_tot / nb_sites

    cst1 <- matrix(0, nrow = nb_var, ncol = nb_var)
    for(i in 1:nb_sites){
      cst1 <- cst1 + data[[i]] %*% t(data[[i]])
    }
    cst1 <- cst1 / nb_times
    cst2 <- mean_tot %*% t(mean_tot)
    cst2 <- cst2 * nb_sites / nb_times

    index <- fmanova_cpp(data, matrix_clusters, cst1, cst2)
    index_LH <- index$LH
    index_W <- index$W
    index_P <- index$P
    index_R <- index$R

    permutations <- permutate(1:nb_sites, MC)
    generation_signif <- list()
    for(g in 1:MC){
      generation_signif[[g]] <- list()
      for(s in 1:nb_sites){
        generation_signif[[g]][[s]] <- data[[permutations[g,s]]]
      }
    }

    num_cores <- detectCores()

    if(nbCPU <=1 | num_cores <2){
      results <- lapply(1:MC, function(i) fmanova_cpp(generation_signif[[i]], matrix_clusters, cst1, cst2))
      results <- purrr::transpose(results)
    }else{
      if(num_cores < nbCPU){
        nbCPU <- num_cores
      }

      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(nbCPU)
        results <- parLapply(cl, 1:MC, function(i) fmanova_cpp(generation_signif[[i]], matrix_clusters, cst1, cst2))
        stopCluster(cl)

      }else{
        results <- mclapply(1:MC, function(i) fmanova_cpp(generation_signif[[i]], matrix_clusters, cst1, cst2), mc.cores = nbCPU)
      }

      results <- purrr::transpose(results)
    }

    if("LH" %in% method){
      results_LH <- results$LH
      results_LH <- matrix(unlist(results_LH), ncol = nb_clusters, byrow = TRUE)
      stat_MC_LH <- rowMaxs(results_LH)
      pvals_LH <- sapply(1:nb_clusters, function(j) (length(which(stat_MC_LH >= index_LH[j]))+1)/(MC+1))
      index_clusters_LH_temp <- which(pvals_LH <= typeI)

      if(length(index_clusters_LH_temp) > 0){
        ordre_LH <- order(index_LH[index_clusters_LH_temp], decreasing = TRUE)
        index_clusters_LH_temp <- index_clusters_LH_temp[ordre_LH]

        # a posteriori filtering
        if(filtering_post == TRUE){
          if(type_minimaxi_post == "sites/indiv"){
            index_clusters_LH <- post_filt_nb_sites(mini_post, maxi_post, nb_sites, index_clusters_LH_temp, matrix_clusters)
          }
          if(type_minimaxi_post == "radius"){
            index_clusters_LH <- post_filt_radius(mini_post, maxi_post, radius, index_clusters_LH_temp)
          }
          if(type_minimaxi_post == "area"){
            if(is.null(matrix_clusters_provided)){
              areas_clusters <- areas
            }else{
              areas_clusters <- area_clusters_provided
            }

            index_clusters_LH <- post_filt_area(mini_post, maxi_post, areas_clusters, index_clusters_LH_temp)
          }
        }else{
          index_clusters_LH <- index_clusters_LH_temp
        }
      }else{
        index_clusters_LH <- index_clusters_LH_temp
      }



      # non overlapping clusters:
      final_index_LH <- non_overlap(index_clusters_LH, matrix_clusters)

      pval_clusters_LH <- pvals_LH[final_index_LH]
      sites_clusters_LH <- lapply(final_index_LH, function(j) which(matrix_clusters[,j]==1))
      if(is.null(matrix_clusters_provided)){
        centres_clusters_LH <- matrix(centres[final_index_LH,], ncol = 2)
        radius_clusters_LH <- radius[final_index_LH]
        if(is.null(sites_areas) == FALSE){
          areas_clusters_LH <- areas[final_index_LH]
        }
      }else{
        if(is.null(area_clusters_provided) == FALSE){
          areas_clusters_LH <- areas[final_index_LH]
        }
      }
    }
    if("W" %in% method){
      results_W <- results$W
      results_W <- matrix(unlist(results_W), ncol = nb_clusters, byrow = TRUE)
      stat_MC_W <- rowMins(results_W)
      pvals_W <- sapply(1:nb_clusters, function(j) (length(which(stat_MC_W <= index_W[j]))+1)/(MC+1))
      index_clusters_W_temp <- which(pvals_W <= typeI)

      if(length(index_clusters_W_temp) > 0){
        ordre_W <- order(index_W[index_clusters_W_temp], decreasing = FALSE)
        index_clusters_W_temp <- index_clusters_W_temp[ordre_W]

        # a posteriori filtering
        if(filtering_post == TRUE){
          if(type_minimaxi_post == "sites/indiv"){
            index_clusters_W <- post_filt_nb_sites(mini_post, maxi_post, nb_sites, index_clusters_W_temp, matrix_clusters)
          }
          if(type_minimaxi_post == "radius"){
            index_clusters_W <- post_filt_radius(mini_post, maxi_post, radius, index_clusters_W_temp)
          }
          if(type_minimaxi_post == "area"){
            if(is.null(matrix_clusters_provided)){
              areas_clusters <- areas
            }else{
              areas_clusters <- area_clusters_provided
            }

            index_clusters_W <- post_filt_area(mini_post, maxi_post, areas_clusters, index_clusters_W_temp)
          }
        }else{
          index_clusters_W <- index_clusters_W_temp
        }
      }else{
        index_clusters_W <- index_clusters_W_temp
      }

      # non overlapping clusters:
      final_index_W <- non_overlap(index_clusters_W, matrix_clusters)

      pval_clusters_W <- pvals_W[final_index_W]
      sites_clusters_W <- lapply(final_index_W, function(j) which(matrix_clusters[,j]==1))
      if(is.null(matrix_clusters_provided)){
        centres_clusters_W <- matrix(centres[final_index_W,], ncol = 2)
        radius_clusters_W <- radius[final_index_W]
        if(is.null(sites_areas) == FALSE){
          areas_clusters_W <- areas[final_index_W]
        }
      }else{
        if(is.null(area_clusters_provided) == FALSE){
          areas_clusters_W <- areas[final_index_W]
        }
      }
    }
    if("P" %in% method){
      results_P <- results$P
      results_P <- matrix(unlist(results_P), ncol = nb_clusters, byrow = TRUE)
      stat_MC_P <- rowMaxs(results_P)
      pvals_P <- sapply(1:nb_clusters, function(j) (length(which(stat_MC_P >= index_P[j]))+1)/(MC+1))
      index_clusters_P_temp <- which(pvals_P <= typeI)

      if(length(index_clusters_P_temp) > 0){
        ordre_P <- order(index_P[index_clusters_P_temp], decreasing = TRUE)
        index_clusters_P_temp <- index_clusters_P_temp[ordre_P]

        # a posteriori filtering
        if(filtering_post == TRUE){
          if(type_minimaxi_post == "sites/indiv"){
            index_clusters_P <- post_filt_nb_sites(mini_post, maxi_post, nb_sites, index_clusters_P_temp, matrix_clusters)
          }
          if(type_minimaxi_post == "radius"){
            index_clusters_P <- post_filt_radius(mini_post, maxi_post, radius, index_clusters_P_temp)
          }
          if(type_minimaxi_post == "area"){
            if(is.null(matrix_clusters_provided)){
              areas_clusters <- areas
            }else{
              areas_clusters <- area_clusters_provided
            }

            index_clusters_P <- post_filt_area(mini_post, maxi_post, areas_clusters, index_clusters_P_temp)

          }
        }else{
          index_clusters_P <- index_clusters_P_temp
        }
      }else{
        index_clusters_P <- index_clusters_P_temp
      }


      # non overlapping clusters:
      final_index_P <- non_overlap(index_clusters_P, matrix_clusters)

      pval_clusters_P <- pvals_P[final_index_P]
      sites_clusters_P <- lapply(final_index_P, function(j) which(matrix_clusters[,j]==1))
      if(is.null(matrix_clusters_provided)){
        centres_clusters_P <- matrix(centres[final_index_P,], ncol = 2)
        radius_clusters_P <- radius[final_index_P]
        if(is.null(sites_areas) == FALSE){
          areas_clusters_P <- areas[final_index_P]
        }
      }else{
        if(is.null(area_clusters_provided) == FALSE){
          areas_clusters_P <- areas[final_index_P]
        }
      }

    }
    if("R" %in% method){
      results_R <- results$R
      results_R <- matrix(unlist(results_R), ncol = nb_clusters, byrow = TRUE)
      stat_MC_R <- rowMaxs(results_R)
      pvals_R <- sapply(1:nb_clusters, function(j) (length(which(stat_MC_R >= index_R[j]))+1)/(MC+1))
      index_clusters_R_temp <- which(pvals_R <= typeI)

      if(length(index_clusters_R_temp) > 0){
        ordre_R <- order(index_R[index_clusters_R_temp], decreasing = TRUE)
        index_clusters_R_temp <- index_clusters_R_temp[ordre_R]

        # a posteriori filtering
        if(filtering_post == TRUE){
          if(type_minimaxi_post == "sites/indiv"){
            index_clusters_R <- post_filt_nb_sites(mini_post, maxi_post, nb_sites, index_clusters_R_temp, matrix_clusters)
          }
          if(type_minimaxi_post == "radius"){
            index_clusters_R <- post_filt_radius(mini_post, maxi_post, radius, index_clusters_R_temp)
          }
          if(type_minimaxi_post == "area"){
            if(is.null(matrix_clusters_provided)){
              areas_clusters <- areas
            }else{
              areas_clusters <- area_clusters_provided
            }

            index_clusters_R <- post_filt_area(mini_post, maxi_post, areas_clusters, index_clusters_R_temp)

          }
        }else{
          index_clusters_R <- index_clusters_R_temp
        }
      }else{
        index_clusters_R <- index_clusters_R_temp
      }


      # non overlapping clusters:
      final_index_R <- non_overlap(index_clusters_R, matrix_clusters)

      pval_clusters_R <- pvals_R[final_index_R]
      sites_clusters_R <- lapply(final_index_R, function(j) which(matrix_clusters[,j]==1))
      if(is.null(matrix_clusters_provided)){
        centres_clusters_R <- matrix(centres[final_index_R,], ncol = 2)
        radius_clusters_R <- radius[final_index_R]
        if(is.null(sites_areas) == FALSE){
          areas_clusters_R <- areas[final_index_R]
        }
      }else{
        if(is.null(area_clusters_provided) == FALSE){
          areas_clusters_R <- areas[final_index_R]
        }
      }

    }
  }else{
      stop("The data must be a list")
  }

  output <- list()

  if(is.null(matrix_clusters_provided)){
    if("LH" %in% method){
      if(is.null(sites_areas)){
        output <- c(output, list(sites_clusters_LH = sites_clusters_LH, pval_clusters_LH = pval_clusters_LH, centres_clusters_LH = centres_clusters_LH, radius_clusters_LH = radius_clusters_LH))
      }else{
        output <- c(output, list(sites_clusters_LH = sites_clusters_LH, pval_clusters_LH = pval_clusters_LH, centres_clusters_LH = centres_clusters_LH, radius_clusters_LH = radius_clusters_LH, areas_clusters_LH = areas_clusters_LH))
      }
    }
    if("W" %in% method){
      if(is.null(sites_areas)){
        output <- c(output, list(sites_clusters_W = sites_clusters_W, pval_clusters_W = pval_clusters_W, centres_clusters_W = centres_clusters_W, radius_clusters_W = radius_clusters_W))
      }else{
        output <- c(output, list(sites_clusters_W = sites_clusters_W, pval_clusters_W = pval_clusters_W, centres_clusters_W = centres_clusters_W, radius_clusters_W = radius_clusters_W, areas_clusters_W = areas_clusters_W))
      }
    }
    if("P" %in% method){
      if(is.null(sites_areas)){
        output <- c(output, list(sites_clusters_P = sites_clusters_P, pval_clusters_P = pval_clusters_P, centres_clusters_P = centres_clusters_P, radius_clusters_P = radius_clusters_P))
      }else{
        output <- c(output, list(sites_clusters_P = sites_clusters_P, pval_clusters_P = pval_clusters_P, centres_clusters_P = centres_clusters_P, radius_clusters_P = radius_clusters_P, areas_clusters_P = areas_clusters_P))
      }
    }
    if("R" %in% method){
      if(is.null(sites_areas)){
        output <- c(output, list(sites_clusters_R = sites_clusters_R, pval_clusters_R = pval_clusters_R, centres_clusters_R = centres_clusters_R, radius_clusters_R = radius_clusters_R))
      }else{
        output <- c(output, list(sites_clusters_R = sites_clusters_R, pval_clusters_R = pval_clusters_R, centres_clusters_R = centres_clusters_R, radius_clusters_R = radius_clusters_R, areas_clusters_R = areas_clusters_R))
      }
    }
    output <- c(output, list(system = system))
  }else{
    if("LH" %in% method){
      if(is.null(area_clusters_provided)){
        output <- c(output, list(sites_clusters_LH = sites_clusters_LH, pval_clusters_LH = pval_clusters_LH))
      }else{
        output <- c(output, list(sites_clusters_LH = sites_clusters_LH, pval_clusters_LH = pval_clusters_LH, areas_clusters_LH = areas_clusters_LH))
      }
    }
    if("W" %in% method){
      if(is.null(area_clusters_provided)){
        output <- c(output, list(sites_clusters_W = sites_clusters_W, pval_clusters_W = pval_clusters_W))
      }else{
        output <- c(output, list(sites_clusters_W = sites_clusters_W, pval_clusters_W = pval_clusters_W, areas_clusters_W = areas_clusters_W))
      }
    }
    if("P" %in% method){
      if(is.null(area_clusters_provided)){
        output <- c(output, list(sites_clusters_P = sites_clusters_P, pval_clusters_P = pval_clusters_P))
      }else{
        output <- c(output, list(sites_clusters_P = sites_clusters_P, pval_clusters_P = pval_clusters_P, areas_clusters_P = areas_clusters_P))
      }
    }
    if("R" %in% method){
      if(is.null(area_clusters_provided)){
        output <- c(output, list(sites_clusters_R = sites_clusters_R, pval_clusters_R = pval_clusters_R))
      }else{
        output <- c(output, list(sites_clusters_R = sites_clusters_R, pval_clusters_R = pval_clusters_R, areas_clusters_R = areas_clusters_R))
      }
    }
  }


  return(output)
}



################################################################
##' @title DFFSS scan procedure
##'
##' @description This function returns the significant clusters and their associated p-value for the DFFSS.
##'
##' @param data matrix. Matrix of the data, the rows correspond to the sites (or to the individuals if the observations are by individuals and not by sites) and each column represents an observation time. The times must be the same for each site/individual.
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
##' @param nbCPU numeric. Number of CPU. If nbCPU > 1 parallelization is done. By default: 1.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' coords <- coordinates(map_sites)
##'
##' res_dffss <- DFFSS(data = funi_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_dffss <- DFFSS(data = funi_data[indices,], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 9)
##' }
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
DFFSS <- function(data, sites_coord = NULL, system = NULL, mini = 1, maxi = nrow(sites_coord)/2, type_minimaxi = "sites/indiv", mini_post = NULL, maxi_post = NULL, type_minimaxi_post = "sites/indiv", sites_areas = NULL, MC = 999, typeI = 0.05, nbCPU = 1){

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

  if(class(data)[1]=="matrix"){
    if(nrow(matrix_clusters)!=nrow(data)){
      stop("The number of sites must be the same in the matrix of clusters and in the data")
    }
    if(ncol(data)==1){
      stop("There must be at least two observation times")
    }

    nb_sites <- nrow(matrix_clusters)

    index <- pointwise_dfree(data, matrix_clusters)

    nb_clusters <- ncol(matrix_clusters)
    nb_sites <- nrow(matrix_clusters)
    permutations <- permutate(1:nb_sites, MC)

    generation_signif <- lapply(1:MC, function(g) data[permutations[g,],])

    num_cores <- detectCores()

    if(nbCPU <=1 | num_cores <2){
      results <- lapply(1:MC, function(i) pointwise_dfree(generation_signif[[i]], matrix_clusters))
    }else{
      if(num_cores < nbCPU){
        nbCPU <- num_cores
      }

      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(nbCPU)
        results <- parLapply(cl, 1:MC, function(i) pointwise_dfree(generation_signif[[i]], matrix_clusters))
        stopCluster(cl)

      }else{
        results <- mclapply(1:MC, function(i) pointwise_dfree(generation_signif[[i]], matrix_clusters), mc.cores = nbCPU)
      }
    }
    results <- matrix(unlist(results), ncol = nb_clusters, byrow = TRUE)
    stat_MC <- rowMaxs(results)
    pvals <- sapply(1:nb_clusters, function(j) (length(which(stat_MC >= index[j]))+1)/(MC+1))

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
    stop("The data must be a matrix")
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
##' @title MDFFSS scan procedure
##'
##' @description This function returns the significant clusters and their associated p-value for the MDFFSS.
##'
##' @param data list of numeric matrices. List of nb_sites (or nb_individuals if the observations are by individuals and not by sites) matrices of the data, the rows correspond to the variables and each column represents an observation time. The times must be the same for each site/individual.
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
##' @param nbCPU numeric. Number of CPU. If nbCPU > 1 parallelization is done. By default: 1.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("fmulti_data")
##' coords <- coordinates(map_sites)
##'
##' res_mdffss <- MDFFSS(data = fmulti_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("fmulti_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_mdffss <- MDFFSS(data = fmulti_data[indices], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 9)
##' }
##'
##'
##' @return The list of the following elements:
##' \itemize{
##' \item sites_clusters: the index of the sites (or individuals) included in the significant clusters. The clusters are listed in their order of detection. The secondary clusters are defined according to Kulldorff.
##' \item pval_clusters: the associated p-values
##' \item centres_clusters: the coordinates of the centres of each cluster
##' \item radius_clusters: the radius of the clusters in km if system = WGS84 or in the user units
##' \item areas_clusters: the areas of the clusters (in the same units as sites_areas). Only if sites_areas is not NULL.
##' \item system: the system of coordinates
##' }
##' @export
##'
##'
MDFFSS <- function(data, sites_coord = NULL, system = NULL, mini = 1, maxi = nrow(sites_coord)/2, type_minimaxi = "sites/indiv", mini_post = NULL, maxi_post = NULL, type_minimaxi_post = "sites/indiv", sites_areas = NULL, MC = 999, typeI = 0.05, nbCPU = 1){

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

  if(class(data)[1]=="list"){
    if(length(data)!=nrow(matrix_clusters)){
      stop("The data must contain the same number of sites than the matrix of clusters")
    }
    nb_rows <- sapply(1:length(data), function(r) nrow(data[[r]]))
    nb_cols <- sapply(1:length(data), function(r) ncol(data[[r]]))

    if(sum((nb_rows!=nb_rows[1])) + sum((nb_cols!=nb_cols[1])) > 0){
      stop("All matrices of the list must be of the same dimensions")
    }
    if(nb_rows[1]==1){
      stop("There must be at least two variables (rows)")
    }
    if(nb_cols[1]==1){
      stop("There must be at least two observation times (columns)")
    }


    nb_sites <- nrow(matrix_clusters)
    nb_clusters <- ncol(matrix_clusters)

    index <- dfree_index_multi(data, matrix_clusters)


    permutations <- permutate(1:nb_sites, MC)
    generation_signif <- list()
    for(g in 1:MC){
      generation_signif[[g]] <- list()
      for(s in 1:nb_sites){
        generation_signif[[g]][[s]] <- data[[permutations[g,s]]]
      }
    }

    num_cores <- detectCores()

    if(nbCPU <=1 | num_cores <2){
      results <- lapply(1:MC, function(i) dfree_index_multi(generation_signif[[i]], matrix_clusters))
    }else{
      if(num_cores < nbCPU){
        nbCPU <- num_cores
      }

      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(nbCPU)
        results <- parLapply(cl, 1:MC, function(i) dfree_index_multi(generation_signif[[i]], matrix_clusters))
        stopCluster(cl)

      }else{
        results <- mclapply(1:MC, function(i) dfree_index_multi(generation_signif[[i]], matrix_clusters), mc.cores = nbCPU)
      }
    }

    results <- matrix(unlist(results), ncol = nb_clusters, byrow = TRUE)
    stat_MC <- rowMaxs(results)
    pvals <- sapply(1:nb_clusters, function(j) (length(which(stat_MC >= index[j]))+1)/(MC+1))
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
    stop("The data must be a list")
  }

  output <- list()

  if(is.null(matrix_clusters_provided)){
    if(is.null(sites_areas)){
      output <- c(output, list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, centres_clusters = centres_clusters, radius_clusters = radius_clusters))
    }else{
      output <- c(output, list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, centres_clusters = centres_clusters, radius_clusters = radius_clusters, areas_clusters = areas_clusters))
    }

    output <- c(output, list(system = system))
  }else{
      if(is.null(area_clusters_provided)){
        output <- c(output, list(sites_clusters = sites_clusters, pval_clusters = pval_clusters))
      }else{
        output <- c(output, list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, areas_clusters = areas_clusters))
      }
  }


  return(output)
}


################################################################
##' @title MRBFSS scan procedure
##'
##' @description This function returns the significant clusters and their associated p-value for the MRBFSS.
##'
##' @param data list of numeric matrices. List of nb_sites (or nb_individuals if the observations are by individuals and not by sites) matrices of the data, the rows correspond to the variables and each column represents an observation time. The times must be the same for each site/individual.
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
##' @param nbCPU numeric. Number of CPU. If nbCPU > 1 parallelization is done. By default: 1.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("fmulti_data")
##' coords <- coordinates(map_sites)
##'
##' res_mrbfss <- MRBFSS(data = fmulti_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("fmulti_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_mrbfss <- MRBFSS(data = fmulti_data[indices], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = 2, MC = 9)
##' }
##'
##'
##' @return The list of the following elements:
##' \itemize{
##' \item sites_clusters: the index of the sites (or individuals) included in the significant clusters. The clusters are listed in their order of detection. The secondary clusters are defined according to Kulldorff.
##' \item pval_clusters: the associated p-values
##' \item centres_clusters: the coordinates of the centres of each cluster
##' \item radius_clusters: the radius of the clusters in km if system = WGS84 or in the user units
##' \item areas_clusters: the areas of the clusters (in the same units as sites_areas). Only if sites_areas is not NULL.
##' \item system: the system of coordinates
##' }
##' @export
##'
##'
MRBFSS <- function(data, sites_coord = NULL, system = NULL, mini = 1, maxi = nrow(sites_coord)/2, type_minimaxi = "sites/indiv", mini_post = NULL, maxi_post = NULL, type_minimaxi_post = "sites/indiv", sites_areas = NULL, MC = 999, typeI = 0.05, nbCPU = 1){

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

  if(class(data)[1]=="list"){
    if(length(data)!=nrow(matrix_clusters)){
      stop("The data must contain the same number of sites than the matrix of clusters")
    }
    nb_rows <- sapply(1:length(data), function(r) nrow(data[[r]]))
    nb_cols <- sapply(1:length(data), function(r) ncol(data[[r]]))

    if(sum((nb_rows!=nb_rows[1])) + sum((nb_cols!=nb_cols[1])) > 0){
      stop("All matrices of the list must be of the same dimensions")
    }
    if(nb_rows[1]==1){
      stop("There must be at least two variables (rows)")
    }
    if(nb_cols[1]==1){
      stop("There must be at least two observation times (columns)")
    }


    new_data <- transform_data(data)

    nb_clusters <- ncol(matrix_clusters)
    nb_sites <- nrow(matrix_clusters)
    nb_times <- ncol(data[[1]])

    index <- pointwise_wmw_multi(new_data, matrix_clusters)


    permutations <- permutate(1:nb_sites, MC)
    generation_signif <- list()
    for(g in 1:MC){
      generation_signif[[g]] <- list()
      for(t in 1:nb_times){
        generation_signif[[g]][[t]] <- new_data[[t]][permutations[g,],]
      }
    }

    num_cores <- detectCores()

    if(nbCPU <=1 | num_cores <2){
      results <- lapply(1:MC, function(i) pointwise_wmw_multi(generation_signif[[i]], matrix_clusters))
    }else{

      if(num_cores < nbCPU){
        nbCPU <- num_cores
      }

      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(nbCPU)
        results <- parLapply(cl, 1:MC, function(i) pointwise_wmw_multi(generation_signif[[i]], matrix_clusters))
        stopCluster(cl)

      }else{
        results <- mclapply(1:MC, function(i) pointwise_wmw_multi(generation_signif[[i]], matrix_clusters), mc.cores = nbCPU)
      }
    }

    results <- matrix(unlist(results), ncol = nb_clusters, byrow = TRUE)
    stat_MC <- rowMaxs(results)
    pvals <- sapply(1:nb_clusters, function(j) (length(which(stat_MC >= index[j]))+1)/(MC+1))
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
    stop("The data must be a list")
  }

  output <- list()

  if(is.null(matrix_clusters_provided)){
    if(is.null(sites_areas)){
      output <- c(output, list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, centres_clusters = centres_clusters, radius_clusters = radius_clusters))
    }else{
      output <- c(output, list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, centres_clusters = centres_clusters, radius_clusters = radius_clusters, areas_clusters = areas_clusters))
    }

    output <- c(output, list(system = system))
  }else{
    if(is.null(area_clusters_provided)){
      output <- c(output, list(sites_clusters = sites_clusters, pval_clusters = pval_clusters))
    }else{
      output <- c(output, list(sites_clusters = sites_clusters, pval_clusters = pval_clusters, areas_clusters = areas_clusters))
    }
  }


  return(output)
}



################################################################
##' @title URBFSS scan procedure
##'
##' @description This function returns the significant clusters and their associated p-value for the URBFSS.
##'
##' @param data matrix. Matrix of the data, the rows correspond to the sites (or to the individuals if the observations are by individuals and not by sites) and each column represents an observation time. The times must be the same for each site/individual.
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
##' @param nbCPU numeric. Number of CPU. If nbCPU > 1 parallelization is done. By default: 1.
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' coords <- coordinates(map_sites)
##'
##' res_urbfss <- URBFSS(data = funi_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2) }
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_urbfss <- URBFSS(data = funi_data[indices,], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 9)
##' }
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
URBFSS <- function(data, sites_coord = NULL, system = NULL, mini = 1, maxi = nrow(sites_coord)/2, type_minimaxi = "sites/indiv", mini_post = NULL, maxi_post = NULL, type_minimaxi_post = "sites/indiv", sites_areas = NULL, MC = 999, typeI = 0.05, nbCPU = 1){

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

  if(class(data)[1]=="matrix"){
    if(nrow(matrix_clusters)!=nrow(data)){
      stop("The number of sites must be the same in the matrix of clusters and in the data")
    }
    if(ncol(data)==1){
      stop("There must be at least two observation times")
    }

    nb_sites <- nrow(matrix_clusters)

    rank_data <- t(colRanks(data, ties.method = "average"))

    index <- pointwise_wmw_uni(rank_data, matrix_clusters)

    nb_clusters <- ncol(matrix_clusters)
    nb_sites <- nrow(matrix_clusters)
    permutations <- permutate(1:nb_sites, MC)

    generation_signif <- lapply(1:MC, function(g) rank_data[permutations[g,],])

    num_cores <- detectCores()

    if(nbCPU <=1 | num_cores <2){
      results <- lapply(1:MC, function(i) pointwise_wmw_uni(generation_signif[[i]], matrix_clusters))
    }else{
      if(num_cores < nbCPU){
        nbCPU <- num_cores
      }

      if(Sys.info()["sysname"] == "Windows"){
        cl <- makeCluster(nbCPU)
        results <- parLapply(cl, 1:MC, function(i) pointwise_wmw_uni(generation_signif[[i]], matrix_clusters))
        stopCluster(cl)

      }else{
        results <- mclapply(1:MC, function(i) pointwise_wmw_uni(generation_signif[[i]], matrix_clusters), mc.cores = nbCPU)
      }
    }
    results <- matrix(unlist(results), ncol = nb_clusters, byrow = TRUE)
    stat_MC <- rowMaxs(results)
    pvals <- sapply(1:nb_clusters, function(j) (length(which(stat_MC >= index[j]))+1)/(MC+1))

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
    stop("The data must be a matrix")
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
