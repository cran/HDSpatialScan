########################################################################################################################################################
##' @title Summary of the clusters
##'
##' @description This function gives a summary of the clusters in a table
##'
##' @param output_clusters list. List of the sites in the clusters: it is the sites_clusters of the output of NPFSS, PFSS, DFFSS, URBFSS, MDFFSS, MRBFSS, MG, MNP, UG or UNP, or the sites_clusters_LH/sites_clusters_W/sites_clusters_P/sites_clusters_R of the MPFSS.
##' @param data list of numeric matrices or a matrix, or a vector. List of nb_sites (or nb_individuals if the observations are by individuals and not by site) matrices of the data, the rows correspond to the variables and each column represents an observation time (multivariate case) ; or Matrix of the data, the rows correspond to the sites (or to the individuals) and each column represents an observation time (univariate case). The times must be the same for each site/individual and they may need to be equally spaced depending on the scan method ; Or vector a the data (univariate non-functional framework) in which each element corresponds to a site or an individual.
##' @param type character. Type of the data: "multi" (for multivariate data) or "funct" (for functional univariate data or functional multivariate data). Only taken into account when the data is a matrix : if the data is a vector (univariate context) the value has no importance.
##' @param type_summ character. "param" or "nparam". "param" gives the mean and the sd for each variable in the clusters, outside, and globally and "nparam" gives the Q25, Q50 and Q75 quantiles for each variables in the clusters, outside, and globally.
##' @param nb_digits integer. Number of decimals in the means, sds and quantiles.
##' @param html logical. If TRUE the table is an HTML table. This allows a better readability but a less good extraction of the information than if html = FALSE: in the latter case the table is simply returned in the console.
##' @param variable_names character vector. The names of the variables if you want the names to appear. By default NULL: the names are var1, var2, etc. Only considered for multivariate functional or multivariate non-functional data.
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
##' plot_summary(output_clusters = res_npfss$sites_clusters, data = funi_data,
##' type = "funct", type_summ = "param")
##' plot_summary(output_clusters = res_npfss$sites_clusters, data = funi_data,
##' type = "funct", type_summ = "nparam")}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_npfss <- NPFSS(data = funi_data[indices,], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 99)
##' if(length(res_npfss$sites_clusters)>0){
##' plot_summary(output_clusters = res_npfss$sites_clusters,
##' data = funi_data[indices,],
##' type = "funct", type_summ = "param", html = FALSE)
##' }
##'
##' }
##'
##'
##' @return No value returned, displays the results in the console (html = FALSE) or in the Viewer (html = TRUE)
##'
##' @export
##'
plot_summary <- function(output_clusters, data, type = "multi", type_summ = "param", nb_digits = 3, html = TRUE, variable_names = NULL){
  if(is.logical(html)==FALSE){
    stop("The argument html must be logical")
  }
  if(class(data)[1]=="matrix"){
    if(!(type%in% c("multi", "funct"))){
      stop("The type must be multi or funct")
    }
    if(type == "multi"){
      final_type <- "multi_nfunct"
    }else{
      final_type <- "uni_funct"
    }
  }else{
    if(class(data)[1]=="list"){
      final_type <- "multi_funct"
    }else{
      if(is.vector(data)){
        final_type <- "uni_nfunct"
      }else{
        stop("The data must be a vector (univariate data), a matrix (multivariate or functional univariate) or a list (functional multivariate)")
      }
    }
  }

  if(final_type == "uni_funct"){

    if(!(type_summ%in% c("param", "nparam"))){
      stop("type_summ must be param or nparam")
    }
    if(type_summ == "param"){
      if(!(class(nb_digits) %in% c("integer", "numeric"))){
        stop("nb_digits must be an integer")
      }else{
        if(as.integer(nb_digits) != nb_digits){
          stop("nb_digits must be an integer")
        }
      }
      data_summary <- data.frame(matrix(nrow = 2*length(output_clusters)+1, ncol = 3))
      colnames(data_summary) <- c("Number of sites", "Mean", "Sd")
      data_summary[1,] <- c(nrow(data), round(mean(data), digits = nb_digits), round(sd(data), digits = nb_digits))

      row_names <- c("Overall")

      for(i in 1:length(output_clusters)){
        data_summary[2*i,] <- c(length(output_clusters[[i]]), round(mean(data[output_clusters[[i]],]), digits = nb_digits), round(sd(data[output_clusters[[i]],]), digits = nb_digits))
        data_summary[2*i+1,] <- c(nrow(data) - length(output_clusters[[i]]), round(mean(data[-output_clusters[[i]],]), digits = nb_digits), round(sd(data[-output_clusters[[i]],]), digits = nb_digits))
        row_names <- c(row_names, paste("Inside cluster ", i, sep = ""), paste("Outside cluster ", i, sep = ""))
      }
    }else{
      data_summary <- data.frame(matrix(nrow = 2*length(output_clusters)+1, ncol = 4))

      colnames(data_summary) <- c("Number of sites", "Q25", "Median", "Q75")

      data_summary[1,] <- c(nrow(data), round(as.numeric(quantile(data, probs = c(0.25,0.5,0.75))), digits = nb_digits))

      row_names <- c("Overall")

      for(i in 1:length(output_clusters)){
        data_summary[2*i,] <- c(length(output_clusters[[i]]), round(as.numeric(quantile(data[output_clusters[[i]],], probs = c(0.25,0.5,0.75))), digits = nb_digits))
        data_summary[2*i+1,] <- c(nrow(data) - length(output_clusters[[i]]), round(as.numeric(quantile(data[-output_clusters[[i]],], probs = c(0.25,0.5,0.75))), digits = nb_digits))
        row_names <- c(row_names, paste("Inside cluster ", i, sep = ""), paste("Outside cluster ", i, sep = ""))
      }
    }
  }

  if(final_type == "multi_funct"){

    if(is.null(variable_names)){
      variable_names <- paste("var", c(1:nrow(data[[1]])), sep = "")
    }else{
      if(length(variable_names) != nrow(data[[1]])){
        stop("There must be the same number of elements in variable_names than the number of variables in data")
      }
    }

    if(!(type_summ%in% c("param", "nparam"))){
      stop("type_summ must be param or nparam")
    }
    if(type_summ == "param"){
      data_summary <- data.frame(matrix(nrow = 2*length(output_clusters)+1, ncol = 1 + 2*nrow(data[[1]])))
      columns <- c("Number of sites")
      for(v in 1:nrow(data[[1]])){
        columns <- c(columns, paste("Mean ", variable_names[v], sep = ""))
        columns <- c(columns, paste("Sd ", variable_names[v], sep = ""))
      }

      colnames(data_summary) <- columns
      data_summary[1,1] <- c(length(data))
      row_names <- c("Overall")
      for(i in 1:length(output_clusters)){
        data_summary[2*i,1] <- c(length(output_clusters[[i]]))
        data_summary[2*i+1,1] <- c(length(data) - length(output_clusters[[i]]))
        row_names <- c(row_names, paste("Inside cluster ", i, sep = ""), paste("Outside cluster ", i, sep = ""))

      }
      for(v in 1:nrow(data[[1]])){
        temp <- matrix(ncol = ncol(data[[1]]), nrow = length(data))
        for(l in 1:length(data)){
          temp[l,] <- data[[l]][v,]
        }
        data_summary[1,c(2*v, 1+2*v)] <- c(round(mean(matrix(temp, ncol = ncol(data[[1]]))), digits = nb_digits), round(sd(matrix(temp, ncol = ncol(data[[1]]))), digits = nb_digits))
        for(i in 1:length(output_clusters)){
          data_summary[2*i,c(2*v, 1+2*v)] <- c(round(mean(matrix(temp[output_clusters[[i]],],ncol = ncol(data[[1]]))), digits = nb_digits), round(sd(matrix(temp[output_clusters[[i]],], ncol = ncol(data[[1]]))), digits = nb_digits))
          data_summary[2*i+1,c(2*v, 1+2*v)] <- c(round(mean(matrix(temp[-output_clusters[[i]],],ncol = ncol(data[[1]]))), digits = nb_digits), round(sd(matrix(temp[-output_clusters[[i]],], ncol = ncol(data[[1]]))), digits = nb_digits))
        }
      }

    }else{

      data_summary <- data.frame(matrix(nrow = 2*length(output_clusters)+1, ncol = 1 + 3*nrow(data[[1]])))
      columns <- c("Number of sites")
      for(v in 1:nrow(data[[1]])){
        columns <- c(columns, paste("Q25 ", variable_names[v], sep = ""))
        columns <- c(columns, paste("Median ", variable_names[v], sep = ""))
        columns <- c(columns, paste("Q75 ", variable_names[v], sep = ""))

      }

      colnames(data_summary) <- columns
      data_summary[1,1] <- c(length(data))
      row_names <- c("Overall")
      for(i in 1:length(output_clusters)){
        data_summary[2*i,1] <- c(length(output_clusters[[i]]))
        data_summary[2*i+1,1] <- c(length(data) - length(output_clusters[[i]]))
        row_names <- c(row_names, paste("Inside cluster ", i, sep = ""), paste("Outside cluster ", i, sep = ""))

      }
      for(v in 1:nrow(data[[1]])){
        temp <- matrix(ncol = ncol(data[[1]]), nrow = length(data))
        for(l in 1:length(data)){
          temp[l,] <- data[[l]][v,]
        }
        data_summary[1,c(3*v-1, 3*v, 3*v+1)] <- round(as.numeric(c(quantile(matrix(temp, ncol = ncol(data[[1]])), probs = c(0.25, 0.5, 0.75)))), digits = nb_digits)
        for(i in 1:length(output_clusters)){
          data_summary[2*i,c(3*v-1, 3*v, 3*v+1)] <- round(as.numeric(c(quantile(matrix(temp[output_clusters[[i]],], ncol = ncol(data[[1]])), probs = c(0.25, 0.5, 0.75)))), digits = nb_digits)
          data_summary[2*i+1,c(3*v-1, 3*v, 3*v+1)] <- round(as.numeric(c(quantile(matrix(temp[-output_clusters[[i]],], ncol = ncol(data[[1]])), probs = c(0.25, 0.5, 0.75)))), digits = nb_digits)
        }
      }
    }
  }

  if(final_type == "multi_nfunct"){

    if(is.null(variable_names)){
      variable_names <- paste("var", c(1:ncol(data)), sep = "")
    }else{
      if(length(variable_names) != ncol(data)){
        stop("There must be the same number of elements in variable_names than the number of columns in data")
      }
    }

    if(!(type_summ%in% c("param", "nparam"))){
      stop("type_summ must be param or nparam")
    }
    if(type_summ == "param"){
      data_summary <- data.frame(matrix(nrow = 2*length(output_clusters)+1, ncol = 1+2*ncol(data)))
      columns <- c("Number of sites")
      for(v in 1:ncol(data)){
        columns <- c(columns, paste("Mean ", variable_names[v], sep = ""))
        columns <- c(columns, paste("Sd ", variable_names[v], sep = ""))
      }

      colnames(data_summary) <- columns

      temp_mean <- round(colMeans(matrix(data, ncol = ncol(data))), digits = nb_digits)
      temp_sd <- round(colSds(matrix(data, ncol=ncol(data))), digits = nb_digits)

      complete_mean_sd <- numeric(2*length(temp_mean))
      for(v in 1:length(temp_mean)){
        complete_mean_sd[c(2*(v-1)+1,2*(v-1)+2)] <- c(temp_mean[v], temp_sd[v])
      }


      data_summary[1,] <- c(nrow(data), complete_mean_sd)
      row_names <- c("Overall")
      for(i in 1:length(output_clusters)){

        temp_mean <- round(colMeans(matrix(data[output_clusters[[i]],], ncol = ncol(data))), digits = nb_digits)
        temp_sd <- round(colSds(matrix(data[output_clusters[[i]],], ncol = ncol(data))), digits = nb_digits)
        complete_mean_sd <- numeric(2*length(temp_mean))
        for(v in 1:length(temp_mean)){
          complete_mean_sd[c(2*(v-1)+1,2*(v-1)+2)] <- c(temp_mean[v], temp_sd[v])
        }
        data_summary[2*i,] <- c(length(output_clusters[[i]]), complete_mean_sd)

        temp_mean <- round(colMeans(matrix(data[-output_clusters[[i]],], ncol = ncol(data))), digits = nb_digits)
        temp_sd <- round(colSds(matrix(data[-output_clusters[[i]],], ncol = ncol(data))), digits = nb_digits)
        complete_mean_sd <- numeric(2*length(temp_mean))
        for(v in 1:length(temp_mean)){
          complete_mean_sd[c(2*(v-1)+1,2*(v-1)+2)] <- c(temp_mean[v], temp_sd[v])
        }

        data_summary[2*i+1,] <- c(nrow(data) - length(output_clusters[[i]]), complete_mean_sd)

        row_names <- c(row_names, paste("Inside cluster ", i, sep = ""), paste("Outside cluster ", i, sep = ""))

      }
    }else{
      data_summary <- data.frame(matrix(nrow = 2*length(output_clusters)+1, ncol = 1+3*ncol(data)))
      columns <- c("Number of sites")
      for(v in 1:ncol(data)){
        columns <- c(columns, paste("Q25 ", variable_names[v], sep = ""))
        columns <- c(columns, paste("Median ", variable_names[v], sep = ""))
        columns <- c(columns, paste("Q75 ", variable_names[v], sep = ""))
      }

      colnames(data_summary) <- columns

      temp_quantiles <- round((colQuantiles(matrix(data, ncol = ncol(data)), probs = c(0.25,0.5,0.75))), digits = nb_digits)

      complete_quantiles <- c()
      for(v in 1:nrow(temp_quantiles)){
        complete_quantiles <- c(complete_quantiles, temp_quantiles[v,])
      }

      data_summary[1,] <- c(nrow(data), complete_quantiles)

      row_names <- c("Overall")

      for(i in 1:length(output_clusters)){

        temp_quantiles <- round((colQuantiles(matrix(data[output_clusters[[i]],], ncol = ncol(data)), probs = c(0.25,0.5,0.75))), digits = nb_digits)

        complete_quantiles <- c()
        for(v in 1:nrow(temp_quantiles)){
          complete_quantiles <- c(complete_quantiles, temp_quantiles[v,])
        }
        data_summary[2*i,] <- c(length(output_clusters[[i]]), complete_quantiles)

        temp_quantiles <- round((colQuantiles(matrix(data[-output_clusters[[i]],], ncol = ncol(data)), probs = c(0.25,0.5,0.75))), digits = nb_digits)

        complete_quantiles <- c()
        for(v in 1:nrow(temp_quantiles)){
          complete_quantiles <- c(complete_quantiles, temp_quantiles[v,])
        }

        data_summary[2*i+1,] <- c(nrow(data) - length(output_clusters[[i]]), complete_quantiles)

        row_names <- c(row_names, paste("Inside cluster ", i, sep = ""), paste("Outside cluster ", i, sep = ""))

      }
    }

  }

  if(final_type == "uni_nfunct"){
    if(!(type_summ%in% c("param", "nparam"))){
      stop("type_summ must be param or nparam")
    }
    if(type_summ == "param"){
      data_summary <- data.frame(matrix(nrow = 2*length(output_clusters)+1, ncol = 1+2))
      columns <- c("Number of sites", "Mean", "Sd")
      colnames(data_summary) <- columns

      temp_mean <- round(mean(data), digits = nb_digits)
      temp_sd <- round(sd(data), digits = nb_digits)

      complete_mean_sd <- c(temp_mean, temp_sd)

      data_summary[1,] <- c(length(data), complete_mean_sd)
      row_names <- c("Overall")
      for(i in 1:length(output_clusters)){

        temp_mean <- round(mean(data[output_clusters[[i]]]), digits = nb_digits)
        temp_sd <- round(sd(data[output_clusters[[i]]]), digits = nb_digits)
        complete_mean_sd <- c(temp_mean, temp_sd)

        data_summary[2*i,] <- c(length(output_clusters[[i]]), complete_mean_sd)

        temp_mean <- round(mean(data[-output_clusters[[i]]]), digits = nb_digits)
        temp_sd <- round(sd(data[-output_clusters[[i]]]), digits = nb_digits)
        complete_mean_sd <- c(temp_mean, temp_sd)

        data_summary[2*i+1,] <- c(length(data) - length(output_clusters[[i]]), complete_mean_sd)

        row_names <- c(row_names, paste("Inside cluster ", i, sep = ""), paste("Outside cluster ", i, sep = ""))

      }
    }else{
      data_summary <- data.frame(matrix(nrow = 2*length(output_clusters)+1, ncol = 1+3))
      columns <- c("Number of sites", "Q25", "Median", "Q75")

      colnames(data_summary) <- columns

      complete_quantiles <- round((quantile(data, probs = c(0.25,0.5,0.75))), digits = nb_digits)

      data_summary[1,] <- c(length(data), complete_quantiles)

      row_names <- c("Overall")

      for(i in 1:length(output_clusters)){

        complete_quantiles <- round(quantile(data[output_clusters[[i]]], probs = c(0.25,0.5,0.75)), digits = nb_digits)

        data_summary[2*i,] <- c(length(output_clusters[[i]]), complete_quantiles)

        complete_quantiles <- round(quantile(data[-output_clusters[[i]]], probs = c(0.25,0.5,0.75)), digits = nb_digits)

        data_summary[2*i+1,] <- c(length(data) - length(output_clusters[[i]]), complete_quantiles)

        row_names <- c(row_names, paste("Inside cluster ", i, sep = ""), paste("Outside cluster ", i, sep = ""))

      }
    }

  }

  rownames(data_summary) <- row_names

  if(html == TRUE){
    DT::datatable(t(data_summary))
  }else{
    output <-
    return(data.frame(t(data_summary), check.names = FALSE))
  }



}



########################################################################################################################################################
##' @title Plots the curves in the clusters (functional data)
##'
##' @description This function plots the curves in the clusters. Only for functional data.
##'
##' @param output_clusters list. List of the sites in the clusters: it is the sites_clusters of the output of NPFSS, PFSS, DFFSS, URBFSS, MDFFSS, MRBFSS, MG, MNP, UG or UNP, or the sites_clusters_LH/sites_clusters_W/sites_clusters_P/sites_clusters_R of the MPFSS.
##' @param data list of numeric matrices or a matrix. List of nb_sites (or nb_individuals if the observations are by individuals and not by site) matrices of the data, the rows correspond to the variables and each column represents an observation time (multivariate case) ; or Matrix of the data, the rows correspond to the sites (or to the individuals) and each column represents an observation time (univariate case). The times must be the same for each site/individual and they may need to be equally spaced depending on the scan method.
##' @param times numeric vector. The times of observations in the data. By default NULL: it supposes equally spaced times of [0;1].
##' @param add_mean boolean. If TRUE it adds the global mean curve in black.
##' @param add_median boolean. If TRUE it adds the global median curve in blue.
##' @param colors character. The colors to plot the clusters' curves. If length(colors)==1 then all the clusters will be plotted in this color. Else there must be the same number of elements in colors than the number of clusters
##' @param variable_names character vector. The names of the variables if you want the names to appear. By default NULL: the names are var1, var2, etc. Only considered for multivariate functional data.
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
##' plot_curves_clusters(output_clusters = res_npfss$sites_clusters, data = funi_data,
##' add_mean = TRUE, add_median = TRUE)}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_npfss <- NPFSS(data = funi_data[indices,], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 99)
##' if(length(res_npfss$sites_clusters)>0){
##' plot_curves_clusters(output_clusters = res_npfss$sites_clusters,
##' data = funi_data[indices,], add_mean = TRUE, add_median = TRUE)
##' }
##'
##' }
##'
##' @return No value returned, plots the curves.
##'
##' @export
##'
plot_curves_clusters <- function(output_clusters, data, times = NULL, add_mean = FALSE, add_median = FALSE, colors = "red", variable_names = NULL){

  if(length(colors) == 1){
    if(colors == "black" | colors == "grey" | colors == "blue"){
      stop("The color provided cannot be black, grey or blue")
    }
  }else{
    if(length(colors) != length(output_clusters)){
      stop("colors must be of length 1 or with as many colors as the length of output_clusters")
    }else{
      if(sum(colors %in% c("black", "blue", "grey")) > 0){
        stop("The colors provided cannot be black, grey or blue")
      }
    }
  }
  if(class(data)[1]=="list"){
    if(is.null(variable_names)){
      variable_names <- paste("var", c(1:nrow(data[[1]])), sep = "")
    }else{
      if(length(variable_names) != nrow(data[[1]])){
        stop("There must be the same number of elements in variable_names than the number of variables in data")
      }
    }

    if(is.null(times)){
      times <- seq(from = 0, to = 1, length.out = ncol(data[[1]]))
    }else{
      if(!(class(times) %in% c("integer", "numeric"))){
        stop("The times must be numeric")
      }else{
        if(length(times) != ncol(data[[1]])){
          stop("There must be the same number of times in times and in the data")
        }
      }
    }

    for(i in 1:length(output_clusters)){
      for(v in 1:nrow(data[[1]])){
        temp <- matrix(ncol = ncol(data[[1]]), nrow = length(data))
        for(l in 1:length(data)){
          temp[l,] <- data[[l]][v,]
        }

        plot(NULL, ylim = range(temp), xlim = range(times), xlab = "Time", ylab = variable_names[v], main = paste("Cluster ", i, sep = ""), cex.lab = 0.8, cex.main = 0.8, cex.axis = 0.8)
        for(elem in 1:nrow(temp)){
          lines(x = times, y = temp[elem,], col = "grey", lwd = 1.2)
        }

        for(elem in output_clusters[[i]]){
          if(length(colors) == 1){
            lines(x = times, y = temp[elem,], col = colors, lwd = 1.2)
          }else{
            lines(x = times, y = temp[elem,], col = colors[i], lwd = 1.2)
          }
        }
        if(length(colors) == 1){
          col_legend <- c(colors, "grey")
        }else{
          col_legend <- c(colors[i], "grey")
        }
        text_legend <- c("Sites included in the cluster", "Sites not included in the cluster")
        lty_legend <- c(1,1)
        lwd_legend <- c(1.2,1.2)

        if(add_mean == TRUE){
          mean_curve <- colMeans(matrix(temp, ncol = ncol(temp)))
          lines(x = times, y = mean_curve, col = "black", lwd = 1.5)
          col_legend <- c(col_legend, "black")
          text_legend <- c(text_legend, "Global mean")
          lty_legend <- c(lty_legend,1)
          lwd_legend <- c(lwd_legend, 1.5)
        }
        if(add_median == TRUE){
          median_curve <- colMedians(matrix(temp, ncol= ncol(temp)))
          lines(x = times, y = median_curve, col = "blue", lwd = 1.5)
          col_legend <- c(col_legend, "blue")
          text_legend <- c(text_legend, "Global median")
          lty_legend <- c(lty_legend,1)
          lwd_legend <- c(lwd_legend, 1.5)
        }
        legend("topleft", legend = text_legend, col = col_legend, lty = lty_legend, lwd = lwd_legend, bty = "n", cex = 0.8)
      }
    }

  }else{
    if(class(data)[1] == "matrix"){

      if(is.null(times)){
        times <- seq(from = 0, to = 1, length.out = ncol(data))
      }else{
        if(!(class(times) %in% c("integer", "numeric"))){
          stop("The times must be numeric")
        }else{
          if(length(times) != ncol(data)){
            stop("There must be the same number of times in times and in the data")
          }
        }
      }

      for(i in 1:length(output_clusters)){
        plot(NULL, ylim = range(data), xlim = range(times), xlab = "Time", ylab = "", main = paste("Cluster ", i, sep = ""), cex.lab = 0.8, cex.main = 0.8, cex.axis = 0.8)
        for(elem in 1:nrow(data)){
          lines(x = times, y = data[elem,], col = "grey", lwd = 1.2)
        }
        for(elem in output_clusters[[i]]){
          if(length(colors)==1){
            lines(x = times, y = data[elem,], col = colors, lwd = 1.2)
          }else{
            lines(x = times, y = data[elem,], col = colors[i], lwd = 1.2)
          }
        }

        if(length(colors) == 1){
          col_legend <- c(colors, "grey")
        }else{
          col_legend <- c(colors[i], "grey")
        }
        text_legend <- c("Sites included in the cluster", "Sites not included in the cluster")
        lty_legend <- c(1,1)
        lwd_legend <- c(1.2,1.2)


        if(add_mean == TRUE){
          mean_curve <- colMeans(matrix(data, ncol = ncol(data)))
          lines(x = times, y = mean_curve, col = "black", lwd = 1.5)
          col_legend <- c(col_legend, "black")
          text_legend <- c(text_legend, "Global mean")
          lty_legend <- c(lty_legend,1)
          lwd_legend <- c(lwd_legend, 1.5)
        }
        if(add_median == TRUE){
          median_curve <- colMedians(matrix(data, ncol = ncol(data)))
          lines(x = times, y = median_curve, col = "blue", lwd = 1.5)
          col_legend <- c(col_legend, "blue")
          text_legend <- c(text_legend, "Global median")
          lty_legend <- c(lty_legend,1)
          lwd_legend <- c(lwd_legend, 1.5)
        }
        legend("topleft", legend = text_legend, col = col_legend, lty = lty_legend, lwd = lwd_legend, bty = "n", cex = 0.8)

      }

    }else{
      stop("The data must be a list of matrices (functional multivariate) or a matrix (functional univariate)")
    }
  }
}


########################################################################################################################################################
##' @title Plots the mean or median curves in the clusters (functional data)
##'
##' @description This function plots the mean or median curves in the clusters. Only for functional data.
##'
##' @param output_clusters list. List of the sites in the clusters: it is the sites_clusters of the output of NPFSS, PFSS, DFFSS, URBFSS, MDFFSS, MRBFSS, MG, MNP, UG or UNP, or the sites_clusters_LH/sites_clusters_W/sites_clusters_P/sites_clusters_R of the MPFSS.
##' @param data list of numeric matrices or a matrix. List of nb_sites (or nb_individuals if the observations are by individuals and not by site) matrices of the data, the rows correspond to the variables and each column represents an observation time (multivariate case) ; or Matrix of the data, the rows correspond to the sites (or to the individuals) and each column represents an observation time (univariate case). The times must be the same for each site/individual and they may need to be equally spaced depending on the scan method.
##' @param times numeric vector. The times of observations in the data. By default NULL: it supposes equally spaced times of [0;1].
##' @param type character. "mean" or "median". If "mean": the mean curves in the clusters are plotted in solid lines, outside the cluster in dots, the global mean curve is in black. If "median": the median curves in the clusters are plotted in solid lines, outside the cluster in dots, the global median curve is in black.
##' @param colors character. The colors to plot the clusters' summary curves. If length(colors)==1 then all the clusters will be plotted in this color. Else there must be the same number of elements in colors than the number of clusters
##' @param variable_names character vector. The names of the variables if you want the names to appear. By default NULL: the names are var1, var2, etc. Only considered for multivariate functional data.
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
##' plot_summary_curves(output_clusters = res_npfss$sites_clusters,
##' data = funi_data, type = "mean")
##' plot_summary_curves(output_clusters = res_npfss$sites_clusters,
##' data = funi_data, type = "median")}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("funi_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_npfss <- NPFSS(data = funi_data[indices,], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 99)
##' if(length(res_npfss$sites_clusters)>0){
##' plot_summary_curves(output_clusters = res_npfss$sites_clusters,
##' data = funi_data[indices,], type = "mean")
##' }
##'
##' }
##'
##' @return No value returned, plots the curves.
##'
##' @export
##'
plot_summary_curves <- function(output_clusters, data, times = NULL, type = "mean", colors = "red", variable_names = NULL){

  if(length(colors) == 1){
    if(colors == "black" | colors == "grey"){
      stop("The color provided cannot be black or grey")
    }
  }else{
    if(length(colors) != length(output_clusters)){
      stop("colors must be of length 1 or with as many colors as the length of output_clusters")
    }else{
      if(sum(colors %in% c("black", "grey")) > 0){
        stop("The colors provided cannot be black or grey")
      }
    }
  }

  if(length(type)!=1){
    stop("type must be of length 1")
  }else{
    if(!(type %in% c("mean", "median"))){
      stop("type must be mean or median")
    }
  }
  if(class(data)[1]=="list"){

    if(is.null(variable_names)){
      variable_names <- paste("var", c(1:nrow(data[[1]])), sep = "")
    }else{
      if(length(variable_names) != nrow(data[[1]])){
        stop("There must be the same number of elements in variable_names than the number of variables in data")
      }
    }

    if(is.null(times)){
      times <- seq(from = 0, to = 1, length.out = ncol(data[[1]]))
    }else{
      if(!(class(times) %in% c("integer", "numeric"))){
        stop("The times must be numeric")
      }else{
        if(length(times) != ncol(data[[1]])){
          stop("There must be the same number of times in times and in the data")
        }
      }
    }

    for(i in 1:length(output_clusters)){
      for(v in 1:nrow(data[[1]])){
        temp <- matrix(ncol = ncol(data[[1]]), nrow = length(data))
        for(l in 1:length(data)){
          temp[l,] <- data[[l]][v,]
        }
        plot(NULL, ylim = range(temp), xlim = range(times), xlab = "Time", ylab = variable_names[v], main = paste("Cluster ", i, sep = ""), cex.lab = 0.8, cex.main = 0.8, cex.axis = 0.8)
        for(elem in 1:nrow(temp)){
          lines(x = times, y = temp[elem,], col = "grey", lwd = 1.2)
        }

        col_legend <- c("grey")
        text_legend <- c("Sites")
        lty_legend <- c(1)
        lwd_legend <- c(1.2)


        if(type == "mean"){
          mean_curve <- colMeans(matrix(temp, ncol = ncol(temp)))
          lines(x = times, y = mean_curve, col = "black", lwd = 1.5)
          if(length(colors)==1){
            lines(x = times, y = colMeans(matrix(temp[output_clusters[[i]],], ncol = ncol(temp))), col = colors, lwd = 1.2)
            lines(x = times, y = colMeans(matrix(temp[-output_clusters[[i]],], ncol = ncol(temp))), col = colors, lwd = 1.2, lty = 2)

            col_legend <- c(col_legend, "black", colors, colors)
            text_legend <- c(text_legend, "Global mean", "Mean inside the cluster", "Mean outside the cluster")
            lty_legend <- c(lty_legend, 1, 1, 2)
            lwd_legend <- c(lwd_legend, 1.5, 1.2, 1.2)

          }else{
            lines(x = times, y = colMeans(matrix(temp[output_clusters[[i]],], ncol = ncol(temp))), col = colors[i], lwd = 1.2)
            lines(x = times, y = colMeans(matrix(temp[-output_clusters[[i]],], ncol = ncol(temp))), col = colors[i], lwd = 1.2, lty = 2)
            col_legend <- c(col_legend, "black", colors[i], colors[i])
            text_legend <- c(text_legend, "Global mean", "Mean inside the cluster", "Mean outside the cluster")
            lty_legend <- c(lty_legend, 1, 1, 2)
            lwd_legend <- c(lwd_legend, 1.5, 1.2, 1.2)
          }
        }else{
          median_curve <- colMedians(matrix(temp, ncol = ncol(temp)))
          lines(x = times, y = median_curve, col = "black", lwd = 1.5)
          if(length(colors)==1){
            lines(x = times, y = colMedians(matrix(temp[output_clusters[[i]],], ncol = ncol(temp))), col = colors, lwd = 1.2)
            lines(x = times, y = colMedians(matrix(temp[-output_clusters[[i]],], ncol = ncol(temp))), col = colors, lwd = 1.2, lty = 2)
            col_legend <- c(col_legend, "black", colors, colors)
            text_legend <- c(text_legend, "Global median", "Median inside the cluster", "Median outside the cluster")
            lty_legend <- c(lty_legend, 1, 1, 2)
            lwd_legend <- c(lwd_legend, 1.5, 1.2, 1.2)
          }else{
            lines(x = times, y = colMedians(matrix(temp[output_clusters[[i]],], ncol = ncol(temp))), col = colors[i], lwd = 1.2)
            lines(x = times, y = colMedians(matrix(temp[-output_clusters[[i]],], ncol = ncol(temp))), col = colors[i], lwd = 1.2, lty = 2)
            col_legend <- c(col_legend, "black", colors[i], colors[i])
            text_legend <- c(text_legend, "Global median", "Median inside the cluster", "Median outside the cluster")
            lty_legend <- c(lty_legend, 1, 1, 2)
            lwd_legend <- c(lwd_legend, 1.5, 1.2, 1.2)
          }
        }
        legend("topleft", legend = text_legend, col = col_legend, lty = lty_legend, lwd = lwd_legend, bty = "n", cex = 0.8)
      }
    }

  }else{
    if(class(data)[1] == "matrix"){

      if(is.null(times)){
        times <- seq(from = 0, to = 1, length.out = ncol(data))
      }else{
        if(!(class(times) %in% c("integer", "numeric"))){
          stop("The times must be numeric")
        }else{
          if(length(times) != ncol(data)){
            stop("There must be the same number of times in times and in the data")
          }
        }
      }

      for(i in 1:length(output_clusters)){
        plot(NULL, ylim = range(data), xlim = range(times), xlab = "Time", ylab = "", main = paste("Cluster ", i, sep = ""), cex.lab = 0.8, cex.main = 0.8, cex.axis = 0.8)
        for(elem in 1:nrow(data)){
          lines(x = times, y = data[elem,], col = "grey", lwd = 1.2)
        }

        col_legend <- c("grey")
        text_legend <- c("Sites")
        lty_legend <- c(1)
        lwd_legend <- c(1.2)

        if(type == "mean"){
          mean_curve <- colMeans(matrix(data, ncol = ncol(data)))
          lines(x = times, y = mean_curve, col = "black", lwd = 1.5)
          if(length(colors)==1){
            lines(x = times, y = colMeans(matrix(data[output_clusters[[i]],], ncol = ncol(data))), col = colors, lwd = 1.2)
            lines(x = times, y = colMeans(matrix(data[-output_clusters[[i]],], ncol = ncol(data))), col = colors, lwd = 1.2, lty = 2)
            col_legend <- c(col_legend, "black", colors, colors)
            text_legend <- c(text_legend, "Global mean", "Mean inside the cluster", "Mean outside the cluster")
            lty_legend <- c(lty_legend, 1, 1, 2)
            lwd_legend <- c(lwd_legend, 1.5, 1.2, 1.2)
          }else{
            lines(x = times, y = colMeans(matrix(data[output_clusters[[i]],], ncol = ncol(data))), col = colors[i], lwd = 1.2)
            lines(x = times, y = colMeans(matrix(data[-output_clusters[[i]],], ncol = ncol(data))), col = colors[i], lwd = 1.2, lty = 2)
            col_legend <- c(col_legend, "black", colors[i], colors[i])
            text_legend <- c(text_legend, "Global mean", "Mean inside the cluster", "Mean outside the cluster")
            lty_legend <- c(lty_legend, 1, 1, 2)
            lwd_legend <- c(lwd_legend, 1.5, 1.2, 1.2)
          }
        }else{
          median_curve <- colMedians(matrix(data, ncol = ncol(data)))
          lines(x = times, y = median_curve, col = "black", lwd = 1.5)
          if(length(colors)==1){
            lines(x = times, y = colMedians(matrix(data[output_clusters[[i]],], ncol = ncol(data))), col = colors, lwd = 1.2)
            lines(x = times, y = colMedians(matrix(data[-output_clusters[[i]],], ncol = ncol(data))), col = colors, lwd = 1.2, lty = 2)
            col_legend <- c(col_legend, "black", colors, colors)
            text_legend <- c(text_legend, "Global median", "Median inside the cluster", "Median outside the cluster")
            lty_legend <- c(lty_legend, 1, 1, 2)
            lwd_legend <- c(lwd_legend, 1.5, 1.2, 1.2)
          }else{
            lines(x = times, y = colMedians(matrix(data[output_clusters[[i]],], ncol = ncol(data))), col = colors[i], lwd = 1.2)
            lines(x = times, y = colMedians(matrix(data[-output_clusters[[i]],], ncol = ncol(data))), col = colors[i], lwd = 1.2, lty = 2)
            col_legend <- c(col_legend, "black", colors[i], colors[i])
            text_legend <- c(text_legend, "Global median", "Median inside the cluster", "Median outside the cluster")
            lty_legend <- c(lty_legend, 1, 1, 2)
            lwd_legend <- c(lwd_legend, 1.5, 1.2, 1.2)
          }
        }
        legend("topleft", legend = text_legend, col = col_legend, lty = lty_legend, lwd = lwd_legend, bty = "n", cex = 0.8)
      }
    }else{
      stop("The data must be a list of matrices (functional multivariate) or a matrix (functional univariate)")
    }
  }
}


########################################################################################################################################################
##' @title Plots the mean or median spider chart of the clusters (multivariate non-functional data)
##'
##' @description This function plots the mean or median spider chart of the clusters. Only for multivariate non-functional data.
##'
##' @param output_clusters list. List of the sites in the clusters: it is the sites_clusters of the output of NPFSS, PFSS, DFFSS, URBFSS, MDFFSS, MRBFSS, MG, MNP, UG or UNP, or the sites_clusters_LH/sites_clusters_W/sites_clusters_P/sites_clusters_R of the MPFSS.
##' @param data matrix. Matrix of the data, the rows correspond to the sites (or to the individuals) and each column represents a variable.
##' @param variable_names character vector. The names of the variables if you want the names to appear. By default NULL: the names are var1, var2, etc.
##' @param type character. "mean" or "median". If "mean": the mean curves in the clusters are plotted in solid lines, outside the cluster in dots, the global mean curve is in black. If "median": the median curves in the clusters are plotted in solid lines, outside the cluster in dots, the global median curve is in black.
##' @param colors character. The colors to plot the clusters' summary curves. If length(colors)==1 then all the clusters will be plotted in this color. Else there must be the same number of elements in colors than the number of clusters
##'
##' @examples
##' \donttest{
##' library(sp)
##' data("map_sites")
##' data("multi_data")
##' coords <- coordinates(map_sites)
##'
##' res_mnp <- MNP(data=multi_data, sites_coord = coords, system = "WGS84",
##' mini = 1, maxi = nrow(coords)/2)
##'
##' plot_summary_chart(output_clusters = res_mnp$sites_clusters,
##' data = multi_data, type = "mean")}
##' \dontshow{
##' library(sp)
##' data("map_sites")
##' data("multi_data")
##' indices <- which(map_sites$NAME_3 == "Lille")
##' coords <- coordinates(map_sites[indices,])
##' res_mnp <- MNP(data=multi_data[indices,], sites_coord = coords,
##' system = "WGS84", mini = 1, maxi = nrow(coords)/2, MC = 99)
##' if(length(res_mnp$sites_clusters)>0){
##' plot_summary_chart(output_clusters = res_mnp$sites_clusters,
##' data = multi_data[indices,], type = "mean")
##' }
##' }
##'
##' @return No value returned, plots the spider chart.
##'
##' @export
##'
plot_summary_chart <- function(output_clusters, data, variable_names = NULL, type = "mean", colors = "red"){

  if(length(colors) == 1){
    if(colors == "black"){
      stop("The color provided cannot be black")
    }
  }else{
    if(length(colors) != length(output_clusters)){
      stop("colors must be of length 1 or with as many colors as the length of output_clusters")
    }else{
      if(sum(colors %in% c("black")) > 0){
        stop("The colors provided cannot be black")
      }
    }
  }

  if(length(type)!=1){
    stop("type must be of length 1")
  }else{
    if(!(type %in% c("mean", "median"))){
      stop("type must be mean or median")
    }
  }

  if(class(data)[1] != "matrix"){
    stop("The data must be a matrix")
  }

  if(is.null(variable_names)){
    variable_names <- paste("var", c(1:ncol(data)), sep = "")
  }else{
    if(length(variable_names) != ncol(data)){
      stop("There must be the same number of elements in variable_names than the number of columns in data")
    }
  }

    if(class(data)[1] == "matrix"){

      if(type == "mean"){
        for(i in 1:length(output_clusters)){
          mean_in <- colMeans(matrix(data[output_clusters[[i]],], ncol = ncol(data)))
          mean_out <- colMeans(matrix(data[-output_clusters[[i]],], ncol = ncol(data)))
          mean_overall <- colMeans(matrix(data, ncol = ncol(data)))

          mini <- colMins(matrix(data, ncol = ncol(data)))
          maxi <- colMaxs(matrix(data, ncol = ncol(data)))

          data_frame <- data.frame(matrix(ncol=ncol(data), nrow=3))

          data_frame[1,] <- mean_in
          data_frame[2,] <- mean_out
          data_frame[3,] <- mean_overall

          colnames(data_frame) <- variable_names
          rownames(data_frame) <- c("mean_in", "mean_out", "mean_overall")

          data_frame <- rbind(maxi, mini, data_frame)

          if(length(colors)==1){
            radarchart(data_frame, pcol = c(colors, colors, "black"), plty = c(1,2,1), title = paste("Spider chart of the means for the cluster", i, sep = " "), vlcex = 0.8, cex.main = 0.8)
            legend("topleft", legend = c("Mean inside the cluster", "Mean outside the cluster", "Global mean"), col = c(colors, colors, "black"), lty = c(1,2,1), bty = "n", cex = 0.8)
          }else{
            radarchart(data_frame, pcol = c(colors[i], colors[i], "black"), plty = c(1,2,1), title = paste("Spider chart of the means for the cluster", i, sep = " "), vlcex = 0.8, cex.main = 0.8)
            legend("topleft", legend = c("Mean inside the cluster", "Mean outside the cluster", "Global mean"), col = c(colors[i], colors[i], "black"), lty = c(1,2,1), bty = "n", cex = 0.8)
          }
        }
      }else{
        for(i in 1:length(output_clusters)){
          med_in <- colMedians(matrix(data[output_clusters[[i]],], ncol = ncol(data)))
          med_out <- colMedians(matrix(data[-output_clusters[[i]],], ncol = ncol(data)))
          med_overall <- colMedians(matrix(data, ncol = ncol(data)))

          mini <- colMins(matrix(data, ncol = ncol(data)))
          maxi <- colMaxs(matrix(data, ncol = ncol(data)))

          data_frame <- data.frame(matrix(ncol=ncol(data), nrow=3))

          data_frame[1,] <- med_in
          data_frame[2,] <- med_out
          data_frame[3,] <- med_overall

          colnames(data_frame) <- variable_names
          rownames(data_frame) <- c("median_in", "median_out", "median_overall")

          data_frame <- rbind(maxi, mini, data_frame)

          if(length(colors)==1){
            radarchart(data_frame, pcol = c(colors, colors, "black"), plty = c(1,2,1), title = paste("Spider chart of the medians for the cluster", i, sep = " "), vlcex = 0.8, cex.main = 0.8)
            legend("topleft", legend = c("Median inside the cluster", "Median outside the cluster", "Global median"), col = c(colors, colors, "black"), lty = c(1,2,1), bty = "n", cex = 0.8)
          }else{
            radarchart(data_frame, pcol = c(colors[i], colors[i], "black"), plty = c(1,2,1), title = paste("Spider chart of the medians for the cluster", i, sep = " "), vlcex = 0.8, cex.main = 0.8)
            legend("topleft", legend = c("Median inside the cluster", "Median outside the cluster", "Global median"), col = c(colors[i], colors[i], "black"), lty = c(1,2,1), bty = "n", cex = 0.8)
          }
        }
      }


    }else{
      stop("The data must be a matrix")
    }

}
