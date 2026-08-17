#' Creates an UpSet Plot based on the cluster distribuition among datasets.
#' Recommeded whe working with a lot of datasets
#' @param collapsed_datasets A dataframe from the collapse_datasets function
#' @param dataset_label A list with the names of the datasets wished to be in the diagram, KEEP THE datasetS' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn".
#' @param empty_intersection If the empty intersectios are to be displayed put "on" or not, default is NULL
#' @param file_name The name for the file in case you want the plot to be saved automatically
#'
#' @return A Venn Diagram based on the cluster distribuition among datasets
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' collapsed_datasets_test <- data.frame(cluster = c(1, 2, 3, 4, 5),
#'  V2 = c("A, B, C", "B, C", "B, A", "A, C", "B"))
#'
#' # 2. Run function with temporary file
#' upset_cd <- dataset_upset_clusters(collapsed_datasets_test)
#' # 3. View result
#' print(upset_cd)
dataset_upset_clusters <- function(collapsed_datasets,
                                   dataset_label = NULL,
                                   color_pallete = "PuBuGn",
                                   empty_intersection = NULL,
                                   file_name = NULL){

  #Create lists for ggVennDiagram
  message("Creating lists for diagram based on the given datasets...")
  cluster_list_dataset <- collapsed_datasets |>
    dplyr::rename_with(~ "cluster", 1) |>
    dplyr::rename_with(~ "datasets", 2) |>
    tidyr::separate_rows(.data$datasets, sep = ", ")

  #Create list where each dataset has a cluster associated
  venn_list <- split(cluster_list_dataset$cluster, cluster_list_dataset$datasets)

  #Count dataset quantity
  num_datasets <- length(venn_list)

  #Check if dataset_label (if given) has the same number of datasets as the datasets in dataset
  if(!is.null(dataset_label)){
    if (num_datasets != length(dataset_label)){
      stop("Error: Number or datasets in dataset and given list in dataset label is not the same")
    }
  }

  #Turn list into the upset plot format
  upset_list <- UpSetR::fromList(venn_list)

  #If dataset labels are given
  if(!is.null(dataset_label)){
    names(upset_list) <- dataset_label
  }
  #Extract collors from color pallete
  palette_colors <- RColorBrewer::brewer.pal(9, color_pallete)
  dark_color <- palette_colors[8]
  mid_color  <- palette_colors[4]


  #Plotting
  message("Generating UpsetPlot...")
  upset_plot <- UpSetR::upset(data = upset_list,
                order.by = "freq",
                nsets = num_datasets,
                point.size = 3.5,
                mainbar.y.label = "Cluster quantity\n",
                sets.x.label = "Clusters per dataset",
                main.bar.color = dark_color,
                sets.bar.color = mid_color,
                matrix.color = dark_color,
                empty.intersections = empty_intersection,
                text.scale=c(1.3, 1.3, 1, 1, 1.3, 1.5)
                )



  return(upset_plot)

}

