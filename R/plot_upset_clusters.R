#' Creates a Upset Plot based on the cluster distribuition among datasets.
#'
#' @param cluster_sharing_relation A list from the cluster_sharing_relation function
#' @param dataset_label A list with the names of the datasets wished to be in the diagram, KEEP THE datasetS' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#' @param empty_intersection If the empty intersectios are to be displayed put "on" or not, default is NULL
#'
#' @return An Upset Plot based on the cluster distribuition among datasets
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' collapsed_datasets_test <- data.frame(cluster = c(1, 2, 3, 4, 5),
#'  V2 = c("A, B, C", "B, C", "B, A", "A, C", "B"))
#' # 2. Run cluster_sharing_relation function to get list
#' dataset_sharing_relation <- cluster_sharing_relation(collapsed_datasets_test)
#' # 3. Run function with temporary file
#' upset_cd <- plot_upset_clusters(dataset_sharing_relation)
#' # 3. View result
#' print(upset_cd)
plot_upset_clusters <- function(cluster_sharing_relation,
                                dataset_label = NULL,
                                color_pallete = "PuBuGn",
                                empty_intersection = NULL){
  #Count dataset quantity
  num_dataset <- length(cluster_sharing_relation)

  #Check if dataset_label (if given) has the same number of datasets as the datasets in dataset
  if(!is.null(dataset_label)){
    if (num_dataset != length(dataset_label)){
      stop("Error: Number or datasets in dataset and given list in dataset label is not the same")
    }
    names(cluster_sharing_relation) <- dataset_label
  }

  #Turn list into the upset plot format
  upset_list <- UpSetR::fromList(cluster_sharing_relation)


  #Extract collors from color pallete
  palette_colors <- RColorBrewer::brewer.pal(9, color_pallete)
  dark_color <- palette_colors[8]
  mid_color  <- palette_colors[4]


  #Plotting
  message("Generating UpsetPlot...")
  upset_plot <- UpSetR::upset(data = upset_list,
                              order.by = "freq",
                              nsets = num_dataset,
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
