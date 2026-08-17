#' Creates a Venn Diagram based on the cluster distribuition among datasets.
#' The suggested maximum number of datasets is 4 but can hold it up to 7 datasets, if there are more markes it c
#' @param collapsed_datasets A dataframe from the collapse_datasets function
#' @param dataset_label A list with the names of the datasets wished to be in the diagram, KEEP THE datasetS' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#' @param plot_title The title for the plot.
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
#' venn_cd <- dataset_venn_clusters(collapsed_datasets_test)
#' # 3. View result
#' print(venn_cd)
dataset_venn_clusters <- function(collapsed_datasets,
                                  dataset_label = NULL,
                                  color_pallete = "PuBuGn",
                                  plot_title = "Cluster distribution across datasets"){

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

  #Check if dataset quantity is the same as in the list
  if(!is.null(dataset_label)){
    if (num_datasets != length(dataset_label)){
      stop("Error: Number or datasets in dataset and given list in dataset label is not the same")
    }
  }

  #Check if dataset quantity is doable
  if(num_datasets > 4 & num_datasets <= 7){
    message("The suggested maximum number of datasets is 4 but can hold it up to 7 datasets.")
  }

  else if(num_datasets > 7){
    message("If number of datasets is larger than 7, it will give an upset plot")
  }

  #If no label is list given, use the name of the datasets present in the dataframe as labels
  if(is.null(dataset_label)){
    dataset_label <- unique(trimws(unlist(stringr::str_split(collapsed_datasets$datasets, ", "))))
  }

  #Create gradient of collors based on given the pallete
  color_gradient <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, color_pallete))(100)

  #Get the upper limit for gradient
  max_val <- max(sapply(venn_list, length))

  #Create venn diagram
  message("Generating Venn Diagram...\n")
  venn <- ggVennDiagram::ggVennDiagram(venn_list,
                                       edge_size = 0.6,
                                       label_alpha = 0,
                                       label_percent_digit = 0,
                                       category.names = dataset_label) +
    ggplot2::scale_fill_gradientn(colors = color_gradient, limits = c(0, max_val)) +
    ggplot2::labs(title = plot_title)



  return(venn)

}
