#' Creates a Euler Diagram based on the cluster distribuition among datasets.
#' The suggested maximum number of datasets is 4 but can hold it up to 7 datasets, if there are more datasets the diagram may not make sense
#' @param cluster_sharing_relation A list from the cluster_sharing_relation function
#' @param dataset_label A list with the names of the datasets wished to be in the diagram, KEEP THE datasetS' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#' @param plot_title The title for the plot.
#'
#' @return An Euler Diagram based on the cluster distribuition among datasets
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
#' euler_cd <- plot_euler_clusters(dataset_sharing_relation)
#' # 3. View result
#' print(euler_cd)
plot_euler_clusters <- function(cluster_sharing_relation,
                               dataset_label = NULL,
                               color_pallete = "PuBuGn",
                               plot_title = "Cluster distribution across datasets"){
  #Count dataset quantity
  num_dataset <- length(cluster_sharing_relation)

  #Check if dataset_label (if given) has the same number of dataset as the dataset in dataset
  if(!is.null(dataset_label)){
    if (num_dataset != length(dataset_label)){
      stop("Error: Number or dataset in dataset and given list in dataset label is not the same")
    }
  }

  #Check if dataset quantity is doable
  if(num_dataset > 4){
    stop("Error: Cannot plot more than 7 dataset in a Venn diagram. Pleace reduce the number of dataset")
  }

  #Fit data to the euler diagram
  euler_list <- eulerr::euler(cluster_sharing_relation)

  #If no label is list given, use the name of the dataset present in the dataframe as labels
  if(is.null(dataset_label)){
    dataset_label <- names(cluster_sharing_relation)
  }

  #Generate color pallete
  message("Generating color palette based on the number of dataset...")

  #Generate num_dataset collors dynamically
  dynamic_colors <- grDevices::hcl.colors(n = num_dataset, palette = color_pallete)

  #Create euler diagram
  message("Generating Euler Diagram...\n")
  euler_cd <- plot(euler_list,
                   quantities = TRUE,
                   labels = dataset_label,
                   fills = list(fill = dynamic_colors, alpha = 0.6),
                   main = plot_title)

  return(euler_cd)

}
