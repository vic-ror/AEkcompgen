#' Creates an Euler Diagram based on the cluster distribuition among dataset
#' The suggested maximum number of dataset is 4.
#' @param collapsed_dataset A dataframe from the collapse_datasets function
#' @param dataset_label A list with the names of the datasets wished to be in the diagram, KEEP THE dataset' ORDER IN MIND
#' @param color_pallete A color pallete native from R, default is "Set 2"
#' @param title_name The title for the plot.
#' @return A Euler diagram based on the cluster distribuition among datasets
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' collapsed_dataset_test <- data.frame(cluster = c(1, 2, 3, 4, 5),
#'  V2 = c("A, B, C", "B, C", "B, A", "A, C", "B"))
#'
#' # 2. Run function with temporary file
#' euler_cd <- dataset_euler_clusters(collapsed_dataset_test)
#' # 3. View result
#' print(euler_cd)
dataset_euler_clusters <- function(collapsed_dataset,
                                   dataset_label = NULL,
                                   color_pallete = "Set 2",
                                   title_name = "Cluster distribution across dataset"){

  #Create a list where each dataset has a cluster associated with
  message("Creating lists for diagram based on the given dataset...")
  cluster_list_dataset <- collapsed_dataset |>
    dplyr::rename_with(~ "cluster", 1) |>
    dplyr::rename_with(~ "dataset", 2) |>
    tidyr::separate_rows(.data$dataset, sep = ", ")

  #Create list where each dataset has a cluster associated
  venn_list <- split(cluster_list_dataset$cluster, cluster_list_dataset$dataset)

  #Count dataset quantity
  num_dataset <- length(venn_list)

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
  euler_list <- eulerr::euler(venn_list)

  #If no label is list given, use the name of the dataset present in the dataframe as labels
  if(is.null(dataset_label)){
    dataset_label <- unique(trimws(unlist(stringr::str_split(collapsed_dataset$dataset, ", "))))
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
                   main = title_name)

  return(euler_cd)

}
