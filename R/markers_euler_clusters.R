#' Creates an Euler Diagram based on the cluster distribuition among markers
#' The suggested maximum number of markers is 4.
#' @param collapsed_markers A dataframe from the collapse_markers function
#' @param marker_label A list with the names of the markers wished to be in the diagram, KEEP THE MARKERS' ORDER IN MIND
#' @param color_pallete A color pallete native from R, default is "Set 2"
#' @return A Euler diagram based on the cluster distribuition among markers
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' collapsed_markers_test <- data.frame(cluster = c(1, 2, 3, 4, 5),
#'  V2 = c("A, B, C", "B, C", "B, A", "A, C", "B"))
#'
#' # 2. Run function with temporary file
#' euler_cd <- markers_euler_clusters(collapsed_markers_test)
#' # 3. View result
#' print(euler_cd)
markers_euler_clusters <- function(collapsed_markers, marker_label = NULL, color_pallete = "Set 2"){

  #Create a list where each markers has a cluster associated with
  message("Creating lists for diagram based on the given markers...")
  cluster_list_marker <- collapsed_markers |>
    dplyr::rename_with(~ "cluster", 1) |>
    dplyr::rename_with(~ "markers", 2) |>
    tidyr::separate_rows(.data$markers, sep = ", ")

  #Create list where each marker has a cluster associated
  venn_list <- split(cluster_list_marker$cluster, cluster_list_marker$markers)

  #Count marker quantity
  num_markers <- length(venn_list)

  #Check if marker_label (if given) has the same number of markers as the markers in dataset
  if(!is.null(marker_label)){
    if (num_markers != length(marker_label)){
      stop("Error: Number or markers in dataset and given list in marker label is not the same")
    }
  }

  #Check if marker quantity is doable
  if(num_markers > 4){
    stop("Error: Cannot plot more than 7 markers in a Venn diagram. Pleace reduce the number of markers")
  }

  #Fit data to the euler diagram
  euler_list <- eulerr::euler(venn_list)

  #If no label is list given, use the name of the markers present in the dataframe as labels
  if(is.null(marker_label)){
    marker_label <- unique(trimws(unlist(stringr::str_split(collapsed_markers$markers, ", "))))
  }

 #Generate color pallete
  message("Generating color palette based on the number of markers...")

  #Generate num_markers collors dynamically
  dynamic_colors <- grDevices::hcl.colors(n = num_markers, palette = color_pallete)

  #Create euler diagram
  message("Generating Euler Diagram...\n")
  euler_cd <- plot(euler_list,
                   quantities = TRUE,
                   labels = marker_label,
                   fills = list(fill = dynamic_colors, alpha = 0.6),
                   main = "Cluster distribution across markers")

  return(euler_cd)

}
