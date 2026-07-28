#' Creates a Venn Diagram based on the cluster distribuition among markers.
#' The suggested maximum number of markers is 4 but can hold it up to 7 markers, if there are more markes it c
#' @param collapsed_markers A dataframe from the collapse_markers function
#' @param marker_label A list with the names of the markers wished to be in the diagram, KEEP THE MARKERS' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#'
#' @return A Venn Diagram based on the cluster distribuition among markers
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' collapsed_markers_test <- data.frame(cluster = c(1, 2, 3, 4, 5),
#'  V2 = c("A, B, C", "B, C", "B, A", "A, C", "B"))
#'
#' # 2. Run function with temporary file
#' venn_cd <- markers_venn_clusters(collapsed_markers_test)
#' # 3. View result
#' print(venn_cd)
markers_venn_clusters <- function(collapsed_markers, marker_label = NULL, color_pallete = "PuBuGn"){

  #Create lists for ggVennDiagram
  message("Creating lists for diagram based on the given markers...")
  cluster_list_marker <- collapsed_markers |>
    dplyr::rename_with(~ "cluster", 1) |>
    dplyr::rename_with(~ "markers", 2) |>
    tidyr::separate_rows(.data$markers, sep = ", ")

  #Create list where each marker has a cluster associated
  venn_list <- split(cluster_list_marker$cluster, cluster_list_marker$markers)

  #Count marker quantity
  num_markers <- length(venn_list)

  #Check if marker quantity is the same as in the list
  if(!is.null(marker_label)){
    if (num_markers != length(marker_label)){
      stop("Error: Number or markers in dataset and given list in marker label is not the same")
    }
  }

  #Check if marker quantity is doable
  if(num_markers > 4 & num_markers <= 7){
    message("The suggested maximum number of markers is 4 but can hold it up to 7 markers.")
  }

  else if(num_markers > 7){
    message("If number of markers is larger than 7, it will give an upset plot")
  }

  #If no label is list given, use the name of the markers present in the dataframe as labels
  if(is.null(marker_label)){
    marker_label <- unique(trimws(unlist(stringr::str_split(collapsed_markers$markers, ", "))))
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
                                       category.names = marker_label) +
    ggplot2::scale_fill_gradientn(colors = color_gradient, limits = c(0, max_val)) +
    ggplot2::labs(title = "Cluster distribution across markers")



  return(venn)

}
