#' Creates an UpSet Plot based on the cluster distribuition among markers.
#' Recommeded whe working with a lot of markers
#' @param collapsed_markers A dataframe from the collapse_markers function
#' @param marker_label A list with the names of the markers wished to be in the diagram, KEEP THE MARKERS' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#' @param empty_intersection If the empty intersectios are to be displayed put "on" or not, default is NULL
#' @param file_name The name for the file in case you want the plot to be saved automatically
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
#' upset_cd <- markers_upset(collapsed_markers_test)
#' # 3. View result
#' print(upset_cd)
markers_upset <- function(collapsed_markers, marker_label = NULL, color_pallete = "PuBuGn", empty_intersection = NULL, file_name = NULL){

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

  #Check if marker_label (if given) has the same number of markers as the markers in dataset
  if(!is.null(marker_label)){
    if (num_markers != length(marker_label)){
      stop("Error: Number or markers in dataset and given list in marker label is not the same")
    }
  }

  #Turn list into the upset plot format
  upset_list <- UpSetR::fromList(venn_list)

  #If marker labels are given
  if(!is.null(marker_label)){
    names(upset_list) <- marker_label
  }
  #Extract collors from color pallete
  palette_colors <- RColorBrewer::brewer.pal(9, color_pallete)
  dark_color <- palette_colors[8]
  mid_color  <- palette_colors[4]


  #Plotting
  message("Generating UpsetPlot...")
  upset_plot <- UpSetR::upset(data = upset_list,
                order.by = "freq",
                nsets = num_markers,
                point.size = 3.5,
                mainbar.y.label = "Cluster quantity\n",
                sets.x.label = "Clusters per Marker",
                main.bar.color = dark_color,
                sets.bar.color = mid_color,
                matrix.color = dark_color,
                empty.intersections = empty_intersection,
                text.scale=c(1.3, 1.3, 1, 1, 1.3, 1.5)
                )


  #OPCAO DE MUDAR NOME DE TAGSf

  return(upset_plot)

}

