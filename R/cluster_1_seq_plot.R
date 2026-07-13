#' Plots the quantity of clusters with 1 k-mer, and the kmer count
#' @param cluster_1_seq_df Dataframe from  the get_clusters_1_seq function.
#' @param marker_label A list with the names of the markers wished to be in the diagram, KEEP THE MARKERS' ORDER IN MIND
#'
#' @return A Venn Diagram based on the cluster distribuition among markers
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' cluster_1_seq_df <- data.frame(cluster = c(1, 2, 3, 4, 5),
#'  V2 = c("A_1_69", "B_1_34", "C_1_35", "B_2_8", "C_2_9"),
#'   V3 = c("A", "B", "C", "B", "C"))
#'
#' # 2. Run function with temporary file
#' plots_1_seq <- cluster_1_seq_plot(cluster_1_seq_df)
#' # 3. View result
#' print(plots_1_seq)
cluster_1_seq_plot <- function(cluster_1_seq_df, marker_label = NULL){
  #Get coverage from kmers
  message("Extracting the k-mer coverage...")
  cov_cd_1_seq <- cluster_1_seq_df |>
    dplyr::rename_with(~"cluster", 1) |>
    dplyr::rename_with(~"id", 2) |>
    dplyr::rename_with(~"marker", 3) |>
    tidyr::separate(.data$id, into = c("marker", "row_n", "cov"), sep = "_") |>
    dplyr::mutate(cov = as.numeric(.data$cov)) |>
    dplyr::mutate(cluster = as.numeric(.data$cluster))

 #Get list of markers
  marker_list <- cov_cd_1_seq |>
    dplyr::select(.data$marker) |>
    dplyr::distinct() |>
    dplyr::pull(.data$marker)

  #Get the quantity of clusters
  message("Counting the quantity of clusters...")
  cov_clusters <- cov_cd_1_seq |>
    dplyr::filter(.data$marker %in% marker_list) |>
    dplyr::group_by(.data$marker, .data$cov) |>
    dplyr::count() |>
    dplyr::rename(cluster_quant = .data$n) |>
    dplyr::ungroup()

  #Set labels
  ##If marker_label is not given uses the markers in the dataframe
  if(is.null(marker_label)){
    corrected_labels = as.list(marker_list)
  }

  #If a list of labels is given it must associate them with the markers
  else if(!is.null(marker_label)){
    #If the number of given labels and markers is not the same give off an error
    if (length(marker_label) != length(marker_list)){
      stop(glue::glue("Error: The quantity of markers given in marker_label must exactly {length(marker_list)}"))

      #Associate the labels with markers
      corrected_labels <- stats::setNames(marker_label, marker_list)
      }
  }

  #Plot
  message("Generating plots...\n")
  plot_1_seq <- ggplot2::ggplot(cov_clusters, ggplot2::aes(x=.data$cluster_quant, y = .data$cov, color = .data$marker)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::facet_wrap(~marker, labeller = ggplot2::as_labeller(corrected_labels)) +
    ggplot2::labs(x = "Cluster quantity", y = "K-mer count", color = "Marker") +
    ggplot2::theme_bw()

  return(plot_1_seq)
}
