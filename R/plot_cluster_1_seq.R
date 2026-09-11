#' Plots the quantity of clusters with 1 k-mer, and the kmer count
#' @param cluster_1_seq_df Dataframe from  the get_clusters_1_seq function.
#' @param dataset_label A list with the names of the labels wished to be in the diagram, KEEP THE labelS' ORDER IN MIND
#' @param free_scale If the plots by labels should be free scale, put TRUE, or in the same scale, default.
#' @param color_palette A color palette from the scale_fill_gradientn function, default is "Set2".
#' @param plot_title Change the title of the plot.
#'
#' @return A Venn Diagram based on the cluster distribuition among labels
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
#' plots_1_seq <- plot_cluster_1_seq(cluster_1_seq_df)
#' # 3. View result
#' print(plots_1_seq)
plot_cluster_1_seq <- function(cluster_1_seq_df,
                               dataset_label = NULL,
                               free_scale = FALSE,
                               color_palette = "Set2",
                               plot_title = "Not clustered k-mer distribution"){
  #Get counterage from kmers
  message("Extracting the k-mer count...")
  count_cd_1_seq <- cluster_1_seq_df |>
    dplyr::rename_with(~"cluster", 1) |>
    dplyr::rename_with(~"id", 2) |>
    dplyr::rename_with(~"label", 3) |>
    tidyr::separate(.data$id, into = c("label", "row_n", "count"), sep = "_") |>
    dplyr::mutate(count = as.numeric(.data$count)) |>
    dplyr::mutate(cluster = as.numeric(.data$cluster))


 #Get list of labels
  label_list <- count_cd_1_seq |>
    dplyr::select(.data$label) |>
    dplyr::distinct() |>
    dplyr::pull(.data$label)


  #Get the quantity of clusters
  message("Counting the quantity of clusters...")
  count_clusters <- count_cd_1_seq |>
    dplyr::filter(.data$label %in% label_list) |>
    dplyr::group_by(.data$label, .data$count) |>
    dplyr::count() |>
    dplyr::rename(cluster_quant = .data$n) |>
    dplyr::ungroup()

  #Set labels
  ##If dataset_label is not given uses the labels in the dataframe
  if(is.null(dataset_label)){
    corrected_labels = as.list(label_list)
  }

  #If a list of labels is given it must associate them with the labels
  else if(!is.null(dataset_label)){
    #If the number of given labels and labels is not the same give off an error
    if (length(dataset_label) != length(label_list)){
      stop(glue::glue("Error: The quantity of labels given in dataset_label must exactly {length(label_list)}"))

      #Associate the labels with labels
      corrected_labels <- stats::setNames(dataset_label, label_list)
      }
  }

  #Plot
  if (isFALSE(free_scale)){
    message("Generating plots...\n")
    plot_1_seq <- ggplot2::ggplot(count_clusters, ggplot2::aes(y =.data$cluster_quant, x = .data$count, color = .data$label)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 0.8) +
      ggplot2::facet_wrap(~label, labeller = ggplot2::as_labeller(corrected_labels)) +
      ggplot2::labs(x = "K-count", y = "K-mer quantity", color = "Label", title = plot_title) +
      ggplot2::theme_bw() +
      ggplot2::scale_color_brewer(palette = color_palette)

  }

  if (isTRUE(free_scale)){
    message("Generating plots...\n")
    plot_1_seq <- ggplot2::ggplot(count_clusters, ggplot2::aes(y =.data$cluster_quant, x = as.numeric(.data$count), color = .data$label)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 0.8) +
      ggplot2::facet_wrap(~label, labeller = ggplot2::as_labeller(corrected_labels), scales = "free_y") +
      ggplot2::labs(x = "K-mer count", y = "K-mer quantity", color = "Label", title = plot_title) +
      ggplot2::theme_bw() +
      ggplot2::scale_color_brewer(palette = color_palette)

  }

  return(plot_1_seq)
}
