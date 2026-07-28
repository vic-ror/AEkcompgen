#' Plots an histogram based on the k-mer count and frequency
#' @param ... Dataframes from the mark_jellyfish histo function
#' @param normalization_value A value to divide the frequency of k-mers, usually based on the coverage of the mapped reads.
#' @param colors The color or list of colors to fill the plot's bars.
#'
#' @return A plot objetct displaying a k-mer frequency hostogram:
#'  \item{X-axis}{K-mer frequency} the number of times a k-mer was counted.
#'  \item{Y-axis}{Count} Quantity of distinct k-mers that appeared the number of times indicated in the x axis.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframes
#' data_frame_test <- data.frame(id = c("1", "2", "3"),
#' seq = c("100", "150", "50"),
#' tag = c("1", "2", "2"))
#' # 2. Run function with temporary file
#' histogram <- plot_histogram_kmer_freq(data_frame_test)
#' # 3. View result
#' print(histogram)

plot_histogram_kmer_freq <- function(..., normalization_value = NULL, colors = NULL){

  #Create list with dataframes
  histo_list <- list(...)

  #Merge dataframes
  binded_histo_list <- dplyr::bind_rows(histo_list) |>
    dplyr::rename_with(~"kmer_freq", 1) |>
    dplyr::rename_with(~"count", 2) |>
    dplyr::rename_with(~"tag", 3)

  #Check if file is empty
  if (length(binded_histo_list) == 0){
    stop("Error: The file(s) given are empty.")
  }

  #If a normalization value is given
  if(!is.null(normalization_value)){
    # Plot histogram
    histogram <- ggplot2::ggplot(binded_histo_list,  ggplot2::aes(x = as.numeric(.data$kmer_freq)/normalization_value, y = as.numeric(.data$count), fill = .data$tag)) +
      ggplot2::geom_area(alpha=0.8) +
      ggplot2::scale_y_log10(labels = scales::label_scientific()) +
      ggplot2::labs(x = "Kmer frequency", y = "Count", title="K-mer count histogram") +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~ .data$tag) +
      ggplot2::guides(fill= ggplot2::guide_legend(title="Marker"))


  }

  else{
    # Plot histogram
    histogram <- ggplot2::ggplot(binded_histo_list,  ggplot2::aes(x = as.numeric(.data$kmer_freq), y = as.numeric(.data$count), fill = .data$tag)) +
      ggplot2::geom_area(alpha=0.8) +
      ggplot2::scale_y_log10(labels = scales::label_scientific()) +
      ggplot2::labs(x = "Kmer frequency", y = "Count", title="K-mer count histogram") +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~ .data$tag) +
      ggplot2::guides(fill= ggplot2::guide_legend(title="Marker"))
  }

  #Apply colors if given any
  if (!is.null(colors)) {
    histogram <- histogram + ggplot2::scale_fill_manual(values = colors)
  }
  return(histogram)
}

