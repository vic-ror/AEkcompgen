#' Compares the kmer frequency between different possible pairs of dataset labels,
#' Takes output of the select_shared_seq and compares the kmer frequency between two labels giving the modulus
#' of the difference.
#' @param shared_df Dataframe from select_shared_seq.
#' @param cut_off A minimal value of k-mer count difference.
#' @param color_palette A color pallete from the scale_fill_gradientn function, default is "Set2".
#' @param free_scale If the plots by labels should be free scale, put TRUE, or in the same scale, default.
#' @param plot_title The title for the plot.
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file.
#'
#' @return A daframe with the frequency modulus difference of different pairs of dataset labels.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' shared_test <- data.frame(seq = c("GAGA", "ATGA", "AAAA"),
#'  label = c("label1,label2", "label1,label2", "label1,label2"),
#'   id_1 = c("label1_3_3", "label1_2_15", "label1_3_25"),
#'   id_2 = c("label2_3_45", "label2_2_50",
#'    "label2_3_30"))
#'
#' # 2. Run function with temporary file
#' plots_freq_comp <- plot_count_diff_shared_kmers(shared_test)
#' # 3. View result
#' print(plots_freq_comp)

plot_count_diff_shared_kmers <- function(shared_df,
                                         cut_off = NULL,
                                         color_palette = "Set2",
                                         free_scale = FALSE,
                                         plot_title = "Modulus of Count Difference",
                                         path_db = NULL) {


    #Get each header in its own row to extract each kmer count
  message("Extracting k-mer count values...")
  kmer_count_info <- shared_df |>
    tidyr::pivot_longer(cols = dplyr::starts_with("id_"),
                        names_to = "header_id",
                        values_to = "header_val",
                        values_drop_na = TRUE) |>
    dplyr::mutate(kmer_count = as.numeric(stringr::str_split_i(.data$header_val, "_", 3)),
                  marker_name = stringr::str_split_i(.data$header_val, '_', 1))


  #Group the sewuences and combine the marker names and the kmer count
  message("Paring different dataset labels combination and calculating the modulus of the difference of the k-mer count...")
  diff_results <- kmer_count_info |>
    dplyr::group_by(.data$seq) |>
    dplyr::summarise(
      pair = list(utils::combn(.data$marker_name, 2, paste, collapse = " - ")),
      diff = list(utils::combn(.data$kmer_count, 2, function(x) abs(x[1] - x[2]))),
      .groups = "drop"
    ) |>
    tidyr::unnest(cols = c(.data$pair, .data$diff))


  if(isTRUE(free_scale)){
    free_scale = "free_y"
  }  else{
    free_scale = "fixed"
  }

  #If cut off was given filter the kmers that the count difference reached the cut off
  if(!is.null(cut_off)){
    #If cut off is not a numeric value
    if(!is.numeric(cut_off)){
      stop("ERROR: cut-off should be a single numeric value, or NULL")
    }

    diff_results <- diff_results |>
      dplyr::filter(.data$diff >= cut_off)
  }

  #If no kmer met the cut off give off an error
  if (nrow(diff_results) == 0) {
    stop("ERROR: No k-mer count difference met the cut_off.")
  }

    message("Making the plot...")
    plot_freq <- ggplot2::ggplot(data = diff_results, ggplot2::aes(x = .data$seq,
                                                                   y = .data$diff,
                                                                   fill = .data$pair)) +
      ggplot2::geom_bar(stat = "identity",
                        width = 0.7) +
      ggplot2::facet_wrap(~.data$pair, scales = free_scale) +
      ggplot2::labs(x = "", y = "Modulus of Count difference", title = plot_title) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "white"),  #Remove background
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::scale_fill_brewer(palette = color_palette)


  return(plot_freq)

}
