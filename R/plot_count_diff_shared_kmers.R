#' Compares the kmer frequency between two markers
#' Takes output of the select_shared_seq and compares the kmer frequency between two markers giving the modulus
#' of the difference.
#' @param shared_df Dataframe from select_shared_seq.
#' @param cut_off A minimal value of kmer_frequency
#' @param separate_markers Separates markers into different plots, using facet_wrap in ggplot.
#' @param colors  A list of two colors to fill the plot's bars.
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file.
#'
#' @return A daframe with the frequency modulus difference.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' shared_test <- data.frame(seq = c("GAGA", "ATGA", "AAAA"),
#'  marker = c("Marker1,Marker2", "Marker1,Marker2", "Marker1,Marker2"),
#'   header_1 = c("Marker1_3_3", "Marker1_2_15", "Marker1_3_25"),
#'   header_2 = c("Marker2_3_45", "Marker2_2_50", "Marker2_3_30"))
#'
#' # 2. Run function with temporary file
#' plots_freq_comp <- plot_count_diff_shared_kmers(shared_test)
#' # 3. View result
#' print(plots_freq_comp)
plot_count_diff_shared_kmers <- function(shared_df, cut_off = NULL, separate_markers = FALSE, colors = c("purple", "green"), path_db = NULL) {

  #If no datablase path is given, creates a temporary file for it
  if (is.null(path_db)) {
    path_db <- tempfile(fileext = ".duckdb")
  }
  #Create connection
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path_db)

  #Guarantee that the conection will be closed after leaving the function
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  #Load table on duckdb
  duckdb::duckdb_register(con, "shared_table", shared_df)

  #Load_shared seqs
  count_diff_info <- dplyr::tbl(con, "shared_table") |>
    dplyr::select(.data$header_1, .data$header_2, .data$seq) |>
    dplyr::mutate(
      marker_1 = dbplyr::sql("split_part(header_1, '_', 1)"),
      count_1 = dbplyr::sql("split_part(header_1, '_', 3)"),
      marker_2 = dbplyr::sql("split_part(header_2, '_', 1)"),
      count_2 = dbplyr::sql("split_part(header_2, '_', 3)")
    ) |>
    dplyr::select(.data$marker_1, .data$marker_2, .data$count_1, .data$count_2, .data$seq) |>
    dplyr::mutate(count_diff = as.numeric(.data$count_1)-as.numeric(.data$count_2))

  #Filter the sequences with a count difference
  count_diff_higher_marker_1 <- count_diff_info |>
    dplyr::filter(.data$count_diff > 0) |>
    dplyr::mutate(higher_marker = .data$marker_1) |>
    dplyr::collect()

  count_diff_higher_marker_2 <- count_diff_info |>
    dplyr::filter(.data$count_diff < 0) |>
    dplyr::mutate(higher_marker = .data$marker_2) |>
    dplyr::collect()


  if(is.null(cut_off)) {
    #Join dataframes to plot them together
    joined_count_diff <- rbind(count_diff_higher_marker_1, count_diff_higher_marker_2) |>
      dplyr::select(.data$seq, .data$count_diff, .data$higher_marker) |>
      dplyr::rename(Marker = .data$higher_marker)
  }

  else if(!is.null(cut_off)) {
    #Filter values by its cut off
    count_diff_higher_marker_1_cut_off <- count_diff_higher_marker_1 |>
      dplyr::filter(.data$count_diff >= cut_off)

    count_diff_higher_marker_2_cut_off <- count_diff_higher_marker_2 |>
      dplyr::filter(.data$count_diff <= -cut_off)

    #Join
    joined_count_diff <- rbind(count_diff_higher_marker_1_cut_off, count_diff_higher_marker_2_cut_off) |>
      dplyr::select(.data$seq, .data$count_diff, .data$higher_marker) |>
      dplyr::rename(Marker = .data$higher_marker)
  }


  #Not eparate markers
  if(!isTRUE(separate_markers)) {
    #Plot
    plot_freq <- ggplot2::ggplot(data = joined_count_diff, ggplot2::aes(x = .data$seq,
                                                                        y = .data$count_diff,
                                                                        fill = .data$Marker)) +
      ggplot2::geom_bar(stat = "identity",
                        width = 0.7) +
      ggplot2::labs(x = "", y = "Modulus of Count difference", title = "Modulus of Count Difference") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "white"),  #Remove background
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::scale_fill_manual(values = colors)
  }

  #If want to separate markers
  else if(isTRUE(separate_markers)) {
    #Make all the values the absolute value
    joined_count_diff_abs <- joined_count_diff |>
      dplyr::mutate(count_diff = abs(.data$count_diff))

    #Plot
    plot_freq <- ggplot2::ggplot(data = joined_count_diff_abs, ggplot2::aes(x = .data$seq,
                                                                            y = .data$count_diff,
                                                                            fill = .data$Marker)) +
      ggplot2::geom_bar(stat = "identity",
                        width = 0.7) +
      ggplot2::facet_wrap(~.data$Marker, scales = "free_x") +
      ggplot2::labs(x = "", y = "Modulus of Count difference", title = "Modulus of Count Difference") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(color = "white"),  #Remove background
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::scale_fill_manual(values = colors)
  }

  return(plot_freq)
}
