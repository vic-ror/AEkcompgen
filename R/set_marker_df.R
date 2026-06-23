#' Add a marker to a dataframe
#'
#' @param df A dataframe.
#' @param marker The tag term that will tag the sequences in this file.
#'
#' @return A dataframe with the collumns "k-mer_id", "sequence" and "tag".
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' dataframe_test <- data.frame(id = c("seq_1", "seq_2"),
#'  seq = c("ATCGATCG", "TCGATCGA"))
#' # 2. Run function with temporary file
#' marked_file <- set_marker_df(dataframe_test, "Leish")
#' # 3. View result
#' print(marked_file)

#Fasta Input
set_marker_df <- function(df, marker) {
  #Add Marker
  kmer_line_marked <- df |>
    dplyr::mutate(tag = marker) |>
    dplyr::select(.data$id, .data$seq, .data$tag)

  #Print dataframe content
  return(kmer_line_marked)
}

