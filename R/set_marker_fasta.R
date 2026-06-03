#' Add a marker to a fasta file
#'
#' @param fasta A fasta file or the path to one.
#' @param marker The tag term that will tag the sequences in this file.
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#'
#' @return A dataframe with the collumns "k-mer_id", "sequence" and "tag".
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary fasta file
#' fasta_temp <- tempfile(fileext = ".fasta")
#' writeLines(
#' text = c(
#' ">seq1", "ATCGATCG",
#' ">seq2", "TCGATCGA"),
#' con = fasta_temp
#' )
#' # 2. Run function with temporary file
#' marked_file <- set_marker_fasta(fasta_temp, "Leish")
#' # 3. View result
#' print(marked_file)
#' # 4. Delete temporary file
#' unlink(fasta_temp)

#Fasta Input
set_marker_fasta <- function(fasta, marker, path_db = NULL) {
  #If no databse path is given, it creates a temporary file
  if (is.null(path_db)) {
    path_db <- tempfile(fileext = ".duckdb")
  }
  # Creates and connects to database
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path_db)

  #Closes the conection after leaving the function
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  #Sets the path of the file
  path_fasta <- fasta

  #Create table on .duckdb
  DBI::dbExecute(
    con,
    paste0("CREATE TABLE fasta_table AS SELECT * FROM read_csv_auto('", path_fasta, "', header = FALSE);")) #Creates dataframe from fasta file

  #List table
  fasta.data <- dplyr::tbl(con, "fasta_table")

  #Add row number
  fasta_data_row <- fasta.data |>
    dplyr::mutate (row_number = dplyr::row_number())

  #Transform fasta into a line file
  #Select header
  kmer_id <- fasta_data_row|>
    dplyr::filter(stringr::str_detect(
      .data$column0, ">")
      ) |>
    dplyr::mutate(column0 = stringr::str_remove(.data$column0, ">")) |>
    dplyr::rename(id = .data$column0)

  #Select sequencia
  kmer_seq <- fasta_data_row |>
    dplyr::filter(stringr::str_detect(
      .data$column0, ">",
      negate = TRUE)
      ) |>
    dplyr::rename(seq = .data$column0) |>
    dplyr::mutate(row_number = .data$row_number - 1)

  #Join header and sequence
  kmer_line <-  dplyr::full_join(kmer_id, kmer_seq, by = "row_number")

  #Add Marker
  kmer_line_marked <- kmer_line |>
    dplyr::mutate(tag = marker) |>
    dplyr::collect() |>
    dplyr::select(.data$id, .data$seq, .data$tag)


  #Print dataframe content
  return(kmer_line_marked)
}

