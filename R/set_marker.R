#' Add a marker to a fasta file
#'
#' @param fasta A fasta file or the path to one.
#' @param marker The tag term that will tag the sequences in this file.
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#'
#' @return A dataframe with the collumns "k-mer_id", "sequence" and "tag".
#' @export
#'
#' @examples
#' marked_file <- set_marker("../fasta.fasta.txt", "Leish")

#Fasta Input
set_marker <- function(fasta, marker, path_db = NULL) {
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
  fasta.data <- tbl(con, "fasta_table")

  #Add row number
  fasta_data_row <- fasta.data |> mutate (row_number = row_number())

  #Transform fasta into a line file
  #Select header
  kmer_id <- fasta_data_row|>
    filter(str_detect(
      column0, ">")
      ) |>
    mutate(column0 = str_remove(column0, ">")) |>
    rename(id = column0)
  #Select sequencia
  kmer_seq <- fasta_data_row |>
    filter(str_detect(
      column0, ">",
      negate = TRUE)
      ) |>
    rename(seq = column0) |>
    mutate(row_number = row_number - 1)
  #Join header and sequence
  kmer_line <- full_join(kmer_id, kmer_seq, by = "row_number")

  #Add Marker
  kmer_line_marked <- kmer_line |>
    mutate(tag = marker) |>
    collect() |>
    select(id, seq, tag)


  #Save dataframe in the environment
  return(kmer_line_marked)
}

