#' Selects the exclusive sequences of two files after the set_marker() function
#'
#' @param  dataset_1 Dataframe from set_marker() to be compared with dataset_2
#' @param dataset_2 Dataframe from set_marker() to be compared with dataset_1
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file

#'
#' @return A dataframe containing the exclusive sequences to each marker
#' @export
#'
#' @examples
#' # 1. Create temporary fasta file
#' fasta_temp1 <- tempfile(fileext = ".fasta")
#' writeLines(c(">seq1", "ATCGATCG",
#' ">seq2", "TCGATCGA",
#' ">seq3", "TCAGTCAG"),
#' con = fasta_temp1
#' )
#'
#' fasta_temp2 <- tempfile(fileext = ".fasta")
#' writeLines(c(">seq1", "GCGCGCGC",
#' ">seq2", "TATATATA",
#' ">seq3", "TCAGTCAG"),
#' con = fasta_temp2
#' )
#' # 2. Run set_marker
#' marked_df_1 <- set_marker(fasta_temp1, "1")
#' marked_df_2 <- set_marker(fasta_temp2, "2")
#'
#' # 3. Run function with temporary file
#' exc_seq <- select_exclusive_seq(marked_df_1, marked_df_2)
#' # 4. View result
#' print(exc_seq)
#' # 5. Delete temporary file
#' unlink(fasta_temp1)
#' unlink(fasta_temp2)

select_exclusive_seq <- function(dataset_1, dataset_2, path_db = NULL) {
  #If no databse path is given, it creates a temporary file
  if (is.null(path_db)) {
    path_db <- tempfile(fileext = ".duckdb")
  }
  # Creates and connects to database
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path_db)

  #Closes the conection after leaving the function
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

#load tables on duckdb
duckdb::duckdb_register(con, "table_1", dataset_1)
duckdb::duckdb_register(con, "table_2", dataset_2)

table_1.data <- dplyr::tbl(con, "table_1")
table_2.data <- dplyr::tbl(con, "table_2")

#Select exclusive sequences
table_1_exc <- table_1.data |>
  dplyr::anti_join(table_2.data,
                   by = "seq")
table_2_exc <- table_2.data |>
  dplyr::anti_join(table_1.data,
                   by = "seq")

tables_merged <- dplyr::union_all(table_1_exc, table_2_exc) |>  dplyr::collect()


#Save dataframe in the environment
return(tables_merged)
}
