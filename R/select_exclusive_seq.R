#' Selects the exclusive sequences of two files after the set_marker() function
#'
#' @param  ...  At least two dataframes containing the header and the seq column
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#'
#' @return A dataframe containing the exclusive sequences to each marker
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary fasta file
#' fasta_temp1 <- data.frame(header = c("seq1", "seq2", "seq3", "seq4", "seq5"),
#'  seq = c("ATCGATCG", "TCGATCGA", "TCAGTCAG", "AATATATA", "TAGGGT"))
#'
#' fasta_temp2 <- data.frame(header = c("seq6", "seq7", "seq8", "seq9", "seq10"),
#'  seq = c("GCGCGCG", "GAGAGAG", "TCTCTCTC", "AATATATA", "TAGGGT"))
#' # 3. Run function with example dataframes
#' exc_seq <- select_exclusive_seq(fasta_temp1, fasta_temp2)
#' # 4. View result
#' print(exc_seq)

select_exclusive_seq <- function(..., path_db = NULL){
  #Load List
  data_frame_list <- list(...)

  #Stop function if at least two dataframes were not given
  if(length(data_frame_list) < 2){
    stop("Error: This function requires at least two dataframes")
  }

  #If no databse path is given, it creates a temporary file
  if (is.null(path_db)) {
    path_db <- tempfile(fileext = ".duckdb")
  }
  # Creates and connects to database
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path_db)

  #Closes the conection after leaving the function
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))


  #Create loop to append each datraframe into DUCKDB
  for (df in seq_along(data_frame_list)){
    table_name <- paste0("temp_df_", df)

    #Register each dataframe
    duckdb::duckdb_register(con, name = table_name, df = data_frame_list[[df]])

    #If its the first dataframe creates the main table
    if (df == 1){
      DBI::dbExecute(con, sprintf("CREATE TABLE merged_table AS SELECT * FROM %s", table_name))
    }

    #Insert other dataframes into the main table
    else{
      DBI::dbExecute(con, sprintf("INSERT INTO merged_table SELECT * FROM %s", table_name))
    }
    #Unregister tables
    duckdb::duckdb_unregister(con, name = table_name)

  }

  exclusive_seqs <- dplyr::tbl(con, "merged_table") |>
    dplyr::group_by(.data$seq) |>
    dplyr::filter(dplyr::n() == 1) |> #Filter sequences that appear only once = Exclusive
    dplyr::ungroup() |>
    dplyr::collect()

  #Unregister table

  return(exclusive_seqs)

}
