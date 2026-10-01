#' Selects the shared sequences of two files after the run_jellyfish or rename_fasta_jf functions based on the dataset label
#'
#' @param  ...  At least two dataframes containing the header and the seq column
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#'
#' @return A dataframe containing the exclusive sequences to each label
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary fasta file
#' fasta_temp1 <- data.frame(id = c("seq1", "seq2", "seq3", "seq4", "seq5"),
#'  seq = c("ATCGATCG", "TCGATCGA", "TCAGTCAG", "AATATATA", "TAGGGT"))
#'
#' fasta_temp2 <- data.frame(id = c("seq6", "seq7", "seq8", "seq9", "seq10"),
#'  seq = c("GAGAGAG", "TAGGGT", "TCTCTCTC", "AATATATA", "ATCACTG"))
#' # 3. Run function with example dataframes
#' exc_seq <- select_shared_seq(fasta_temp1, fasta_temp2)
#' # 4. View result
#' print(exc_seq)
select_shared_seq <- function(..., path_db = NULL) {
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
  shared_seqs_label <- dplyr::tbl(con, "merged_table") |>
    dplyr::group_by(.data$seq) |>
    dplyr::filter(dplyr::n() > 1) |> #Filter sequences that appear more than once
    dplyr::mutate(label = stringr::str_remove(.data$id, "_.*")) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$seq) |>
    dplyr::summarise(label = stringr::str_flatten(.data$label, collapse = ","))

  shared_seqs_header <- dplyr::tbl(con, "merged_table") |>
    dplyr::group_by(.data$seq) |>
    dplyr::filter(dplyr::n() > 1) |> #Filter sequences that appear more than once
    dplyr::ungroup() |>
    dplyr::group_by(.data$seq) |>
    dplyr::summarise(id = stringr::str_flatten(.data$id, collapse = ","))

  shared_seqs <- dplyr::full_join(shared_seqs_label, shared_seqs_header, by = "seq") |>
    dplyr::collect() |>
    tidyr::separate_longer_delim(.data$id, delim = ",") |>
    dplyr::mutate(extracted_label = stringr::str_remove(.data$id, "_.*")) |>
    dplyr::group_by(.data$seq) |>
    dplyr::arrange(.data$extracted_label, .by_group = TRUE) |>
    dplyr::mutate(col_id = paste0("id_", dplyr::row_number())) |>
    dplyr::mutate(label = stringr::str_flatten(.data$extracted_label, collapse = ",")) |>
    dplyr::select(-.data$extracted_label) |>
    dplyr::distinct() |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = .data$col_id, values_from = .data$id)

  #Information on quantity of clusters
  shared_seqs_quantity <- shared_seqs |>
    nrow()

  message(glue::glue("{shared_seqs_quantity} shared k-mers found!"))


  #Reordering the labels so they are organized in alphabetic order
  return(shared_seqs)
}

