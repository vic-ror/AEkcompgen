#' Retrieves the sequence information from a dataframe containing only the id and cluster information.
#'
#' @param clustered_df A dataframe containing the Cluster and id, k-mer identification, collumns.
#' @param seq_df A dataframe containing the seq, sequence of the k-mer, and id, kmer identification, collumns.
#' @param fasta_file_name If want to save a fasta file containing the k-mers and their cluster info, give the name of the file.
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file.
#'
#' @return A dataframe containing the Cluster, id and seq information.
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # 1. Create temporary fasta file
#' data_frame_cluster <- data.frame(cluster = c(1, 1, 2),
#'  id = c("dataset1_1_28", "dataset2_1_33", "dataset1_2_2"))
#' data_frame_seq <- data.frame(id = c("dataset1_1_28", "dataset2_1_33", "dataset1_2_2"),
#' seq = c("GACAGGTACAAGAAGGAGTA", "AGGGCGACCTTCGATTCGGA", "TTTACACACTCTCCTTGGAC"))
#' # 2. Run function with temporary file
#' seq_info_cluster_df <- retrieve_sequence(data_frame_cluster, data_frame_seq)
#' # 3. View resulting line dataframe
#' print(seq_info_cluster_df)
retrieve_sequence <- function(clustered_df,
                              seq_df,
                              fasta_file_name = NULL,
                              path_db = NULL){

  #Creates a temporary path to save connection if no path is given
  if (is.null(path_db)) {
    path_db <- tempfile(fileext = ".duckdb")
  }
  #Creates and connects to database
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path_db)

  #Guarantee exiting the connection after leaving the function
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  #Load table on duckdb
  duckdb::duckdb_register(con, "cluster_table", clustered_df) # clustered df
  duckdb::duckdb_register(con, "seq_table", seq_df) # seq df

  #Join data
  joined_df <- dplyr::full_join(dplyr::tbl(con, "cluster_table"), dplyr::tbl(con, "seq_table"), by = "id")

  if(!is.null(fasta_file_name)){
    fasta_cluster <- joined_df |>
      dplyr::mutate(header=paste0(">", .data$id, " | ", .data$cluster)) |>
      dplyr::select(.data$header, .data$seq) |>
      dplyr::collect() |>
      tidyr::pivot_longer(cols = c(.data$header, .data$seq),
                               names_to = "type",
                               values_to = "fasta_line") |>
      dplyr::select(.data$fasta_line)

    utils::write.table(fasta_cluster, file = glue::glue("{fasta_file_name}.fasta"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }

  complete_df <- joined_df |>
    dplyr::collect()

  return(complete_df)
  }
