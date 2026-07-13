#' Selects the clusters that contain only one sequence
#'
#' @param dataframe_cd A dataframe from the run_cd_hit function
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#'
#' @return A dataframe containing the collumns "cluster", "id", "tag".
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' cdhit_out <- data.frame(V1 = c(1, 1, 2, 3, 3),
#'  V2 = c("A_1", "B_13", "C_69", "A_4", "A_2"))
#'
#' # 2. Run function with  dataframe
#' cluster_1_seq <- get_clusters_1_seq(cdhit_out)
#' # 3. View result
#' print(cluster_1_seq)
get_clusters_1_seq <- function(dataframe_cd,
                                   path_db = NULL){
  #If no datablase path is given, creates a temporary file for it
  if (is.null(path_db)) {
    path_db <- tempfile(fileext = ".duckdb")
  }
  #Create connection
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path_db)

  #Guarantee that the conection will be closed after leaving the function
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  #Load table on duckdb
  duckdb::duckdb_register(con, "input_table", dataframe_cd)

  #Work with the given table
  input_id_marker <- dplyr::tbl(con, "input_table") |>
    dplyr::rename_with(~"cluster", 1) |>
    dplyr::rename_with(~"id", 2) |>
    dplyr::mutate(marker = stringr::str_remove(.data$id, "_.*"))


    only_one_seq <- input_id_marker |>
      dplyr::group_by(.data$cluster) |>
      dplyr::tally() |>
      dplyr::filter(.data$n == 1)

    only_1_seq_info <- dplyr::inner_join(input_id_marker,
                                            only_one_seq,
                                            by = "cluster") |>
      dplyr::select(.data$cluster, .data$id, .data$marker) |>
      dplyr::arrange(.data$cluster) |>
      dplyr::collect()

  return(only_1_seq_info)
}
