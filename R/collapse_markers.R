#' Joins different markers in a single line, based on the cluster they belong to.
#'
#' @param dataframe_cd A dataframe from the run_cd_hit function
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#'
#' @return A dataframe containing the collumns "cluster" and "markers"
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' cdhit_out <- data.frame(V1 = c(1, 1, 2, 3, 3),
#'  V2 = c("infantum_1_28", "major_1_33", "infantum_2_2", "infantum_3_19", "major_2_49"))
#'
#' # 2. Run function with temporary file
#' collapsed_file <- collapse_markers(cdhit_out)
#' # 3. View result
#' print(collapsed_file)

collapse_markers <- function(dataframe_cd, path_db = NULL){
  #Creates a temporary path to save connection if no path is given
  if (is.null(path_db)) {
    path_db <- tempfile(fileext = ".duckdb")
  }
  #Creates and connects to database
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path_db)

  #Guarantee exiting the connection after leaving the function
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  #Load table on duckdb
  duckdb::duckdb_register(con, "input_table", dataframe_cd)

  #Create a dataframe collapsing the markers
  message("Collapsing markers...\n")
  cluster_marker_joined <- dplyr::tbl(con, "input_table") |>
    dplyr::rename_with(~ "cluster", 1) |>
    dplyr::rename_with(~ "id", 2) |>
    dplyr::mutate(id = stringr::str_remove(.data$id, "_.*")) |>
    dplyr::distinct() |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(markers = stringr::str_flatten(.data$id, collapse = ", ")) |>
    dplyr::arrange(.data$cluster) |>
    dplyr::collect()

  #Reordering the markers so they are organized in alphabetic order
  cluster_marker_joined_reordered <- cluster_marker_joined |>
    dplyr::mutate(
      markers = purrr::map_chr(
        stringr::str_split(.data$markers, ", "),
        ~ stringr::str_c(sort(.x), collapse = ", ")
      )
    )

  return(cluster_marker_joined_reordered)
}
