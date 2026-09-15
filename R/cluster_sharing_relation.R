#' Joins different datasets labels based on the cluster they belong to, and returns them in list form to generate plots about the sharing relation
#'
#' @param dataframe_cd A dataframe from the run_cd_hit function.
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file.
#' 
#' @return A list of lists containing the the clusters and in which the dataset is present
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#'  collapsed_datasets_test <- data.frame(cluster = c(1, 2, 3, 4, 5),
#'  V2 = c("A, B, C", "B, C", "B, A", "A, C", "B"))
#' # 2. Run function with temporary file
#' cluster_sharing_list <- cluster_sharing_relation(collapsed_datasets_test)
#' # 3. View result
#' print(cluster_sharing_list)
#'
cluster_sharing_relation <- function(dataframe_cd, path_db = NULL){
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
  
  #Create a dataframe collapsing the labels
  message("Collapsing labels...\n")
  cluster_label_joined <- dplyr::tbl(con, "input_table") |>
    dplyr::rename_with(~ "cluster", 1) |>
    dplyr::rename_with(~ "id", 2) |>
    dplyr::mutate(id = stringr::str_remove(.data$id, "_.*")) |>
    dplyr::distinct() |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(label = stringr::str_flatten(.data$id, collapse = ", ")) |>
    dplyr::arrange(.data$cluster) |>
    dplyr::collect()
  
  #Reordering the labels so they are organized in alphabetic order
  cluster_label_joined_reordered <- cluster_label_joined |>
    dplyr::mutate(
      label = purrr::map_chr(
        stringr::str_split(.data$label, ", "),
        ~ stringr::str_c(sort(.x), collapse = ", ")
      )
    )
  
  #Create list 
  message("Creating lists for diagram based on the given datasets...")
  cluster_list_dataset <- cluster_label_joined_reordered |>
    dplyr::rename_with(~ "cluster", 1) |>
    dplyr::rename_with(~ "datasets", 2) |>
    tidyr::separate_rows(.data$datasets, sep = ", ")
  
  #Create list where each dataset has a cluster associated
  cluster_relation_list <- split(cluster_list_dataset$cluster, cluster_list_dataset$datasets)
  
  return(cluster_relation_list)
}
