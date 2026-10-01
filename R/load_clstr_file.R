#' Load clstr file from CD-hit in environment
#' @param file A clstr file from CD-hit.
#'
#' @return A dataframe containing the clustered groups
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # 1. Load clstr file path
#' clstr_path <- system.file("extdata",
#'                           "git_hub_test",
#'                           "clustered_datasets.clstr",
#'                            package = "AEkcompgen")
#' # 2. Run load_clstr_file
#' clstr_df <- load_clstr_file(clstr_path)
#' # 3. View resulting line dataframe
#' print(clstr_df)
#'

load_clstr_file <- function(file){
  #Turn clster output file into a R dataframe
  message("Turning .clstr file into a dataframe...")

  #Load .clstr file
  clstr_file <- readr::read_lines(file, lazy = FALSE)

  if (length(clstr_file) == 0) {
    stop("Error: The clstr file is empty.")
  }

  #Identify cluster and kmers
  clstr_df <- tibble::tibble(clstr_file) |>
    dplyr::mutate(is_cluster_id = stringr::str_detect(.data$clstr_file, "^>")) |>
    dplyr::mutate(group_id = cumsum(.data$is_cluster_id))

  cluster <- clstr_df |>
    dplyr::filter(.data$is_cluster_id == TRUE) |>
    dplyr::rename(cluster = .data$clstr_file) |>
    dplyr::select(.data$cluster, .data$group_id) |>
    dplyr::mutate(cluster = stringr::str_remove(.data$cluster, ">"))

  kmer <- clstr_df |>
    dplyr::filter(.data$is_cluster_id == FALSE) |>
    dplyr::rename(id = .data$clstr_file) |>
    dplyr::mutate(id = stringr::str_remove(.data$id, ".*>")) |>
    dplyr::mutate(id = stringr::str_remove(.data$id, "\\.\\.\\..*")) |>
    dplyr::select(.data$id, .data$group_id)

  joined_info_clst <- dplyr::full_join(cluster, kmer, by = "group_id", multiple = "all") |>
    dplyr::select(.data$cluster, .data$id)

  return(joined_info_clst)
}
