#' Gives info of the run_cd_hit output.
#' Such as quantity of clusters, mean, median, largest cluster and smallest cluster.
#'
#' @param dataframe_cd A dataframe from the run_cd_hit function
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#' @param out The name or path of a file to save the information on a .txt file
#'
#' @return A dataframe containing the collumns "cluster", "id", "tag".
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframe
#' cdhit_out <- data.frame(V1 = c(1, 1, 2, 3, 3),
#'  V2 = c("infantum_1_28", "major_1_33", "infantum_2_2", "infantum_3_19", "major_2_49"))
#'
#' # 2. Run function with dataframe
#' marked_file <- cd_hit_info(cdhit_out)
#' # 3. View result
#' print(marked_file)


cd_hit_info <- function(dataframe_cd, path_db = NULL, out = NULL){
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

  input.data <- dplyr::tbl(con, "input_table") |> ##
    dplyr::rename_with(~ "cluster", 1) |>
    dplyr::rename_with(~ "id", 2)

  kmers_per_cluster <- input.data |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::arrange(.data$cluster) |>
    dplyr::collect()


  #Calculating stats of cd_hit_output
  #number of clusters
  number_clusters <- nrow(kmers_per_cluster)

  #Mean of the quantity of kmers in clusters
  kmers_per_cluster_mean <- mean(kmers_per_cluster$n)

  #Median of the quantity of kmers in cluster
  kmers_per_cluster_median <- stats::median(kmers_per_cluster$n)

  #Largest cluster
  #Largest_quantity_kmers
  largest_kmers <- max(kmers_per_cluster$n)
  largest_clusters <- kmers_per_cluster |>
    dplyr::filter(.data$n >= largest_kmers) |>
    dplyr::select(.data$cluster) |>
    as.list()

  #Smallest cluster
  smallest_kmers <- min(kmers_per_cluster$n)
  smallest_clusters <- kmers_per_cluster |>
    dplyr::filter(.data$n <= smallest_kmers) |>
    dplyr::select(.data$cluster) |>
    as.list()

  info_cdhit <- glue::glue("Number of clusters: {number_clusters}
                     Mean quantity of k-mers per cluster: {kmers_per_cluster_mean}
                     Median quantity of k-mers per cluster: {kmers_per_cluster_median}
                     Largest cluster(s) quantity of k-mers: {largest_kmers}\
                     Smallest cluster(s) quantity of k-mer(s): {smallest_kmers}")

  if(!is.null(out)){
    writeLines(info_cdhit,
               glue::glue("{out}.txt"))
  }

  return(info_cdhit)
}
