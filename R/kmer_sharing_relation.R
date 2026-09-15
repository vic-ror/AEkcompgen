#' Creates a Venn Diagram based on the kmer sequence distribuition among datasets.
#' The suggested maximum number of datasets is 4 but can hold it up to 7 datasets, if there are more markes it c
#' @param ... A list of at least two dataframes from containing renamed headers and sequences
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file.
#'
#' @return A list of lists containing the the sequence and in which datasets the kmer is present
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframes
#' data_1 <- data.frame(id = c("tag1_1_1", "tag1_2_2", "tag1_3_3", "tag1_4_4", "tag1_5_5"),
#'  seq = c("AAAG", "AGAG", "CTGC", "AAAA", "GAGA"))
#' data_2 <- data.frame(id = c("tag2_1_1", "tag2_2_2", "tag2_3_3", "tag2_4_4", "tag2_5_5"),
#'  seq = c("AATG", "AGAG", "CAGC", "AATA", "GAGA"))
#'
#' # 2. Run function with temporary file
#' dataset_sharing_relation <- kmer_sharing_relation(data_1, data_2)
#' # 3. View result
#' print(dataset_sharing_relation)
kmer_sharing_relation <- function(...,
                                  path_db = NULL){
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

  #Detect datasets
  message("Extracting datasets...")
  exclusive_seqs <- dplyr::tbl(con, "merged_table") |>
    dplyr::mutate(dataset = stringr::str_remove(.data$id, "_.*")) |>
    dplyr::select(.data$seq, .data$dataset) |>
    dplyr::collect()

  #Creating list for venn
  message("Creating lists...")
  #Create list where each dataset has a cluster associated
  dataset_sharing_list <- split(exclusive_seqs$seq, exclusive_seqs$dataset)

  if(length(data_frame_list) > 4){
    message("For better visualization upset plot is better recommended \n Venn and Euler diagrams can get confusing with more than 4 datasets")
  }

  return(dataset_sharing_list)
}
