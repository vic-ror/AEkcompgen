#' Creates an Euler Diagram based on the kmer sequence distribuition among datasets.
#' The suggested maximum number of datasets is 4 but can hold it up to 7 datasets, if there are more markes it c
#' @param ... A list of at least two dataframes from containing renamed headers and sequences
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#' @param dataset_label A list with the names of the datasets wished to be in the diagram, KEEP THE datasetS' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#' @param title_name The title for the plot.
#'
#' @return An Euler Diagram based on the kmer distribuition among datasets
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
#' venn_k <- dataset_euler_kmers(data_1, data_2)
#' # 3. View result
#' print(venn_k)
dataset_euler_kmers <- function(...,
                                path_db = NULL,
                                dataset_label = NULL,
                                color_pallete = "PuBuGn",
                                title_name = "Kmer distribution across datasets"){
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
  venn_list <- split(exclusive_seqs$seq, exclusive_seqs$dataset)

  #Count dataset quantity
  num_datasets <- length(venn_list)

  #Check if dataset_label (if given) has the same number of datasets as the datasets in dataset
  if(!is.null(dataset_label)){
    if (num_datasets != length(dataset_label)){
      stop("Error: Number or datasets in dataset and given list in dataset label is not the same")
    }
  }

  #Check if dataset quantity is doable
  if(num_datasets > 4){
    stop("Error: Cannot plot more than 7 datasets in a Venn diagram. Pleace reduce the number of datasets")
  }

  #Fit data to the euler diagram
  euler_list <- eulerr::euler(venn_list)

  #If no label is list given, use the name of the datasets present in the dataframe as labels
  if(is.null(dataset_label)){
    dataset_label <- unique(trimws(unlist(stringr::str_split(exclusive_seqs$dataset, ", "))))
  }

  #Generate color pallete
  message("Generating color palette based on the number of datasets...")

  #Generate num_datasets collors dynamically
  dynamic_colors <- grDevices::hcl.colors(n = num_datasets, palette = color_pallete)

  #Create euler diagram
  message("Generating Euler Diagram...\n")
  euler_cd <- plot(euler_list,
                   quantities = TRUE,
                   labels = dataset_label,
                   fills = list(fill = dynamic_colors, alpha = 0.6),
                   main = title_name)

  return(euler_cd)

}
