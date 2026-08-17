#' Creates a Venn Diagram based on the kmer sequence distribuition among datasets.
#' The suggested maximum number of datasets is 4 but can hold it up to 7 datasets, if there are more markes it c
#' @param ... A list of at least two dataframes from containing renamed headers and sequences
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#' @param dataset_label A list with the names of the datasets wished to be in the diagram, KEEP THE datasets' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#' @param plot_title The title for the plot.
#'
#' @return A Venn Diagram based on the kmer distribuition among datasets
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
#' venn_k <- dataset_venn_kmers(data_1, data_2)
#' # 3. View result
#' print(venn_k)
dataset_venn_kmers <- function(...,
                               path_db = NULL,
                               dataset_label = NULL,
                               color_pallete = "PuBuGn",
                               plot_title = "K-mer distribution across datasets"){
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

    #Check if dataset quantity is the same as in the list
    if(!is.null(dataset_label)){
      if (num_datasets != length(dataset_label)){
        stop("Error: Number or datasets in dataset and given list in dataset label is not the same")
      }
    }

    #Check if dataset quantity is doable
    if(num_datasets > 4 & num_datasets <= 7){
      message("The suggested maximum number of datasets is 4 but can hold it up to 7 datasets.")
    }

    else if(num_datasets > 7){
      message("If number of datasets is larger than 7, it will give an upset plot")
    }

    #If no label is list given, use the name of the datasets present in the dataframe as labels
    if(is.null(dataset_label)){
      dataset_label <- unique(trimws(unlist(stringr::str_split(exclusive_seqs$dataset, ", "))))
    }

    #Create gradient of collors based on given the pallete
    color_gradient <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, color_pallete))(100)

    #Get the upper limit for gradient
    max_val <- max(sapply(venn_list, length))

    #Create venn diagram
    message("Generating Venn Diagram...\n")
    venn <- ggVennDiagram::ggVennDiagram(venn_list,
                                         edge_size = 0.6,
                                         label_alpha = 0,
                                         label_percent_digit = 0,
                                         category.names = dataset_label) +
      ggplot2::scale_fill_gradientn(colors = color_gradient, limits = c(0, max_val)) +
      ggplot2::labs(title = plot_title)



    return(venn)
}
