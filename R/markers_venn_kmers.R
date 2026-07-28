#' Creates a Venn Diagram based on the kmer sequence distribuition among markers.
#' The suggested maximum number of markers is 4 but can hold it up to 7 markers, if there are more markes it c
#' @param ... A list of at least two dataframes from containing renamed headers and sequences
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#' @param marker_label A list with the names of the markers wished to be in the diagram, KEEP THE MARKERS' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#'
#' @return A Venn Diagram based on the kmer distribuition among markers
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframes
#' data_1 <- data.frame(header = c("tag1_1_1", "tag1_2_2", "tag1_3_3", "tag1_4_4", "tag1_5_5"),
#'  seq = c("AAAG", "AGAG", "CTGC", "AAAA", "GAGA"))
#' data_2 <- data.frame(header = c("tag2_1_1", "tag2_2_2", "tag2_3_3", "tag2_4_4", "tag2_5_5"),
#'  seq = c("AATG", "AGAG", "CAGC", "AATA", "GAGA"))
#'
#' # 2. Run function with temporary file
#' venn_k <- markers_venn_kmers(data_1, data_2)
#' # 3. View result
#' print(venn_k)
markers_venn_kmers <- function(..., path_db = NULL, marker_label = NULL, color_pallete = "PuBuGn"){
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

  #Detect markers
  message("Extracting markers...")
  exclusive_seqs <- dplyr::tbl(con, "merged_table") |>
    dplyr::mutate(marker = stringr::str_remove(.data$header, "_.*")) |>
    dplyr::select(.data$seq, .data$marker) |>
    dplyr::collect()

  #Creating list for venn
  message("Creating lists...")
    #Create list where each marker has a cluster associated
    venn_list <- split(exclusive_seqs$seq, exclusive_seqs$marker)

    #Count marker quantity
    num_markers <- length(venn_list)

    #Check if marker quantity is the same as in the list
    if(!is.null(marker_label)){
      if (num_markers != length(marker_label)){
        stop("Error: Number or markers in dataset and given list in marker label is not the same")
      }
    }

    #Check if marker quantity is doable
    if(num_markers > 4 & num_markers <= 7){
      message("The suggested maximum number of markers is 4 but can hold it up to 7 markers.")
    }

    else if(num_markers > 7){
      message("If number of markers is larger than 7, it will give an upset plot")
    }

    #If no label is list given, use the name of the markers present in the dataframe as labels
    if(is.null(marker_label)){
      marker_label <- unique(trimws(unlist(stringr::str_split(exclusive_seqs$marker, ", "))))
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
                                         category.names = marker_label) +
      ggplot2::scale_fill_gradientn(colors = color_gradient, limits = c(0, max_val)) +
      ggplot2::labs(title = "K-mer distribution across markers")



    return(venn)
}
