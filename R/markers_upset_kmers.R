#' Creates an Upset Plot based on the kmer sequence distribuition among markers.
#' The suggested maximum number of markers is 4 but can hold it up to 7 markers, if there are more markes it c
#' @param ... A list of at least two dataframes from containing renamed headers and sequences
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#' @param marker_label A list with the names of the markers wished to be in the diagram, KEEP THE MARKERS' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#' @param empty_intersection If the empty intersectios are to be displayed put "on" or not, default is NULL
#' @param file_name The name for the file in case you want the plot to be saved automatically
#'
#' @return An Upset plot based on the kmer distribuition among markers
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
#' venn_k <- markers_upset_kmers(data_1, data_2)
#' # 3. View result
#' print(venn_k)
markers_upset_kmers <- function(..., path_db = NULL, marker_label = NULL, color_pallete = "PuBuGn", empty_intersection = NULL, file_name = NULL){
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


  #Check if marker_label (if given) has the same number of markers as the markers in dataset
  if(!is.null(marker_label)){
    if (num_markers != length(marker_label)){
      stop("Error: Number or markers in dataset and given list in marker label is not the same")
    }
  }

  #Turn list into the upset plot format
  upset_list <- UpSetR::fromList(venn_list)

  #If marker labels are given
  if(!is.null(marker_label)){
    names(upset_list) <- marker_label
  }
  #Extract collors from color pallete
  palette_colors <- RColorBrewer::brewer.pal(9, color_pallete)
  dark_color <- palette_colors[8]
  mid_color  <- palette_colors[4]


  #Plotting
  message("Generating UpsetPlot...")
  upset_plot <- UpSetR::upset(data = upset_list,
                              order.by = "freq",
                              nsets = num_markers,
                              point.size = 3.5,
                              mainbar.y.label = "Kmer quantity\n",
                              sets.x.label = "Kmers per Marker",
                              main.bar.color = dark_color,
                              sets.bar.color = mid_color,
                              matrix.color = dark_color,
                              empty.intersections = empty_intersection,
                              text.scale=c(1.3, 1.3, 1, 1, 1.3, 1.5)
  )


    return(upset_plot)

}
