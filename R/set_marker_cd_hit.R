#' Add a marker to a cd-hit output
#'
#' @param input A file from cd_hit or formated like one.
#' @param markers A list of tag term that will tag the clusters within this file.
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#'
#' @return A dataframe containing the collumns "cluster", "id", "tag".
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary fasta file
#' cdhit_temp <- tempfile(fileext = ".txt")
#' writeLines(
#' text = c(
#' "1  infantum_1",
#' "1  major_1",
#' "2  infantum_2",
#' "3  infantum_3",
#' "3  major_2"),
#' con = cdhit_temp
#' )
#' # 2. Run function with temporary file
#' marked_file <- set_markers_cd_hit(cdhit_temp, c("infantum", "major"))
#' # 3. View result
#' print(marked_file)
#' # 4. Delete temporary file
#' unlink(cdhit_temp)

#Set function
set_markers_cd_hit <- function(input,
                        markers = character(),
                        path_db = NULL){

  #If no databse path is given, it creates a temporary file
  if (length(markers) == 0) {
    stop("This function requires at least one marker")
  }
  if (is.null(path_db)) {
    path_db <- tempfile(fileext = ".duckdb")
  }
  #Creates and connects to database
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path_db)

  #Closes the conection after leaving the function
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  #Creates table from cd_hit file
  test_cd_hit <- utils::read.table(input)

  #load tables on duckdb
 duckdb::duckdb_register(con, "input_table", test_cd_hit)

input.data <- dplyr::tbl(con, "input_table") |>
  dplyr::mutate(tag = NA_character_)


  #Loop for each marker for the given input
  for (marker in markers){
    #Add markers to a file
       input.data <- input.data |>
      dplyr::mutate(tag= dplyr::case_when(
        stringr::str_detect(.data$V2, marker) ~ marker, #When V2 contains the marker, it adds it to the tag
        !is.na(.data$tag) ~ .data$tag,                              #When it was already assigned a marker keep it
        TRUE ~ NA_character_                            #When it still is na keeps na
      ))
  }

 input_data_marked <- input.data |>
   dplyr::collect() |>
   dplyr::mutate(tag = tidyr::replace_na(.data$tag, "OTHER"))

  return(input_data_marked)
  }
