#' Calculates the average sequencing depth of a bam/sam file.
#'
#' @param input A sorted bam/sam file, (run_bwa already sorts its output)
#' @param keep_csv If wanted to keep csv file with each read depth.
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file.
#'
#' @return The average sequencing depth value of the input bam/sam file.
#' @export
#' @importFrom rlang .data
#'
#' @examplesIf Sys.which("samtools") != ""
#' #1. Load .bam test file
#' temp_bam <- system.file("extdata", "test_mapped.bam", package = "AEkcompgen")
#' #2. Run caclulate depth function
#' calculate_depth(input =  temp_bam)

calculate_depth <- function(input,
                            keep_csv = FALSE,
                            path_db = NULL){
  if (Sys.which("samtools") == "") {
    stop(
      "ERROR: The program 'samtools' was not found in the system\n
      Please check if it's installed and add it to your system's path or environment\n
      Recommendation: Install it via Conda",
      call. = FALSE
    )
  }

  #output file name
  file_name <- stringr::str_remove(glue::glue("{input}"), "\\..*")

  #Run samtools depth to calculate read depth
  system2(command = "samtools", args = c("depth", input), stdout = glue::glue("{file_name}_depth.csv"))

  #Creates a temporary path to save connection if no path is given
  if (is.null(path_db)) {
    path_db <- tempfile(fileext = ".duckdb")
  }
  #Creates and connects to database
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path_db)

  #Guarantee exiting the connection after leaving the function
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  #Load table on duckdb
  DBI::dbExecute(
    con,
    paste0("CREATE TABLE input_table AS SELECT * FROM read_csv_auto('", file_name, "_depth.csv", "', header = FALSE);")) #Creates dataframe from file


  #Load dataframe with depth info
  depth_df <- dplyr::tbl(con, "input_table")

  #Calculate average genome depth
  average_depth <- depth_df |>
    dplyr::summarise(mean_depth = mean(.data$column2, na.rm = TRUE)) |>
    dplyr::pull(.data$mean_depth)

  #If not keep csv remove the created csv
  if(isFALSE(keep_csv)){
    unlink(glue::glue("{file_name}_depth.csv"))
  }

  return(message(glue::glue("Average sequencing depth: {average_depth}")))
  }
