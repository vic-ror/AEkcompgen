#' Subsample the given input bam/sam file to the wanted average sequencing depth.
#' @param input A sorted bam/sam file, (run_bwa already sorts its output)
#' @param current_average_depth The current average sequencing depth of the input file.
#' @param wanted_average_depth The wanted average sequencing depth.
#' @param output The path for the subsampled bam/sam file.
#' @param output_format Format for the output of the mapping file (bam/sam), default is "bam' for binary file.
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file.
#'
#' @return Saves the subsampled file in the output given path, and prints the current average depth calculated
#' @export
#' @importFrom rlang .data
#'
#' @examplesIf Sys.which("samtools") != ""
#' #1. Load .bam test file
#' temp_bam <- system.file("extdata", "test_mapped.bam", package = "AEkcompgen")
#' #2. Run subsample_sequencing function
#' subsample_sequencing(input =  temp_bam,
#'  current_average_depth = 1,
#'  wanted_average_depth = 0.8,
#'  output = "test_mapped_subsampled",
#'  output_format = "bam")
#' #3. Delete file
#' unlink("test_mapped_subsampled.bam")
subsample_sequencing <- function(input,
                              current_average_depth,
                              wanted_average_depth,
                              output,
                              output_format = "bam",
                              path_db = NULL){
  #Check if samtools is installed
  if (Sys.which("samtools") == "") {
    stop(
      "ERROR: The program 'samtools' was not found in the system\n
      Please check if it's installed and add it to your system's path or environment\n
      Recommendation: Install it via Conda",
      call. = FALSE
    )
  }

  #Check if the output_format is bam or sam
  if(output_format != "sam"  && output_format !="bam"){
    stop("ERROR: Output file extension given must be 'sam' or 'bam'.")
  }

  #Check if given depht are numeric
  if(!is.numeric(current_average_depth) &&  !is.numeric(wanted_average_depth)){
    stop("ERROR: Depths given should be numeric values.")
  }

  #Calculate relation value
  message("Calculating relation value to subsample...")
  relation_value <- wanted_average_depth/current_average_depth

  #Run samtools collate to collate paired reads and samtools view to do subsample
  #-O to rpeserve the order of the reads
  #-h to keep header
  #-s to subsample
  if(output_format == "sam"){
    #If sam
    message(glue::glue("Collating paired reads and subsampling for depth coverage equal to {wanted_average_depth}.."))
    system(glue::glue("samtools collate -O {input} | samtools view -h -s {relation_value} | samtools sort -o {output}.sam"))
  } else if (output_format == "bam"){   #If waanted output is bam
    message(glue::glue("Collating paired reads and subsampling for depth coverate equal to {wanted_average_depth}.."))
    #-b to output bam
    system(glue::glue("samtools collate -O {input} | samtools view -h -b -s {relation_value} | samtools sort -o {output}.bam"))
  }

  #Check if subsample was sucessfull
  message("Calculating new file depth to verify if subsample was sucessfull..")
  #Run samtools depth to calculate read depth
  system2(command = "samtools", args = c("depth", glue::glue("{output}.bam")), stdout = glue::glue("{input}_temp_depth.csv"))

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
    paste0("CREATE TABLE input_table AS SELECT * FROM read_csv_auto('", input, "_temp_depth.csv", "', header = FALSE);")) #Creates dataframe from file


  #Load dataframe with depth info
  depth_df <- dplyr::tbl(con, "input_table")

  #Calculate average genome depth
  average_depth <- depth_df |>
    dplyr::summarise(mean_depth = mean(.data$column2, na.rm = TRUE)) |>
    dplyr::pull(.data$mean_depth)

  #Remove csv file
  unlink(glue::glue("{input}_temp_depth.csv"))


  message(glue::glue("New average sequencing depth: {average_depth}"))
  }
