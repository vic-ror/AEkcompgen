#' Runs jellyfish and gives it's output
#' The id for the sequence contains the row of the kmer in the fasta file and the the count of the k-mer
#'
#' @param fasta_file A fasta file or the path to one.
#' @param length The length of kmers that will be counted
#' @param marker The marker that will be used in the header of sequence and file names
#' @param lower_count Don't output k-mers with count lower than this
#' @param upper_count Don't output k-mers with count higher than this
#' @param histo Make histogram based on kmer frequency
#' @param save_fasta If you want to save the intermediate files or not
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#' 
#' @return A dataframe containing the counted kmers and their headers with the frequency information
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # 1. Create temporary fasta file
#' fasta_temp <- tempfile(fileext = ".fasta")
#' writeLines(
#' text = c(
#' ">seq1", "GACAGGTACAAGAAGGAGTA",
#' ">seq2", "AGGGCGACCTTCGATTCGGA",
#' ">seq3", "TTTACACACTCTCCTTGGAC",
#' ">seq4", "TGTGAACTTTTAAATTCGAT",
#' ">seq5", "CACTTAAGGCTTGAAAACTA",
#' ">seq6", "GACGAGACACACTCTTCAAG",
#' ">seq7", "AGCTGGCTAAGCGCGCGCGC",
#' ">seq8", "GTTAACGGCAGGTGCAACCC",
#' ">seq9", "CCGTCAGCGGCCAGGATAGG",
#' ">seq10", "TCTGCGGGAGCTCTGCGGG"
#' ),
#' con = fasta_temp
#' )
#' # 2. Run function with temporary file
#' line_file <- run_jellyfish(fasta_temp, path_db =,length =  10, marker = "test", histo = FALSE, save_fasta = FALSE)
#' # 3. View resulting line dataframe
#' print(line_file)
#' # 4. Delete temporary file
#' unlink(fasta_temp)

run_jellyfish <- function(fasta_file, length, marker, lower_count = NULL, upper_count = NULL, histo = TRUE, path_db = NULL, save_fasta = FALSE){
  #Check if jellyfish and seqkit is installed
  if (Sys.which("jellyfish") == "") {
    stop(
      "ERROR: The program 'jellyfish' was not found in the system/n
      Please check if it installed and add it to your system's path",
      call. = FALSE
    )
  }
  if (Sys.which("seqkit") == "") {
    stop(
      "ERROR: The program 'seqkit' was not found in the system/n
      Please check if it installed and add it to your system's path",
      call. = FALSE
    )
  }
  
  #Run_jellyfish
  #With lower count filter = Prints only kmers with count higher or equal to the filter
  if (!is.null(lower_count) && is.null(upper_count)){
    message(glue::glue("Counting k-mers with count higher or equal to {lower_count}..."))
    system(glue::glue("jellyfish count -m {length} -L {lower_count} -s 275M -t 10 -C -o {marker}.jf {fasta_file}"))
  }
  #With higher count filter = Prints only kmers with count smaller or equal to the filter
  else if (is.null(lower_count) && !is.null(upper_count)){
    message(glue::glue("Counting k-mers with count smaller or equal to {upper_count}..."))
    system(glue::glue("jellyfish count -m {length} -U {upper_count} -s 275M -t 10 -C -o {marker}.jf {fasta_file}"))
  }
  #With higher and lower count filter = Prints only kmers with count higher, smaller or equal to the filter
  else if (!is.null(lower_count) && !is.null(upper_count)){
    message(glue::glue("Counting k-mers with count smaller or equal to {upper_count} and higher or equal to {lower_count}..."))
    system(glue::glue("jellyfish count -m {length} -U {upper_count} -L {lower_count} -s 275M -t 10 -C -o {marker}.jf {fasta_file}"))
  }
  #Without count filter
  else {
    message("Counting k-mers without any filters...")
    system(glue::glue("jellyfish count -m {length} -s 275M -t 10 -C -o {marker}.jf {fasta_file}"))
  }
  #If histo is true creates a file to generate histogram
  if(histo){
    message("Creating file to generate histogram...")
  if (!dir.exists("analyses")) {
    dir.create("analyses", recursive = TRUE) #Create directory
  }
    system(glue::glue("jellyfish histo {marker}.jf > analyses/{marker}.histo")) #Histo archive
    }

  message("Creating readable fasta file...")
  system(glue::glue("jellyfish dump {marker}.jf > temp"))# turn output into fasta
  message("Renaming fasta file headers to contain the row number and the count of the kmer...")
  system(glue::glue("seqkit replace -p '(.+)' -r '{marker}_mf_{{nr}}_$1' temp > {marker}_kmers.fasta")) #Rename from fasta headers
  system("rm temp")
  
  #Turn fasta into line
    #If no databse path is given, it creates a temporary file
  message("Creating dataframe...")
  if (is.null(path_db)) {
    path_db <- tempfile(fileext = ".duckdb")
  }
  
  # Creates and connects to database
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path_db)

  #Closes the conection after leaving the function
  on.exit({
    DBI::dbDisconnect(con, shutdown = TRUE)
    
    if (!save_fasta) {
      if (file.exists(glue::glue("{marker}_kmers.fasta"))) file.remove(glue::glue("{marker}_kmers.fasta"))
    }
    })

  #Sets the path of the file
  path_fasta <- glue::glue("{marker}_kmers.fasta")

  #Create table on .duckdb
  DBI::dbExecute(
    con,
    paste0("CREATE TABLE fasta_table AS SELECT * FROM read_csv_auto('", path_fasta, "', header = FALSE);")) #Creates dataframe from fasta file

  #Closes the conection after leaving the function
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  
  #List table
  fasta.data <- dplyr::tbl(con, "fasta_table")

  #Add row number
  fasta_data_row <- fasta.data |>
    dplyr::mutate (row_number = dplyr::row_number())

  #Transform fasta into a line file
  #Select header
  kmer_id <- fasta_data_row|>
    dplyr::filter(stringr::str_detect(
      .data$column0, ">")) |>
    dplyr::mutate(column0 = stringr::str_remove(.data$column0, ">")) |>
    dplyr::rename(id = .data$column0)

  #Select sequencia
  kmer_seq <- fasta_data_row |>
    dplyr::filter(stringr::str_detect(
      .data$column0, ">",
    negate = TRUE)) |>
    dplyr::rename(seq = .data$column0) |>
    dplyr::mutate(row_number = .data$row_number - 1)

  #Join header and sequence
  kmer_line <-  dplyr::full_join(kmer_id, kmer_seq, by = "row_number") |>
    dplyr::collect() |>
    dplyr::select(.data$id, .data$seq)

  #Print dataframe content
  return(kmer_line)
}
