#' Runs the Jellyfish histo function
#' @param jf_file The .jf output file from jellyfish
#' @param output The name of the file in which the histogram will be saved
#'
#' @return A histogram containing the k-mer frequency and count of the k-mers in the .jf file
#' @export
#' @importFrom rlang .data
#'
#' @examplesIf Sys.which("jellyfish") != ""
#' # 1. Load .jf test file
#' temp_jf <- system.file("extdata", "test.jf", package = "AEkcompgen")
#' # 2. Run function with temporary file
#' histo_file <- run_jellyfish_histo(temp_jf)
#' # 3. View result
#' print(histo_file)
run_jellyfish_histo <- function(jf_file, output = NULL){
  #Check if jellyfish is installed
  if (Sys.which("jellyfish") == "") {
    stop(
      "ERROR: The program 'jellyfish' was not found in the system/n
      Please check if it's installed and add it to your system's path",
      call. = FALSE
    )
  }
  #Check if .jf file is available in the given path
  #Check if file exists in the given path
  if(!file.exists(jf_file)){
    stop("ERROR: .jf file not found in the given path.")
  }

  #Running jellyfish histo
  system(glue::glue("jellyfish histo {jf_file} > temp.histo"))

  #Loading the histo file and adding label
  histo_path <- ("temp.histo")

  histo_file <- readr::read_lines(histo_path, lazy = FALSE)

  histo_df <- tibble::tibble(histo_file) |>
    tidyr::separate(histo_file, sep = " ", into = c("kmer_freq", "count"))

  #If no name for output file is given delete it
  if(is.null(output)){
      if(file.exists("temp.histo")) unlink("temp.histo")
  }

  #If output name is given, rename the output as the file name
  else{
    if(file.exists("temp.histo")) file.rename(from = "temp.histo", to = glue::glue("{output}"))
  }

  return(histo_df)
}
