#' Adds a marker collumn to a dataframe, or loaded file of a histogram
#' @param histo_file A histogram file from jellyfish histo
#' @param histo_df A loaded dataframe from the run_jellyfish histo function
#' @param marker The marker tag that will be added to the dataframe
#' @param out The name of the file in which the histogram will be saved
#'
#' @return A marked histogram containing the k-mer frequency and count of the k-mers in the .jf file
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframes
#' temp_file <- tempfile(fileext = ".histo")
#' writeLines(
#' text = c(
#' "1 100",
#' "2 150",
#' "3 50"),
#' con = temp_file
#' )
#'
#' # 2. Run function with temporary file
#' marked_histo <- mark_jellyfish_histo(histo_file = temp_file, marker = "temp")
#' # 3. View result
#' print(marked_histo)
#' unlink(temp_file)
mark_jellyfish_histo <- function(histo_file = NULL, histo_df = NULL, marker, out = NULL){
  #Check if files were correctly provided
  #If both files are provided
  if (!is.null(histo_file) && !is.null(histo_df)){
    stop("Please provide EITHER a .histo file OR a loaded .histo in a dataframe, not both.")
  }

  #If no file is provided
  else if (is.null(histo_file) && is.null(histo_df)){
    stop("Please provide A .histo file OR a loaded .histo in a dataframe.")
  }

  #If a .jf file is given
  else if (!is.null(histo_file) && is.null(histo_df)){
    histo <- readr::read_lines(histo_file, lazy = FALSE)

    renamed_histo_df <- tibble::tibble(histo) |>
      tidyr::separate(histo, sep = " ", into = c("kmer_freq", "count")) |>
      dplyr::mutate(tag = marker)

  }
  #If a dataframe is given add a marker collumn
  else if (is.null(histo_file) && !is.null(histo_df)){

    renamed_histo_df <- histo_df |>
      dplyr::mutate(tag = marker)

  }

  #Saving renamed histo
  if(!is.null(out))
    utils::write.table(renamed_histo_df, file = glue::glue("{out}"), quote = FALSE, row.names = FALSE, col.names = FALSE)

  #Return a marked histogram
  return(renamed_histo_df)
}
