#' Loads jellyfish dump fasta output into a dataframe, and renames the header by adding the marker and the rownumber
#'
#' @param fasta_file A fasta output from jellyfish dump contaning the kmer count in the header
#' @param marker A marker to be added to the header of the fasta, identifying it
#' @param out The name of a path to save the renamed fasta file
#'
#' @return A dataframe containing the counted kmers and their headers with the frequency information
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # 1. Create temporary fasta file
#' fasta_temp <- tempfile(fileext = ".fasta")
#' writeLines(
#' text = c(">10", "ATCGATCG",
#' ">13", "TCGATCGA"),
#' con = fasta_temp
#' )
#' # 2. Run function with temporary file
#' line_file <- rename_fasta_jf(fasta_file = fasta_temp, marker = "tag")
#' # 3. View resulting line dataframe
#' print(line_file)
#' # 4. Delete temporary file
#'unlink(fasta_temp)

rename_fasta_jf <- function(fasta_file, marker, out = NULL){

  #Load fasta file
  #Read it with readr
   fasta <- readr::read_lines(fasta_file, lazy = FALSE)

   #Check if given file is empty
  if (length(fasta) == 0) {
    stop("Error: The given fasta file is empty.")
  }

  fasta_df <- tibble::tibble(fasta) |>
    dplyr::mutate(is_header = stringr::str_detect(.data$fasta, "^>")) |>
    dplyr::mutate(group_id = cumsum(.data$is_header)) #identify the sequence with the header

  header <- fasta_df |> dplyr::filter(stringr::str_detect(.data$fasta, ">")) |>
    dplyr::rename(header = .data$fasta) |>
    dplyr::select(.data$header, .data$group_id)

  seq <- fasta_df |> dplyr::filter(!stringr::str_detect(.data$fasta, ">")) |>
    dplyr::rename(seq = .data$fasta) |>
    dplyr::select(.data$seq, .data$group_id)

  joined_fasta <- dplyr::full_join(header, seq, by = "group_id", multiple = "all") |>
    dplyr::mutate(header = stringr::str_remove(.data$header, ">")) |>
    dplyr::mutate(nrow = dplyr::row_number()) |>
    dplyr::mutate(header = stringr::str_glue("{marker}_{nrow}_{header}")) |>
    dplyr::select(.data$header, .data$seq)

  #Saving renamed fasta
  if(!is.null(out))
    utils::write.table(joined_fasta, file = glue::glue("{out}"), quote = FALSE, row.names = FALSE, col.names = FALSE)

return(joined_fasta)

}
