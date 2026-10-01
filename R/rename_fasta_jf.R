#' Loads jellyfish dump fasta output into a dataframe, and renames the header by adding the dataset_label and the rownumber
#'
#' @param fasta_file A fasta output from jellyfish dump contaning the kmer count in the header
#' @param dataset_label A label to be added to the header of the fasta, identifying it
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
#' line_file <- rename_fasta_jf(fasta_file = fasta_temp, dataset_label = "tag")
#' # 3. View resulting line dataframe
#' print(line_file)
#' # 4. Delete temporary file
#'unlink(fasta_temp)

rename_fasta_jf <- function(fasta_file, dataset_label, out = NULL){

  #Load fasta file
  #Read it with readr
   fasta <- readr::read_lines(fasta_file, lazy = FALSE)

   #Check if given file is empty
  if (length(fasta) == 0) {
    stop("Error: The given fasta file is empty.")
  }

   #Check if given dataset label has a _ in it
   #_ would give an error later on
   if(isTRUE(stringr::str_detect(dataset_label, "_"))){
     stop("ERROR: dataset_label must not contain '_', since it will create problems in later functions, please rename it.")
   }

  #Remove the > and identify elements in the header/sequence
  fasta_df <- tibble::tibble(fasta) |>
    dplyr::mutate(is_header = stringr::str_detect(.data$fasta, "^>")) |>
    dplyr::mutate(group_id = cumsum(.data$is_header)) #identify the sequence with the header

  header <- fasta_df |> dplyr::filter(stringr::str_detect(.data$fasta, ">")) |>
    dplyr::rename(id = .data$fasta) |>
    dplyr::select(.data$id, .data$group_id)

  seq <- fasta_df |> dplyr::filter(!stringr::str_detect(.data$fasta, ">")) |>
    dplyr::rename(seq = .data$fasta) |>
    dplyr::select(.data$seq, .data$group_id)

  #Rename header to contain elements
  joined_fasta <- dplyr::full_join(header, seq, by = "group_id", multiple = "all") |>
    dplyr::mutate(id = stringr::str_remove(.data$id, ">")) |>
    dplyr::mutate(nrow = dplyr::row_number()) |>
    dplyr::mutate(id = stringr::str_glue("{dataset_label}_{nrow}_{id}")) |>
    dplyr::select(.data$id, .data$seq)

  #Saving renamed fasta
  if(!is.null(out))
    utils::write.table(joined_fasta, file = glue::glue("{out}"), quote = FALSE, row.names = FALSE, col.names = FALSE)

return(joined_fasta)

}
