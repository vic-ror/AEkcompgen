#' Converts fasta file to line
#'
#' @param fasta_file A fasta file or the path to one.
#'
#' @return A dataframe with the collumns "k-mer_id", "sequence", organized like a line file
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # 1. Create temporary fasta file
#' fasta_temp <- tempfile(fileext = ".fasta")
#' writeLines(
#' text = c(
#' ">seq1_1_2", "ATCGATCG",
#' ">seq2_3_4", "TCGATCGA"),
#' con = fasta_temp
#' )
#' # 2. Run function with temporary file
#' fasta <- load_fasta(fasta_temp)
#' # 3. View resulting line dataframe
#' print(fasta)
#' # 4. Delete temporary file
#' unlink(fasta_temp)
load_fasta <- function(fasta_file){

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
    dplyr::rename(id = .data$fasta) |>
    dplyr::select(.data$id, .data$group_id)

  seq <- fasta_df |> dplyr::filter(!stringr::str_detect(.data$fasta, ">")) |>
    dplyr::rename(seq = .data$fasta) |>
    dplyr::select(.data$seq, .data$group_id)

  joined_fasta <- dplyr::full_join(header, seq, by = "group_id", multiple = "all") |>
    dplyr::mutate(id = stringr::str_remove(.data$id, ">")) |>
    dplyr::select(.data$id, .data$seq)

  if(stringr::str_count(joined_fasta$id[1], "_") != 2){
    stop("ERROR: The dataset_label must not contain a '_', since it will create problems in later functions, please rename it.")
  }

  return(joined_fasta)
}
