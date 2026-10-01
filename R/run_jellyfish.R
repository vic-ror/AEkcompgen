#' Runs jellyfish and give it's output, both in dataframe form and saving it in the directory given in the output argument
#' The id for the sequence contains the row of the kmer in the fasta file and the the count of the k-mer
#'
#' @param data_frame A dataframe containing the ID and sequence collumns
#' @param fasta_file A fasta file or the path to one.
#' @param length The length of kmers that will be counted
#' @param dataset_label The dataset_label that will be used in the header of sequence and file names
#' @param lower_count Don't output k-mers with count lower than this
#' @param upper_count Don't output k-mers with count higher than this
#' @param output Name for output fasta file
#' @return A dataframe containing the counted kmers and their headers with the frequency information
#' @importFrom rlang .data
#' @export
#'
#' @examplesIf Sys.which("jellyfish") != ""
#' # 1. Create temporary fasta file
#' data_frame_test <- data.frame(id = c("seq1", "seq2", "seq3"),
#' seq = c("GACAGGTACAAGAAGGAGTA", "AGGGCGACCTTCGATTCGGA", "TTTACACACTCTCCTTGGAC"))
#' # 2. Run function with temporary file
#' line_file <- run_jellyfish(data_frame = data_frame_test,
#'  length =  10,
#'  dataset_label = "test")
#' # 3. View resulting line dataframe
#' print(line_file)
#' # 4. Delete temporary file
#' unlink("test_output.fasta")
#'unlink("temp.fasta")
#'unlink("test.jf")
#'unlink("temp_jellyfish_dump")
run_jellyfish <- function(data_frame = NULL,
                          fasta_file = NULL,
                          length,
                          dataset_label,
                          lower_count = NULL,
                          upper_count = NULL,
                          output = NULL){
  #Check if jellyfish is installed
  if (Sys.which("jellyfish") == "") {
    stop(
      "ERROR: The program 'jellyfish' was not found in the system/n
      Please check if it's installed and add it to your system's path",
      call. = FALSE
    )
  }

  if (!is.null(data_frame) && is.null(fasta_file)){
    #If only a dataframe is provided create fasta file to run jellyfish
    message("Creating fasta file to run jellyfish...")
    temp_fasta <- data_frame |> dplyr::mutate(id = stringr::str_c(">", .data$id)) |>
      tidyr::pivot_longer(cols = c(.data$id, .data$seq), names_to = "col_name") |>
      dplyr::select(.data$value)

    utils::write.table(temp_fasta, file = glue::glue("{dataset_label}_temp.fasta"), quote = FALSE, row.names = FALSE, col.names = FALSE)

    on.exit(
      if(file.exists(glue::glue("{dataset_label}_temp.fasta"))) unlink(glue::glue("{dataset_label}_temp.fasta")), add = TRUE
    )


    fasta_file = glue::glue("{dataset_label}_temp.fasta")
  }

  #If no dataframe or fasta file was provided
  else if (is.null(data_frame) && is.null(fasta_file)){
    stop("A dataframe containing sequence id and sequence OR a fasta file path must be provided")
  }

  #If both files are provided
  else if (!is.null(data_frame) && !is.null(fasta_file)){
    stop("Please provide EITHER a dataframe OR a fasta file, not both.")
  }

  #Check if fasta_file is available
  #Check if file exists in the given path
  if(!file.exists(fasta_file)){
    stop("ERROR: Input fasta file not found in the given path.")
  }

  #Check if given dataset label has a _ in it
  #_ would give an error later on
  if(isTRUE(stringr::str_detect(dataset_label, "_"))){
    stop("ERROR: dataset_label must not contain '_', since it will create problems in later functions, please rename it.")
  }

  #Run_jellyfish
  #With lower count filter = Prints only kmers with count higher or equal to the filter
  if (!is.null(lower_count) && is.null(upper_count)){
    message(glue::glue("Counting k-mers with count higher or equal to {lower_count}..."))
    system(glue::glue("jellyfish count -m {length} -L {lower_count} -s 275M -t 10 -C -o {dataset_label}.jf {fasta_file}"))
  }
  #With higher count filter = Prints only kmers with count smaller or equal to the filter
  else if (is.null(lower_count) && !is.null(upper_count)){
    message(glue::glue("Counting k-mers with count smaller or equal to {upper_count}..."))
    system(glue::glue("jellyfish count -m {length} -U {upper_count} -s 275M -t 10 -C -o {dataset_label}.jf {fasta_file}"))
  }
  #With higher and lower count filter = Prints only kmers with count higher, smaller or equal to the filter
  else if (!is.null(lower_count) && !is.null(upper_count)){
    message(glue::glue("Counting k-mers with count smaller or equal to {upper_count} and higher or equal to {lower_count}..."))
    system(glue::glue("jellyfish count -m {length} -U {upper_count} -L {lower_count} -s 275M -t 10 -C -o {dataset_label}.jf {fasta_file}"))
  }
  #Without count filter
  else {
    message("Counting k-mers without any filters...")
    system(glue::glue("jellyfish count -m {length} -s 275M -t 10 -C -o {dataset_label}.jf {fasta_file}"))
  }

  on.exit(
    if(file.exists(glue::glue("{dataset_label}.jf"))) unlink(glue::glue("{dataset_label}.jf")), add = TRUE
  )


  message("Creating readable fasta file...")
  system(glue::glue("jellyfish dump {dataset_label}.jf > {dataset_label}_temp_jellyfish_dump"))# turn output into fasta

  on.exit(
    if(file.exists(glue::glue("{dataset_label}_temp_jellyfish_dump"))) unlink(glue::glue("{dataset_label}_temp_jellyfish_dump")), add = TRUE
  )

  message("Renaming fasta file headers to contain the row number and the count of the k-mer...")

  #Load temporary fasta file
  #Read it with readr
  fasta <- readr::read_lines(glue::glue("{dataset_label}_temp_jellyfish_dump"), lazy = FALSE)

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
    dplyr::mutate(nrow = dplyr::row_number()) |>
    dplyr::mutate(id = stringr::str_glue("{dataset_label}_{nrow}_{id}")) |>
    dplyr::select(.data$id, .data$seq)


  #If an output file is given save it into a .fasta file
  if(!is.null(output)){
    message("Saving renamed fasta file...")

    joined_fasta_file <- joined_fasta |>
      dplyr::mutate(id = glue::glue(">{id}")) |>
      tidyr::pivot_longer(cols = c(.data$id, .data$seq),
                          names_to = "type",
                          values_to = "fasta_line") |>
      dplyr::select(.data$fasta_line)

    utils::write.table(joined_fasta_file, file = glue::glue("{output}.fasta"), quote = FALSE, row.names = FALSE, col.names = FALSE)

    if(file.exists(glue::glue("{dataset_label}.jf"))) file.rename(from = glue::glue("{dataset_label}.jf"), to = glue::glue("{output}.jf"))

  }


  #Print dataframe content
  return(joined_fasta)

}
