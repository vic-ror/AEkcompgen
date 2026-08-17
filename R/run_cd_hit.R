#' CD-HIT-est for the given dataframe or fasta_file
#'
#' @param data_frame A dataframe containing id, seq and label columns.
#' @param fasta_file A fasta file or a path to one.
#' @param identity The identity degree of similarity that will be used to cluster the sequences.
#' @param out If you want to save the output file, the give the name for it.
#'
#' @return A dataframe containing the clustered groups
#' @importFrom rlang .data
#' @export
#'
#' @examplesIf Sys.which("cd-hit") != ""
#' # 1. Create temporary dataframe
#' data_frame_test <- data.frame(id = c("seq1", "seq2", "seq3"),
#' seq = c("GACAGGTACAAGAAGGAGTA", "AGGGCGACCTTCGATTCGGA", "TTTACACACTCTCCTTGGAC")
#' )
#' # 2. Run function with temporary file
#' cd_hit_test <- run_cd_hit(data_frame = data_frame_test, identity = 0.9)
#' # 3. View resulting line dataframe
#' print(cd_hit_test)
#'


run_cd_hit <- function(data_frame = NULL, fasta_file = NULL, identity, out = NULL) {
  #Check if cd-hit is installed
  if (Sys.which("cd-hit") == "") {
    stop(
      "ERROR: The program'cd-hit' was not found in the system \n
      Please check if it's intalled and add it to your system's path",
      call. = FALSE
    )
  }
  #Check if identity is a numerical number and ir it is of the correct value
  if (!is.numeric(identity) || length(identity) != 1){
    stop("ERROR: Identity must be a numeric value between 0.75 and 1.")
  }

  if (identity < 0.75 || identity > 1){
    stop("ERROR: Identity must be between 0.75 and 1.")
  }

  if (is.null(data_frame) && !is.null(fasta_file)) {
    fasta_head <- readLines(fasta_file, n = 2)

    kmer_length <- stringr::str_length(fasta_head[2])

    if (kmer_length < 11){
      stop("CD-HIT est can only work with sequences greater than 10 nucleotides")
    }
  }

  else if (!is.null(data_frame) && is.null(fasta_file)){
    #If only a dataframe is provided create fasta file to run cd-hit
    #If a dataframe from the shared_seqs function is given
    if("id_2" %in% colnames(data_frame)){
      message("Shared seqs dataframe given, merging ids...")
      data_frame <- data_frame |>
        tidyr::unite("id", dplyr::starts_with("id_"), sep = "|", remove = FALSE)
    }

    message("Creating fasta file to run cd-hit...")
    temp_fasta <- data_frame |> dplyr::mutate(id = stringr::str_c(">", .data$id)) |>
      tidyr::pivot_longer(cols = c(.data$id, .data$seq), names_to = "col_name") |>
      dplyr::select(.data$value)

    #Check if sequence size is compatible (greater than 10 nucleotides)
    kmer_length <- stringr::str_length(data_frame$seq[1])

    if (kmer_length < 11){
      stop("CD-HIT est can only work with sequences greater than 10 nucleotides")
    }

    utils::write.table(temp_fasta, file = "temp.fasta", quote = FALSE, row.names = FALSE, col.names = FALSE)

    fasta_file = "temp.fasta"

    on.exit(
      if(file.exists("temp.fasta")) unlink("temp.fasta"), add = TRUE
    )

  }

  #If no dataframe or fasta file was provided
  else if (is.null(data_frame) && is.null(fasta_file)){
    stop("A dataframe containing sequence id and sequence OR a fasta file path must be provided")
  }

  #If both files are provided
  else if (!is.null(data_frame) && !is.null(fasta_file)){
    stop("Please provide EITHER a dataframe OR a fasta file, not both.")
  }


  #Running CD-HIT-est
  ##Necessary to change the word-size depending on the identity given
  message(glue::glue("Clustering sequences with an identity of {identity}..."))
  if(identity >= 0.95 && identity <= 1){
    system(glue::glue("cd-hit-est -i {fasta_file} -o temp_file_cd_hit -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identity} -n 10 -b 1"))
  }
  else if(identity >= 0.9 && identity < 0.95){
    system(glue::glue("cd-hit-est -i {fasta_file} -o temp_file_cd_hit -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identity} -n 8 -b 1"))
  }
  else if(identity >= 0.88 && identity < 0.9){
    system(glue::glue("cd-hit-est -i {fasta_file} -o temp_file_cd_hit -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identity} -n 7 -b 1"))
  }
  else if(identity >= 0.85 && identity < 0.88){
    system(glue::glue("cd-hit-est -i {fasta_file} -o temp_file_cd_hit -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identity} -n 6 -b 1"))
  }
  else if(identity >= 0.8 && identity < 0.85){
    system(glue::glue("cd-hit-est -i {fasta_file} -o temp_file_cd_hit -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identity} -n 5 -b 1"))
  }
  else if(identity >= 0.75 && identity < 0.8){
    system(glue::glue("cd-hit-est -i {fasta_file} -o temp_file_cd_hit -d 0 -T 16 -g 0 -M 75000 -aL 0.97 -aS 0.97 -c {identity} -n 4 -b 1"))
  }

  on.exit(
    if(file.exists("temp_file_cd_hit.clstr")) unlink("temp_file_cd_hit.clstr"), add = TRUE
  )

  #Check if CD-HIT-est worked
  if(!file.exists("temp_file_cd_hit.clstr")){
    stop("ERROR: CD-HIT-est failed. No .clstr file was generated.")
  }

  #Turn clster output file into a R dataframe
  message("Turning .clstr file into a dataframe...")

  #Load .clstr file
  clstr_file <- readr::read_lines("temp_file_cd_hit.clstr", lazy = FALSE)

  if (length(clstr_file) == 0) {
    stop("Error: The clstr file is empty.")
  }

  #Identify cluster and kmers
  clstr_df <- tibble::tibble(clstr_file) |>
    dplyr::mutate(is_cluster_id = stringr::str_detect(.data$clstr_file, "^>")) |>
    dplyr::mutate(group_id = cumsum(.data$is_cluster_id))

  cluster <- clstr_df |>
    dplyr::filter(.data$is_cluster_id == TRUE) |>
    dplyr::rename(cluster = .data$clstr_file) |>
    dplyr::select(.data$cluster, .data$group_id) |>
    dplyr::mutate(cluster = stringr::str_remove(.data$cluster, ">"))

  kmer <- clstr_df |>
    dplyr::filter(.data$is_cluster_id == FALSE) |>
    dplyr::rename(id = .data$clstr_file) |>
    dplyr::mutate(id = stringr::str_remove(.data$id, ".*>")) |>
    dplyr::mutate(id = stringr::str_remove(.data$id, "\\.\\.\\..*")) |>
    dplyr::select(.data$id, .data$group_id)

  joined_info_clst <- dplyr::full_join(cluster, kmer, by = "group_id", multiple = "all") |>
    dplyr::select(.data$cluster, .data$id)

  #If a dataframe from the shared_seqs function is given
  if("id_2" %in% colnames(data_frame)){
    id_cols <- data_frame |>
      dplyr::select(dplyr::starts_with("id")) |>
      colnames()

      joined_info_clst <- joined_info_clst |>
        tidyr::separate_wider_delim(
          cols = .data$id,
          delim = "|",
          names = id_cols
        )
  }

  #Remove the temporary files
  #If an output is given save the .clstr file from CD-HIT
  if(!is.null(out)){
    message("Saving output .clstr file...")
    file.rename(from = "temp_file_cd_hit.clstr", to = glue::glue("{out}.clstr"))
  }

  #Remove other files
  message("Removing temporary files...")
  if(file.exists("temp_file_cd_hit")) unlink("temp_file_cd_hit")

  return(joined_info_clst)

}


