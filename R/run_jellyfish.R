#' Runs jellyfish and give it's output, both in dataframe form and saving it in the directory given in the output argument
#' The id for the sequence contains the row of the kmer in the fasta file and the the count of the k-mer
#'
#' @param data_frame A dataframe containing the ID and sequence collumns
#' @param fasta_file A fasta file or the path to one.
#' @param length The length of kmers that will be counted
#' @param marker The marker that will be used in the header of sequence and file names
#' @param lower_count Don't output k-mers with count lower than this
#' @param upper_count Don't output k-mers with count higher than this
#' @param out Name for output fasta file
#' @return A dataframe containing the counted kmers and their headers with the frequency information
#' @importFrom rlang .data
#' @export
#'
#' @examplesIf Sys.which("jellyfish") != ""
#' # 1. Create temporary fasta file
#' data_frame_test <- data.frame(id = c("seq1", "seq2", "seq3"),
#' seq = c("GACAGGTACAAGAAGGAGTA",
#' "AGGGCGACCTTCGATTCGGA", 
#' "TTTACACACTCTCCTTGGAC"))
#' # 2. Run function with temporary file
#' line_file <- run_jellyfish(data_frame = data_frame_test,
#'  length =  10,
#'  marker = "test",
#'  out = "test_output")
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
                          marker,
                          lower_count = NULL,
                          upper_count = NULL,
                          out){
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
    temp_fasta <- data_frame |> dplyr::mutate(id = stringr::str_c(">", id)) |>
      tidyr::pivot_longer(cols = c(id, seq), names_to = "col_name") |>
      dplyr::select(value)

    write.table(temp_fasta, file = "temp.fasta", quote = FALSE, row.names = FALSE, col.names = FALSE)


    fasta_file = "temp.fasta"
  }

  #If no dataframe or fasta file was provided
  else if (is.null(data_frame) && is.null(fasta_file)){
    stop("A dataframe containing sequence id and sequence OR a fasta file path must be provided")
  }

  #If both files are provided
  else if (!is.null(data_frame) && !is.null(fasta_file)){
    stop("Please provide EITHER a dataframe OR a fasta file, not both.")
  }

  #Get the output path for the output file
  out_dir <- dirname(out)

  #Run_jellyfish
  #With lower count filter = Prints only kmers with count higher or equal to the filter
  if (!is.null(lower_count) && is.null(upper_count)){
    message(glue::glue("Counting k-mers with count higher or equal to {lower_count}..."))
    system(glue::glue("jellyfish count -m {length} -L {lower_count} -s 275M -t 10 -C -o {out_dir}/{marker}.jf {fasta_file}"))
  }
  #With higher count filter = Prints only kmers with count smaller or equal to the filter
  else if (is.null(lower_count) && !is.null(upper_count)){
    message(glue::glue("Counting k-mers with count smaller or equal to {upper_count}..."))
    system(glue::glue("jellyfish count -m {length} -U {upper_count} -s 275M -t 10 -C -o {out_dir}/{marker}.jf {fasta_file}"))
  }
  #With higher and lower count filter = Prints only kmers with count higher, smaller or equal to the filter
  else if (!is.null(lower_count) && !is.null(upper_count)){
    message(glue::glue("Counting k-mers with count smaller or equal to {upper_count} and higher or equal to {lower_count}..."))
    system(glue::glue("jellyfish count -m {length} -U {upper_count} -L {lower_count} -s 275M -t 10 -C -o {out_dir}/{marker}.jf {fasta_file}"))
  }
  #Without count filter
  else {
    message("Counting k-mers without any filters...")
    system(glue::glue("jellyfish count -m {length} -s 275M -t 10 -C -o {out_dir}/{marker}.jf {fasta_file}"))
  }

  message("Creating readable fasta file...")
  system(glue::glue("jellyfish dump {out_dir}/{marker}.jf > temp_jellyfish_dump"))# turn output into fasta

  message("Renaming fasta file headers to contain the row number and the count of the k-mer...")

  #Function to rename the header
  rename_kmer_header <- function(jf_readable_file){
    file <- readLines(jf_readable_file) #Read the file

    #Create arguments
    current_kmer_id = NULL # Create a string to store the current k-mer id
    line_number = 0       #Create a argument to store the line number
    kmer_id_list = c()    #Create a list to store k-mer ids
    kmer_seq_list = c()   #Create a list to store k-mer sequences


    for(line in file) {
      if (startsWith(line, ">")){           #If it is the header
        line_number = line_number + 1       #Add one to kmer_number
        current_kmer_id <- sub(">", glue::glue(">{marker}_{line_number}_"), line)  #Store information of marker and line number to the header
      }
      else{
        kmer_id_list = c(kmer_id_list, current_kmer_id)       #Add kmer id to the list
        kmer_seq_list = c(kmer_seq_list, line)  #If it is a sequence add it to the list of sequences
      }
    }
    #Save lists in a dataframe
    renamed_jf_dump_dataframe <- data.frame(id = kmer_id_list,
                                            seq = kmer_seq_list,
                                            stringsAsFactors = FALSE)
    return(renamed_jf_dump_dataframe)
  }

  message("Loading fasta in dataframe...")
  renamed_jf_df <- rename_kmer_header("temp_jellyfish_dump")

  jf_dataframe <- renamed_jf_df |>
    dplyr::mutate(id = stringr::str_remove(.data$id, ">"))

  message("Removing intermediate files...")
  #Remove intermediate files
  if(file.exists("temp_jellyfish_dump")) unlink("temp_jellyfish_dump")
  if(file.exists("temp.fasta")) unlink("temp.fasta")

  #Save renamed fasta file
  message("Saving renamed fasta file...")
  renamed_jf_fasta <- renamed_jf_df |>
    tidyr::pivot_longer(cols = c(id, seq), names_to = "col_name") |>
    dplyr::select(value)

  utils::write.table(renamed_jf_fasta, file = glue::glue("{out}.fasta"), quote = FALSE, row.names = FALSE, col.names = FALSE)


  #Print dataframe content
  return(jf_dataframe)

}
