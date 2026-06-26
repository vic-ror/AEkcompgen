#' CD-HIT-est for the given dataframe or fasta_file
#'
#' @param data_frame A dataframe containing id, seq and marker columns.
#' @param fasta_file A fasta file or a path to one.
#' @param identity The identity degree of similarity that will be used to cluster the sequences.
#' @param output If you want to save the output file, the give the name for it.
#'
#' @return A dataframe containing the counted kmers and their headers with the frequency information
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # 1. Create temporary dataframe
#' data_frame_test <- data.frame(id = c("seq1", "seq2", "seq3"),
#' seq = c("GACAGGTACAAGAAGGAGTA", "AGGGCGACCTTCGATTCGGA", "TTTACACACTCTCCTTGGAC")
#' )
#' # 2. Run function with temporary file
#' cd_hit_test <- run_cd_hit(data_frame = data_frame_test, identity = 0.9)
#' # 3. View resulting line dataframe
#' print(cd_hit_test)
#'


run_cd_hit <- function(data_frame = NULL, fasta_file = NULL, identity, output = NULL) {
  #Check if cd-hit is installed
  if (Sys.which("cd-hit") == "") {
    stop(
      "ERROR: The program'cd-hit' was not found in the system \n
      Please check if it's intalled and add it to your system's path",
      call. = FALSE
    )
  }

  if (!is.null(data_frame) && is.null(fasta_file)){
    #If only a dataframe is provided create fasta file to run cd-hit
    message("Creating fasta file to run cd-hit...")
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
  else if(identity < 0.75){
    stop("CD-HIT-est is not capable of clustering sequences with an identity smaller than 0.75")
  }

  #Turn clster output file into a R dataframe
  message("Turning .clstr file into a dataframe...")
  read_clstr <- function(clstr_file){
    lines <- readLines(clstr_file)  #Read clstr file
    current_cluster <- NULL         #Create a variable to store cluster number
    current_seq <- NULL             #Create a variable to store current seq id
    cluster_list <- c()             #Create list to store clusters
    seq_list <- c()                 #Create list to store seq_ids

    #Loop for lines in the file
    for(line in lines){
      if (startsWith(line, ">")){                          #If starts with > its a cluster identification
        current_cluster <- sub(">Cluster", "", line)       #Remove >cluster from the string
      }
      else{                                                #If it doesnt start with > it is the seq id
        current_seq <- sub(".*>", "", line)                #Remove everything before id
        current_seq <- sub("\\.\\.\\..*", "", current_seq) #Remove everything after the id
        seq_list <- c(seq_list, current_seq)               #Add seq to list
        cluster_list <- c(cluster_list, current_cluster)   #Add cluster to list
      }
    }
    #Save lists in dataframe form
    df_clusters <- data.frame(seq_id = seq_list,
                              cluster = cluster_list,
                              stringsAsFactors = FALSE)
    return(df_clusters)
  }
  #Run function that will read the file in the output from cdhit
  output_df <- read_clstr("temp_file_cd_hit.clstr")

  #Remove the temporary files
  #If an output is given save the .clstr file from CD-HIT
  if(!is.null(output)){
    message("Saving output .clstr file...")
    file.rename(from = "temp_file_cd_hit.clstr", to = glue::glue("{output}.clstr"))
  }
  #If it isnt true remove it
  else {
    if(file.exists("temp_file_cd_hit.clstr")) unlink("temp_file_cd_hit.clstr")
  }

  #Remove other files
  message("Removing temporary files...")
  if(file.exists("temp_file_cd_hit")) unlink("temp_file_cd_hit")
  if(file.exists("temp.fasta")) unlink("temp.fasta")

  return(output_df)

  }


