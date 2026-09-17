#' Runs trimmomatic to verify read quality using the quality score of phred 33, standard for Illumina data, trim and separate paired and unpaired reads.
#' Works with short reads from the Illumina sequencing platform.
#'
#' @param mode Which mode to run.
#' \describe{
#'  \item{SE}{Single-end data, requires only one input FASTQ.}
#'  \item{PE}{Paired-end data, requires two input FASTQ, one containing the forward and reverse sequences.}
#'  }
#' @param input Input FASTQ files/paths.
#' \describe{
#'  \item{1 file}{For single-end mode, give only one FASTQ file.}
#'  \item{2 files}{If paired-end mode, give the forward and reverse sequen FASTQ.}
#'  }
#' @param output Output FASTQ files names/paths(4 files if pair ended mode, 1 file for single ended mode)
#' \describe{
#'  \item{1 output file}{If single-end mode, output contains only the reads that met the given phred/window_size requirements.}
#'    \item{4 output files}{If paired-end mode, saving what was paired and unpaired, in this order, between forward and reverse sequences and met the given phred/window_size requirements.}
#'  }
#' @param minimun_read_length Minimun read length to kept, default is 36.
#' @param sliding_window_size The sliding window size to verify average phred value, default is 4.
#' @param average_phred The average phred value for the sliding window, default is 25, if the phred value drops below it the read is cut.
#' @param adapter_path The path for the adapter used in sequencing to be removed from reads.
#' \describe{Default is TruSeq3(Ideal for Hiseq/Miseq), Check trimmomatic documentation for more options}

#' @importFrom rlang .data
#' @export
#'
#' @examplesIf Sys.which("trimmomatic") != ""
#' #1. Load .fastq test file
#' temp_fastq <- system.file("extdata", "test.fastq", package = "AEkcompgen")
#' #2. Run run_trimmomatic function
#' run_trimmomatic(mode = "SE", input = temp_fastq, output = "test_trimmed.fastq")
#' #3. Unlink ouput file
#' unlink("test_trimmed.fastq")
#'
run_trimmomatic <- function(mode,
                            input,
                            output,
                            minimun_read_length = 36,
                            sliding_window_size = 4,
                            average_phred = 25,
                            adapter_path = NULL){
  #Verify if trimmomatic is intalled and in users path
  if (Sys.which("trimmomatic") == "") {
    stop(
      "ERROR: The program 'trimmomatic' was not found in the system\n
      Please check if it's installed and add it to your system's path or environment\n
      Recommendation: Install it via Conda",
      call. = FALSE
    )
  }

  #Collapse input and output values into a single line
  input_col <- paste(input, collapse = " ")
  output_col <- paste(output, collapse = " ")




  #Check if quantity of files is compatible with chosen mode
  if(mode == "PE"){
    #Check if the quantity of files/file names is compatible with the chosen mode
    if(length(input) != 2 || length(output) != 4){
      stop("ERROR: Pair ended mode requires 2 input files, and 4 ouput file names.",
           call. = FALSE)

    }
    #Check if file exists in the given path
    if(!all(file.exists(input))) {
        stop("ERROR: One or more Input .fastq files not found in the given path.")
      }

    message("Running Trimmomatic in paired-end mode...")
    #Put adapters to remove
    if(is.null(adapter_path)){
        adapter_path = "TruSeq3-PE.fa"}

    }

  else if(mode == "SE"){
    #Check if the quantity of files/file names is compatible with the chosen mode
    if(length(input) != 1 || length(output) != 1){
      stop("Single ended mode requires 1 input file and 1 output file",
               call = FALSE)
      }
    #Check if file is available in the given path
    if(!file.exists(input)){
      stop("ERROR: Input .fastq file not found in the given path.")
      }

    message("Running Trimmomatic in single-end mode...")

    #Put adapters to remove
    if(is.null(adapter_path)){
      adapter_path = "TruSeq3-SE.fa"}

    } else{
      stop(
        "ERROR: The paramater for mode given is invalid.\n
        Please chose between 'PE' for paired-end mode or 'SE' for single end mode"
      )
    }
  #Run trimmomatic
  system(glue::glue("trimmomatic {mode} -phred33 {input_col} {output_col} ILLUMINACLIP:{adapter_path}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:{sliding_window_size}:{average_phred} MINLEN:{minimun_read_length}"))


}
