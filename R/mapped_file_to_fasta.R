#' Turns bam/sam file into a fasta file.
#'
#' @param input A sorted bam/sam file, (run_bwa already sorts its output)
#' @param output The path to save the .duckdb file, if none given it creates a temporary file.
#'
#' @return Saves in the current directory a fasta file containing the reads information from sam/bam file.
#' @export
#' @importFrom rlang .data
#'
#' @examplesIf Sys.which("samtools") != ""
#' #1. Load .bam test file
#' temp_bam <- system.file("extdata", "test_mapped.bam", package = "AEkcompgen")
#' #2. Run mapped_file_to_fasta function
#' mapped_file_to_fasta(input =  temp_bam, output = "temp_bam")
#' #3. Remove temporary file
#' unlink("temp_bam.fasta")
mapped_file_to_fasta <- function(input,
                         output){
  #Check if samtools is installed and in the users  path
  if (Sys.which("samtools") == "") {
    stop(
      "ERROR: The program 'samtools' was not found in the system\n
      Please check if it's installed and add it to your system's path or environment\n
      Recommendation: Install it via Conda",
      call. = FALSE
    )
  }

  system(glue::glue("samtools sort -n {input}| samtools fasta > {output}.fasta"))

}
