#' Runs bwa mem to remove contaminant based on what was mapped in the contaminant reference gennome.
#' Works with short reads from the Illumina sequencing platform.
#'
#' @param input Input FASTQ files/paths.
#' @param output Mapped reads in the given reference genome.
#' @param ref Reference contaminant genome where the reads will be mapped.
#' @param output_format Format for the output of the mapping file (bam/sam), default is "bam' for binary file.
#' @param keep_index Optional: Default is true, put false if you dont want to keep the index files in the directory for faster performance.
#'
#' @importFrom rlang .data
#' @export
#'
#' @examplesIf Sys.which("bwa") != "" && Sys.which("samtools") != ""
#' #1. Load .fastq test file
#' temp_fastq <- system.file("extdata", "test.fastq", package = "AEkcompgen")
#' temp_ref_fasta <- system.file("extdata", "test_reference.fasta", package = "AEkcompgen")
#' #2. Run run_bwa_remove_contaminant function
#' run_bwa_remove_contaminant(input = temp_fastq,
#'  output = "test_mapped",
#'  ref = temp_ref_fasta,
#'  keep_index = FALSE)
#' #3. Unlink ouput file
#' unlink("test_mapped.bam")
#'
run_bwa_remove_contaminant <- function(input,
                    output,
                    ref,
                    output_format = "bam",
                    keep_index = TRUE
                    ){
  #Check if needed programs are installed
  #Verify if bwa is intalled and in users path
  if (Sys.which("bwa") == "") {
    stop(
      "ERROR: The program 'bwa' was not found in the system\n
      Please check if it's installed and add it to your system's path or environment\n
      Recommendation: Install it via Conda",
      call. = FALSE
    )
  }

  if (Sys.which("samtools") == "") {
    stop(
      "ERROR: The program 'samtools' was not found in the system\n
      Please check if it's installed and add it to your system's path or environment\n
      Recommendation: Install it via Conda",
      call. = FALSE
    )
  }
  #Check if files are available in the path given.
  #Check if file exists in the given path
  if(!all(file.exists(c(input, ref)))){
    stop("ERROR: Input or reference file were not found in the given path.")
  }

  input_col <- paste(input, collapse = " ")

  #Check if the output_format is bam or sam
  if(output_format != "sam"  && output_format !="bam"){
    stop("ERROR: Output file extension given must be 'sam' or 'bam'")
  }

  #List of index extensions
  extensions <- c("amb", "ann", "bwt", "pac", "sa")
  #Create list of index files
  index_files_ref <- paste0(ref, ".", extensions)

  #Remove contaminant
  #Index if no index files were found
  if(!all(file.exists(index_files_ref))){
    message("No indexed files found, indexing contaminant reference file...")

    system(glue::glue("bwa index {ref}"))
  } else{
    message("Contaminant reference indexed files found, proceeding to read mapping...")
  }

  #If keep index is false remove the index files
  if(isFALSE(keep_index)){
    on.exit(unlink(glue::glue("{ref}.{c('amb', 'ann', 'bwt', 'pac', 'sa')}")), add = TRUE)
  }

  #Mapping to remove ref
  message("Mapping reads to contaminant reference file...")
  system(glue::glue("bwa mem -t 10 -v 2 -M {ref} {input_col} > {input[1]}_mapped_cont.sam"))

  #Select what was not mapped
  message("Selecting reads that were not mapped in ref reference file...")
  #Check if it was given one or two files
  if(length(input) == 1){
    #Samtools with the -4 filter to get the single end reads
    system(glue::glue("samtools view -b -f 4 {input[1]}_mapped_cont.sam > {input}_removed_cont.bam"))
    #Remove sam output
    unlink(glue::glue("{input[1]}_mapped_cont.sam"))

    #Sort to generate fastq with pairs ordered
    system(glue::glue("samtools sort -n {input}_removed_cont.bam > {input}_removed_cont_sorted.bam"))

    message("Creating fastaq file from with reads that were not mapped in contaminant...")
    system(glue::glue("samtools fastq -o {input}_rm_ref.fastq.gz {input}_removed_cont_sorted.bam"))
    input_rm_cont <- glue::glue("{input}_rm_ref.fastq.gz")

    on.exit(unlink(glue::glue("{input}_rm_ref.fastq.gz")), add = TRUE)
  }

  else if(length(input) == 2){
    #Samtools with -f 12 to get the paired reads
    system(glue::glue("samtools view -b -f 12 {input[1]}_mapped_cont.sam > {input[1]}_removed_cont.bam"))
    #Remove sam output
    unlink(glue::glue("{input[1]}_mapped_cont.sam"))

    #Sort to generate fastq with pairs ordered
    system(glue::glue("samtools sort -n {input[1]}_removed_cont.bam > {input[1]}_removed_cont_sorted.bam"))

    message("Creating fastaq file from with reads that were not mapped in contaminant...")
    system(glue::glue("samtools fastq -1 {input[1]}_rm_ref.fastq.gz -2 {input[2]}_rm_ref.fastq.gz {input[1]}_removed_cont_sorted.bam"))
    input_rm_cont <- glue::glue("{input[1]}_rm_ref.fastq.gz {input[2]}_rm_ref.fastq.gz")

    on.exit(unlink(glue::glue("{input[1]}_rm_ref.fastq.gz")), add = TRUE)
    on.exit(unlink(glue::glue("{input[2]}_rm_ref.fastq.gz")), add = TRUE)

  }else{
    stop("ERROR: The input should be one or two files")
  }

  #Remove temporary files
  on.exit(unlink(glue::glue("{input[1]}_removed_cont.bam")), add = TRUE)
  on.exit(unlink(glue::glue("{input[1]}_removed_cont_sorted.bam")), add = TRUE)
}
