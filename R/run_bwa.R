#' Runs bwa mem to map the given sequencing files into a reference genome.
#' If given the variable, also runs bwa mem and samtools to remove possible contaminated reads from a contaminant_ref reference file.
#' Also if given another variable, runs samtools to verify mapping quality and filter reads based on it.
#' Works with short reads from the Illumina sequencing platform.
#'
#' @param input Input FASTQ files/paths.
#' @param output Mapped reads in the given reference genome.
#' @param ref Reference genome where the reads will be mapped.
#' @param contaminant_ref Optional: Reference contaminant_ref genome, where reads that have matched with will be removed.
#' @param mapping_quality Optional: Value from 0 to 60, to filter reads with mapping quality equal or above the given value.
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
#' #2. Run run_bwa function
#' run_bwa(input = temp_fastq,
#'  output = "test_mapped",
#'  ref = temp_ref_fasta,
#'  keep_index = FALSE)
#' #3. Unlink ouput file
#' unlink("test_mapped.bam")
#'
run_bwa <- function(input,
                    output,
                    ref,
                    contaminant_ref = NULL,
                    mapping_quality = NULL,
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

  # Remove contaminant_ref ------------------------------------------------------
  #If contaminant_ref file was given
  if(!is.null(contaminant_ref)){
    #Create list of index files
    index_files_contaminant_ref <- paste0(contaminant_ref, ".", extensions)
    #Check if file exists in the given path
    if(!file.exists(contaminant_ref)){
      stop("ERROR: Contaminant reference file not found in the given path.")
    }
      #Index if no index files were found
      if(!all(file.exists(index_files_contaminant_ref))){
        message("No indexed files found, indexing contaminant reference file...")
        system(glue::glue("bwa index {contaminant_ref}"))
      } else{
        message("Contaminant reference indexed files found, proceeding to read mapping...")
      }

      #If keep index is false remove the index files
      if(isFALSE(keep_index)){
        on.exit(unlink(glue::glue("{contaminant_ref}.{c('amb', 'ann', 'bwt', 'pac', 'sa')}")), add = TRUE)
      }

      #Mapping to remove contaminant_ref
      message("Mapping reads to contaminant_ref reference file...")
      system(glue::glue("bwa mem -t 10 -v 2 -M {contaminant_ref} {input_col} > {input[1]}_mapped_cont.sam"))

      #Select what was not mapped
      message("Selecting reads that were not mapped in contaminant_ref reference file...")
      #Check if it was given one or two files
      if(length(input) == 1){
        #Samtools with the -4 filter to get the single end reads
        system(glue::glue("samtools view -b -f 4 {input[1]}_mapped_cont.sam > {input}_removed_cont.bam"))
        #Remove sam output
        unlink(glue::glue("{input[1]}_mapped_cont.sam"))

        #Sort to generate fastq with pairs ordered
        system(glue::glue("samtools sort -n {input}_removed_cont.bam > {input}_removed_cont_sorted.bam"))

        message("Creating fastaq file from with reads that were not mapped in contaminant...")
        system(glue::glue("samtools fastq -o {input}_rm_contaminant_ref.fastq.gz {input}_removed_cont_sorted.bam"))
        input_rm_cont <- glue::glue("{input}_rm_contaminant_ref.fastq.gz")

        on.exit(unlink(glue::glue("{input}_rm_contaminant_ref.fastq.gz")), add = TRUE)
      }

      else if(length(input) == 2){
        #Samtools with -f 12 to get the paired reads
        system(glue::glue("samtools view -b -f 12 {input[1]}_mapped_cont.sam > {input[1]}_removed_cont.bam"))
        #Remove sam output
        unlink(glue::glue("{input[1]}_mapped_cont.sam"))

        #Sort to generate fastq with pairs ordered
        system(glue::glue("samtools sort -n {input[1]}_removed_cont.bam > {input[1]}_removed_cont_sorted.bam"))

        message("Creating fastaq file from with reads that were not mapped in contaminant...")
        system(glue::glue("samtools fastq -1 {input[1]}_rm_contaminant_ref.fastq.gz -2 {input[2]}_rm_contaminant_ref.fastq.gz {input[1]}_removed_cont_sorted.bam"))
        input_rm_cont <- glue::glue("{input[1]}_rm_contaminant_ref.fastq.gz {input[2]}_rm_contaminant_ref.fastq.gz")

        on.exit(unlink(glue::glue("{input[1]}_rm_contaminant_ref.fastq.gz")), add = TRUE)
        on.exit(unlink(glue::glue("{input[2]}_rm_contaminant_ref.fastq.gz")), add = TRUE)

      }else{
        stop("ERROR: The input should be one or two files")
      }

      #Remove temporary files
      on.exit(unlink(glue::glue("{input[1]}_removed_cont.bam")), add = TRUE)
      on.exit(unlink(glue::glue("{input[1]}_removed_cont_sorted.bam")), add = TRUE)
  }


  # Mapping -----------------------------------------------------------------
  #Index if any index file was not found
  if(!all(file.exists(index_files_ref))){
    message("No indexed files found, indexing reference file...")
    system(glue::glue("bwa index {ref}"))
  } else{
    message("Reference indexed files found, proceeding to read mapping...")
  }
  #Remove index files if keep index is false
  if(isFALSE(keep_index)){
    on.exit(unlink(glue::glue("{ref}.{c('amb', 'ann', 'bwt', 'pac', 'sa')}")), add = TRUE)
  }

  #Mapping with bwa mem
  message("Mapping reads to reference file...")
  if(!is.null(contaminant_ref)){ #If removed contaminant_ref run from files with contaminant_ref removed
    system(glue::glue("bwa mem -M {ref} {input_rm_cont} -t 10 > {input[1]}_mapped.sam"))
    unlink(glue::glue("{input}_rm_contaminant_ref.fastq.gz"))
    unlink(glue::glue("{input[1]}_rm_contaminant_ref.fastq.gz"))
    unlink(glue::glue("{input[2]}_rm_contaminant_ref.fastq.gz"))

  }else{ #If contaminant_refs were not removed run from given files
    system(glue::glue("bwa mem -v 2 -M {ref} {input_col} -t 10 > {input[1]}_mapped.sam"))
  }

  #If filter mapping quality, filter with the given mapping quality
  if(!is.null(mapping_quality)){
    #Check if variable is a number
    if(!is.numeric(mapping_quality)){
      stop("ERROR: The mapping quality should be a numeric variable, between 0 and 60")
    }else{
      if(mapping_quality < 0 || mapping_quality > 60){
        stop("ERROR: The mapping quality value should be between 0 and 60")
      }
    }
    #Run samtools view to filter mapping quality
    message(glue::glue("Sorting and filtering mapping quality with quality of {mapping_quality}..."))
    #If wanted output is sam save the output as sam file
    system(glue::glue("samtools view -h -q {mapping_quality} {input[1]}_mapped.sam | samtools sort -o {output}.{output_format}"))

    }else{
      #Sorting
      message("Sorting and ranaming output file...")
      system(glue::glue("samtools sort {input[1]}_mapped.sam -o {output}.{output_format}"))
      }

  on.exit(unlink(glue::glue("{input[1]}_mapped.sam")), add = TRUE)
  message("Removing temporary files...")
}
