#' Creates an Euler Diagram based on the kmer sequence distribuition among datasets.
#' The suggested maximum number of datasets is 4 but can hold it up to 7 datasets, if there are more markes it c
#' @param kmer_sharing_relation A li from the kmer_sharing_relation function
#' @param path_db The path to save the .duckdb file, if none given it creates a temporary file
#' @param dataset_label A list with the names of the datasets wished to be in the diagram, KEEP THE datasetS' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#' @param title_name The title for the plot.
#'
#' @return An Euler Diagram based on the kmer distribuition among datasets
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframes
#' data_1 <- data.frame(id = c("tag1_1_1", "tag1_2_2", "tag1_3_3", "tag1_4_4", "tag1_5_5"),
#'  seq = c("AAAG", "AGAG", "CTGC", "AAAA", "GAGA"))
#' data_2 <- data.frame(id = c("tag2_1_1", "tag2_2_2", "tag2_3_3", "tag2_4_4", "tag2_5_5"),
#'  seq = c("AATG", "AGAG", "CAGC", "AATA", "GAGA"))
#' # 2. Run plot_venn_kmers function to get list
#' dataset_sharing_relation <- kmer_sharing_relation(data_1, data_2)
#' # 3. Run function with temporary file
#' euler_k <- plot_euler_kmers(dataset_sharing_relation)
#' # 4. View result
#' print(euler_k)
plot_euler_kmers <- function(kmer_sharing_relation,
                                path_db = NULL,
                                dataset_label = NULL,
                                color_pallete = "PuBuGn",
                                title_name = "K-mer distribution across datasets"){

  #Checking the number of datasets in the relation list
  #Count dataset quantity
  num_datasets <- length(kmer_sharing_relation)

  #Check if dataset_label (if given) has the same number of datasets as the datasets in dataset
  if(!is.null(dataset_label)){
    if (num_datasets != length(dataset_label)){
      stop("Error: Number or datasets in dataset and given list in dataset label is not the same")
    }
  }

  #Check if dataset quantity is doable
  if(num_datasets > 4){
    stop("Error: Cannot plot more than 7 datasets in a Venn diagram. Pleace reduce the number of datasets")
  }

  #Fit data to the euler diagram
  euler_list <- eulerr::euler(kmer_sharing_relation)

  #If no label is list given, extract the name from the list
  if(is.null(dataset_label)){
    dataset_label <- names(kmer_sharing_relation)
  }

  #Generate color pallete
  message("Generating color palette based on the number of datasets...")

  #Generate num_datasets collors dynamically
  dynamic_colors <- grDevices::hcl.colors(n = num_datasets, palette = color_pallete)

  #Create euler diagram
  message("Generating Euler Diagram...\n")
  euler_cd <- plot(euler_list,
                   quantities = TRUE,
                   labels = dataset_label,
                   fills = list(fill = dynamic_colors, alpha = 0.6),
                   main = title_name)

  return(euler_cd)

}
