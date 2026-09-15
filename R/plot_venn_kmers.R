#' Creates a Venn Diagram based on the kmer sequence distribuition among datasets.
#' The suggested maximum number of datasets is 4 but can hold it up to 7 datasets, if there are more datasets the diagram may not make sense
#' @param kmer_sharing_relation A list from the kmer_sharing_relation function
#' @param dataset_label A list with the names of the datasets wished to be in the diagram, KEEP THE datasets' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#' @param plot_title The title for the plot.
#'
#' @return A Venn Diagram based on the kmer distribuition among datasets
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary dataframes
#' data_1 <- data.frame(id = c("tag1_1_1", "tag1_2_2", "tag1_3_3", "tag1_4_4", "tag1_5_5"),
#'  seq = c("AAAG", "AGAG", "CTGC", "AAAA", "GAGA"))
#' data_2 <- data.frame(id = c("tag2_1_1", "tag2_2_2", "tag2_3_3", "tag2_4_4", "tag2_5_5"),
#'  seq = c("AATG", "AGAG", "CAGC", "AATA", "GAGA"))
#' # 2. Run kmer_sharing_relation function to get list
#' dataset_sharing_relation <- kmer_sharing_relation(data_1, data_2)
#' # 3. Run function with temporary file
#' venn_k <- plot_venn_kmers(dataset_sharing_relation)
#' # 4. View result
#' print(venn_k)
plot_venn_kmers <- function(kmer_sharing_relation,
                               dataset_label = NULL,
                               color_pallete = "PuBuGn",
                               plot_title = "K-mer distribution across datasets"){
  #Checking the number of datasets in the relation list
  #Count dataset quantity
  num_datasets <- length(kmer_sharing_relation)

  #Check if dataset quantity is the same as in the list
  if(!is.null(dataset_label)){
    if (num_datasets != length(dataset_label)){
      stop("Error: Number or datasets in dataset and given list in dataset label is not the same")
    }
  }

  #Check if dataset quantity is doable
  if(num_datasets > 4 & num_datasets <= 7){
    message("The suggested maximum number of datasets is 4 but can hold it up to 7 datasets.")
  }

  else if(num_datasets > 7){
    message("If number of datasets is larger than 7, it will give an upset plot")
  }

  #If no label is list given, extract the name from the list
  if(is.null(dataset_label)){
    dataset_label <- names(kmer_sharing_relation)
  }

  #Create gradient of collors based on given the pallete
  color_gradient <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, color_pallete))(100)

  #Get the upper limit for gradient
  max_val <- max(sapply(kmer_sharing_relation, length))

  #Create venn diagram
  message("Generating Venn Diagram...\n")
  venn <- ggVennDiagram::ggVennDiagram(kmer_sharing_relation,
                                       edge_size = 0.6,
                                       label_alpha = 0,
                                       label_percent_digit = 0,
                                       category.names = dataset_label) +
    ggplot2::scale_fill_gradientn(colors = color_gradient, limits = c(0, max_val)) +
    ggplot2::labs(title = plot_title)



  return(venn)
}
