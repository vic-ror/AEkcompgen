#' Creates an Upset Plot based on the kmer sequence distribuition among datasets.
#' The suggested maximum number of datasets is 4 but can hold it up to 7 datasets, if there are more markes it c
#' @param kmer_sharing_relation A li from the kmer_sharing_relation function
#' @param dataset_label A list with the names of the datasets wished to be in the diagram, KEEP THE datasetS' ORDER IN MIND
#' @param color_pallete A color pallete from the scale_fill_gradientn function, default is "PuBuGn"
#' @param empty_intersection If the empty intersectios are to be displayed put "on" or not, default is NULL
#'
#' @return An Upset plot based on the kmer distribuition among datasets
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
#' upset_k <- plot_upset_kmers(dataset_sharing_relation)
#' # 4. View result
#' print(upset_k)

plot_upset_kmers <- function(kmer_sharing_relation,
                             dataset_label = NULL,
                             color_pallete = "PuBuGn",
                             empty_intersection = NULL){

  #Checking the number of datasets in the relation list
  #Count dataset quantity
  num_datasets <- length(kmer_sharing_relation)


  #Check if dataset_label (if given) has the same number of datasets as the datasets in dataset
  if(!is.null(dataset_label)){
    if (num_datasets != length(dataset_label)){
      stop("Error: Number or datasets in dataset and given list in dataset label is not the same")
    }
    #Rename dataset_labes
    names(upset_list) <- dataset_label
  }

  #Turn list into the upset plot format
  upset_list <- UpSetR::fromList(kmer_sharing_relation)

  #Extract collors from color pallete
  palette_colors <- RColorBrewer::brewer.pal(9, color_pallete)
  dark_color <- palette_colors[8]
  mid_color  <- palette_colors[4]


  #Plotting
  message("Generating UpsetPlot...")
  upset_plot <- UpSetR::upset(data = upset_list,
                              order.by = "freq",
                              nsets = num_datasets,
                              point.size = 3.5,
                              mainbar.y.label = "K-mer quantity\n",
                              sets.x.label = "K-mers per dataset",
                              main.bar.color = dark_color,
                              sets.bar.color = mid_color,
                              matrix.color = dark_color,
                              empty.intersections = empty_intersection,
                              text.scale=c(1.3, 1.3, 1, 1, 1.3, 1.5)
  )


    return(upset_plot)

}
