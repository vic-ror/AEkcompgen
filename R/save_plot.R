#' Saves plots.
#' In pdf or png format.
#' @param plot A plot.
#' @param file_name Name to save the plot.
#' @param format Format to save the file, between "pdf" and "png".
#' @param width Width of the pdf or of the png file
#' @param height Height of the pdf or of the png file
#' @return A Euler diagram based on the cluster distribuition among markers
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # 1. Create temporary plot
#' plot_test <- plot(1, type = "n",
#'     xlim = c(0, 10), ylim = c(0, 10),
#'     xlab = "X Axis Test", ylab = "Y Axis Test",
#'     main = "Empty Test Plot")
#' # 2. Run function with temporary file
#' save_plot(plot_test, file_name = "plot_test", format = "pdf")
#' # 4. Delete temporary file
#' unlink("plot_test.pdf")
save_plot <- function(plot, file_name, format = "pdf", width = NULL, height = NULL){
  if(format == "pdf"){
    if(is.null(width)){
      width = 8
    }
    if(is.null(height)){
      height = 5.5
    }
    width_value = width
    height_value = height
    grDevices::pdf(file = glue::glue("{file_name}.pdf"), width = width_value, height = height_value)
    print({plot})
    grDevices::dev.off()
  }

  else if (format == "png"){
    if(is.null(width)){
      width = 800
    }
    if(is.null(height)){
      height = 600
    }
    width_value = width
    height_value = height
    grDevices::png(filename = glue::glue("{file_name}.png"), width = width_value, height = height_value, res = 150)
    print({plot})
    grDevices::dev.off()
  }
  else{
  message("Please give a format to save the plot (pdf/png).")
  }
}
