select_exclusive_seq <- function(file_1, file_2, path_db = NULL) {
  #If no databse path is given, it creates a temporary file
  if (is.null(path_db)) {
    path_db <- tempfile(fileext = ".duckdb")
  }
  # Creates and connects to database
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path_db)

  #Closes the conection after leaving the function
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

#load tables on duckdb
duckdb::duckdb_register(con, "table_1", file_1)
duckdb::duckdb_register(con, "table_2", file_2)

table_1.data <- dplyr::tbl(con, "table_1")
table_2.data <- dplyr::tbl(con, "table_2")

#Select exclusive sequences
tables_merged <- table_1.data |>
  dplyr::union_all(table_2.data) |>
  dplyr::group_by(seq) |>
  dplyr::filter(dplyr::n() == 1) |>
  dplyr::ungroup() |>
  dplyr::collect()


#Save dataframe in the environment
return(tables_merged)
}
