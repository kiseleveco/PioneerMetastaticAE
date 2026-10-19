# Load required libraries
library(purrr)

merge_mcf=function(results_folder)
{
# Step 1: List zip files
zip_files <- list.files(path = results_folder, pattern = "*.zip", full.names = TRUE)

merged_mcf=list()
# Step 2: Process each zip file
for (zip_file in zip_files) {
  # Create a temporary folder to extract files
  temp_extracted_folder <- tempfile()
  dir.create(temp_extracted_folder)
  # Step 2a: Extract zip files
  unzip(zip_file, exdir = temp_extracted_folder)
  # Step 2b: List RDS files within the extracted folder
  rds_files <- list.files(temp_extracted_folder, pattern = "mcf", ignore.case = TRUE, full.names = TRUE)
  # Step 2c: Read RDS files and merge them by filename
  counts=read.csv(paste0(temp_extracted_folder,"/counts.csv"))
  database_name= unique(counts$databaseId)
  merged_mcf[[database_name]]=list()
  for (rds_file in rds_files) {
    data <- arrange_mcf_objects(readRDS(rds_file))
    file_name <- tools::file_path_sans_ext(basename(rds_file))
    merged_mcf[[database_name]][[file_name]]=data
  }
# Clean up extracted folder
  unlink(temp_extracted_folder, recursive = TRUE)
}
return(merged_mcf)
}
