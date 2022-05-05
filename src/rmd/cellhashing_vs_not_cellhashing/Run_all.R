rmd_list <- list.files()[grep(".Rmd", list.files())]

r_list <- gsub(pattern = ".Rmd", replacement = ".R", x = rmd_list)
dir.create("R_files")

for (file in 1:length(r_list)) {
  print(rmd_list[file])
  print(r_list[file])
  knitr::purl(rmd_list[file], output = paste0("R_files/", r_list[file]))
}

source(paste0("R_files/", r_list[length(r_list)]))
sapply(paste0("R_files/", r_list[1:length(r_list) - 1]), source)
