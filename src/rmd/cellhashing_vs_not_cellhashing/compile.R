#!/usr/bin/env Rscript

file.remove("_main.Rmd")
# source("styler.R")
bookdown::render_book("index.Rmd")
