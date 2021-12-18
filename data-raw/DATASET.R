## code to prepare `DATASET` dataset goes here

# batch_deseq ------------------------------------------------------------------

curve_data <- read.csv("data-raw/curve_data.csv",header = FALSE)
time_Elapsed <- read.csv("data-raw/TimeElapsed.csv",header = FALSE)

usethis::use_data(curve_data,time_Elapsed,overwrite = TRUE)

# testdata <- R.matlab::readMat('./data-raw/wormdata.mat')
