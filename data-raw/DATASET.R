## code to prepare `DATASET` dataset goes here

# batch_deseq ------------------------------------------------------------------

curve_data <- read.csv("data-raw/curve_data.csv", header = FALSE)
curve_data <- data.matrix(curve_data) # 将dataframe转换为矩阵

time_Elapsed <- read.csv("data-raw/TimeElapsed.csv", header = FALSE)
time_Elapsed <- time_Elapsed$V1


w <- read.csv("data-raw/w.csv", header = FALSE)
w <- data.matrix(w)
w <- t(w)


usethis::use_data(curve_data, compress = "xz", overwrite = TRUE)
usethis::use_data(time_Elapsed, compress = "xz", overwrite = TRUE)
usethis::use_data(w, compress = "xz", overwrite = TRUE)
