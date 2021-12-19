usethis::use_package("Rcpp")
usethis::use_package("microbenchmark")
usethis::use_package("boot")
usethis::use_package("bootstrap")
usethis::use_package("stats")
usethis::use_package("base")
usethis::use_package("utils")
usethis::use_package("FNN")
usethis::use_package("energy")
usethis::use_package("Ball")




usethis::use_rcpp()
usethis::use_mit_license("Rongkang Xiong")
usethis::use_readme_md()
usethis::use_code_of_conduct("earth@mail.ustc.edu.cn")

# Vignette workflow
usethis::use_vignette("homework")
usethis::use_vignette("intro")



# 导出vignettes
devtools::build_vignettes()


# CRAN
# 拼写检查
devtools::spell_check()
# 常规本地测试
devtools::check()
# Windows平台测试
devtools::check_win_devel()

# 检查
devtools::release()



# 导入包运行
devtools::load_all()
