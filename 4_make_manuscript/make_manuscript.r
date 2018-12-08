
if (scaffold) {
  rm(list = ls())
  source("../project_support.r")
}

dir_init("./temp")

file.copy('./tex/rtdice_cliodynamics.tex', './temp')

files <- list.files('./assets', full.names=TRUE)
file.copy(files, './temp')

setwd('./temp')
system("pdflatex rtdice_cliodynamics.tex")
system("pdflatex rtdice_cliodynamics.tex")
setwd('..')

file.copy('./tex/rtdice_cliodynamics_esm.tex', './temp')

files <- list.files('./assets', full.names=TRUE)
file.copy(files, './temp')

setwd('./temp')
system("pdflatex rtdice_cliodynamics_esm.tex")
system("pdflatex rtdice_cliodynamics_esm.tex")
setwd('..')

dir_init('./output')
file.copy('./temp/rtdice_cliodynamics.pdf', './output')
file.copy('./temp/rtdice_cliodynamics_esm.pdf', './output')

if(!save_temp) unlink('./temp', recursive=TRUE)
