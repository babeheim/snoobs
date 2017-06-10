
rm(list=ls())

source('./code/project_functions.r')

dir_init('./temp')

file.copy('./code/rtdice_cliodynamics.tex', './temp')
files <- list.files('./inputs', full.names=TRUE)
file.copy(files, './temp')

setwd('./temp')
system("pdflatex rtdice_cliodynamics.tex")
system("pdflatex rtdice_cliodynamics.tex")
setwd('..')

file.copy('./code/rtdice_cliodynamics_esm.tex', './temp')
files <- list.files('./inputs', full.names=TRUE)
file.copy(files, './temp')

setwd('./temp')
system("pdflatex rtdice_cliodynamics_esm.tex")
system("pdflatex rtdice_cliodynamics_esm.tex")
setwd('..')

dir_init('./output')
file.copy('./temp/rtdice_cliodynamics.pdf', './output')
file.copy('./temp/rtdice_cliodynamics_esm.pdf', './output')

if(!save_temp) unlink('./temp', recursive=TRUE)