
rm(list=ls())

source('./code/project_functions.r')
module_init('./simulate_data')
file.copy('./code/simulate_data.r', './simulate_data/code')
file.copy('./code/project_functions.r', './simulate_data/code')
setwd('./simulate_data')
source('./code/simulate_data.r')
setwd('..')

source('./code/project_functions.r')
module_init('./decompose_data')
file.copy('./code/decompose_data.r', './decompose_data/code')
file.copy('./code/project_functions.r', './decompose_data/code')
file.copy('./simulate_data/output/raw_data.csv', './decompose_data/inputs')
setwd('./decompose_data')
source('./code/decompose_data.r')
setwd('..')

source('./code/project_functions.r')
module_init('./analyze_data')
file.copy('./code/analyze_data.r', './analyze_data/code')
file.copy('./code/project_functions.r', './analyze_data/code')
file.copy('./simulate_data/output/raw_data.csv', './analyze_data/inputs')
file.copy('./decompose_data/output/prepped_data.csv', './analyze_data/inputs')
setwd('./analyze_data')
source('./code/analyze_data.r')
setwd('..')

source('./code/project_functions.r')
dir_init('./output')
files <- list.files('./simulate_data/output', full.names=TRUE)
file.copy(files, './output')
files <- list.files('./decompose_data/output', full.names=TRUE)
file.copy(files, './output')
files <- list.files('./analyze_data/output', full.names=TRUE)
file.copy(files, './output')

if(!save_temp){
	unlink( 
		c('./simulate_data',
		'./decompose_data',
		'./analyze_data'),
		recursive=TRUE
	)
}