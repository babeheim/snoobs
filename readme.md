snoobs
============

Analysis script for "Evolutionary Decomposition and the Mechanisms of Cultural Change", by Bret Beheim and Ryan Baldini, in Cliodynamics, Volume 3: 217–233, available at https://escholarship.org/uc/item/5w49c6wt

Requirements:
- R (3.3.1 or greater) https://cran.r-project.org/
- rethinking package (v1.59 or greater), http://xcelab.net/rm/software/
- bbmle package, available on CRAN
- LaTeX, https://www.latex-project.org/

Instructions:

In R, set the working directory to that containing this readme file. For example, on a Mac or Linux machine, you might say

```
    setwd('~/Desktop/snoobs')
```

if the folder containing the project is named 'snoobs' and on your Desktop. You can tell if you are in the right place by typing in `dir()` and seeing the folders 'code' and this readme.txt file.

The analysis itself is broken up into independent modules that pass outputs to each other. The whole process runs by typing one command into R,

```
    source('./code/run_project.r')
```

with the project folder as the working directory. If all goes well, each step of the analysis will execute in sequence, and write the paper's tables and figures into an 'output' folder.

By default the analysis will delete all temporary files and folders, but if you want to see all intermediate steps you can disable this by flipping the `save_temp` variable in 'project_variables.r' from `FALSE` to `TRUE`.

The total time until completion will vary by machine. It takes about 30 minutes for me.

The project is maintained by Bret Beheim (beheim@gmail.com) and is hosted at https://github.com/babeheim.