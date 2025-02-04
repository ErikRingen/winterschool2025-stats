To run this code, you will need to [install R](https://www.r-project.org/) and the following R packages:
  
```
install.packages(c(
  "tidyverse",
  "patchwork",
  "ggsci",
  "RColorBrewer",
  "cubelyr",
  "proxy",
  "foreach",
  "doParallel")
```

You will also require the `rethinking` package, which you can install with the following code (see [here](https://github.com/rmcelreath/rethinking) for more details):
  
  ```
# from https://github.com/rmcelreath/rethinking
install.packages(c("coda", "mvtnorm", "devtools", "loo", "dagitty"))
devtools::install_github("rmcelreath/rethinking")
```

Analyses from this manuscript were run in R version 4.2.1., with models fit using RStan version 2.26.13. Please be aware that this pipeline takes a very long time to run, somewhere between 7 and 10 days depending on your system.

### Executing code

1. Set your working directory in R to this code repository
2. Load the `targets` package with `library(targets)`
3. To run all analyses, run `tar_make()`
4. To load individual targets into your environment, run `tar_load(targetName)`

To visualize how the various code and functions work together in the pipeline, run `targets::tar_visnetwork()`.
