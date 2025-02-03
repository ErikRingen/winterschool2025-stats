# _targets.R file
library(targets)

tar_source() # load the R scripts in my dir
tar_option_set(packages = c("tidyverse", "marginaleffects")) # load any packages used
# define pipeline
list(
    tar_target(data_file, "data/cognitive_flexibility.csv", format = "file"),
    tar_target(data, read_csv(data_file)),
    tar_target(model, fit_model(data)),
    tar_target(treatment_effect, marginaleffects(model, newdata = data, variables = "condition")),
    tar_target(slopes_df, as.data.frame(slopes))
)
