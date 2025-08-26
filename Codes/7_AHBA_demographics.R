#### Load Libraries ####
library(readxl)
library(table1)
library(dplyr)
library(ggplot2)
library(stats)
library(grid)

#### Helper Functions ####

# Factor assignment helper
set_factors <- function(df, factor_levels) {
     df %>%
          mutate(across(names(factor_levels), ~factor(.x, levels = factor_levels[[cur_column()]])))
}

# Reusable table1 helper
make_table1 <- function(vars, group, data, mean_sd = TRUE, overall = TRUE) {
     if (!is.null(group)) {
          formula <- as.formula(paste("~", paste(vars, collapse = " + "), "|", group))
     } else {
          formula <- as.formula(paste("~", paste(vars, collapse = " + ")))
     }
     
     if(overall) overall <- "Overall"
     
     if (mean_sd) {
          table1(formula,
                 data = data,
                 overall = overall,
                 render.continuous = function(x) sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)))
     } else {
          table1(formula,
                 data = data,
                 overall = overall)
     }
}

#### 1. AHBA_donors ####

###### Load data ####
AHBA_donors <- read_excel("Input_Tables/AHBA_donors.xlsx")

# Apply factor levels
factor_levels <- list(
     Sex = c("Female", "Male")
)
AHBA_donors <- set_factors(AHBA_donors, factor_levels)

##### Table ####
make_table1(vars = c("Age_y", "Sex", "PMI"), 
            group = NULL, 
            data = AHBA_donors, 
            mean_sd = FALSE, 
            overall = TRUE)
