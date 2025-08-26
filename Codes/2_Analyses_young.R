#### Load Library ####
library(readxl)
library(table1)
library(dplyr)
library(stats)

#### Helper Functions ####

# Factor assignment helper
set_factors <- function(df, factor_levels) {
     df %>%
          mutate(across(names(factor_levels), ~factor(.x, levels = factor_levels[[cur_column()]])))
}

# Reusable table1 helper
make_table1 <- function(vars, group, data, mean_sd = TRUE, overall = TRUE) {
     formula <- as.formula(
          paste("~", paste(vars, collapse = " + "), "|", group)
     )
     
     if(overall) {
          overall <- "Overall"
     }
     
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

# Cutpoint helper: mean + 2*sd
cutpoint <- function(x) {
     mean(x, na.rm = TRUE) + 2 * sd(x, na.rm = TRUE)
}

#### 1. TSPO Young ####

###### Load data ####
TSPO_young <- read_excel("Input_Tables/TSPO_young.xlsx")

# Apply factor levels
factor_levels <- list(
     sex = c("F", "M"),
     apoee4_status = c("Noncarrier", "Carrier"),
     DX2 = c("CU", "CI")
)
TSPO_young <- set_factors(df = TSPO_young, factor_levels = factor_levels)

##### Table ####
demographic_vars_tspo <- c("age_at_mri", "sex", "years_of_education", 
                           "apoee4_status", "plasmaGFAPpgmL_normalized", 
                           "NeoctxAZD_SUVR", "PBR_PCC", "MMSE")

make_table1(vars = demographic_vars_tspo, group = "DX2", data = TSPO_young,
            mean_sd = FALSE, overall = TRUE)

# sTREM2 subset
TSPO_young_sTREM2 <- TSPO_young %>%
     filter(!is.na(sTREM2)) %>%
     mutate(sTREM2 = as.numeric(sTREM2))

mean(TSPO_young_sTREM2$sTREM2, na.rm = TRUE)
sd(TSPO_young_sTREM2$sTREM2, na.rm = TRUE)

# Cutpoints (TSPO PET measures)
list(
     PCC = cutpoint(TSPO_young$PBR_PCC), # PCC
     Difference = cutpoint(TSPO_young$PBR_DKT_difference_composite), # Difference-derived composite brain region
     T_composite = cutpoint(TSPO_young$PBR_DKT_T_composite) # T-derived composite brain region
)

#### 2. sTREM2 Young (TRIAD) ####

###### Load data ####
sTREM2_young_TRIAD <- read_excel("Input_Tables/sTREM2_young_TRIAD.xlsx")
sTREM2_young_TRIAD <- set_factors(df = sTREM2_young_TRIAD, factor_levels = factor_levels)

##### Table ####
demographic_vars_strem2 <- c("age_at_mri", "sex", "years_of_education", 
                             "apoee4_status", "plasmaGFAP_pgmL_normalized", 
                             "NeoctxAZD_SUVR", "sTREM2", "MMSE")
make_table1(vars = demographic_vars_strem2, group = "DX2",
            data = sTREM2_young_TRIAD, mean_sd = FALSE, overall = TRUE)

# PBR PCC distribution for TSPO young sTREM2
mean(TSPO_young_sTREM2$PBR_PCC, na.rm = TRUE)
sd(TSPO_young_sTREM2$PBR_PCC, na.rm = TRUE)

# Cutpoint for sTREM2
cutpoint(sTREM2_young_TRIAD$sTREM2)
