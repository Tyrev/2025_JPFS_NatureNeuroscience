#### Load Libraries ####
library(readxl)
library(table1)
library(dplyr)
library(ggplot2)
library(stats)
library(grid)
library(tidyr)

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

# Scatterplot function
scatter_plot <- function(data, group_var, x = "centiloid", y = "plasmaGFAP_z",
                         xlab = "Aβ-PET centiloid", ylab = "Plasma GFAP (z-scored)",
                         colors = c("MA-" = "#B28C36", "MA+" = "#2A3345"),
                         xlim = c(-15, 155), ylim = c(-3, 6)) {
     
     ggplot(data, aes_string(x = x, y = y, color = group_var, fill = group_var)) +
          geom_point(size = 4) +
          geom_smooth(method = "lm", se = TRUE) +
          scale_color_manual(values = colors) +
          scale_fill_manual(values = colors) +
          labs(x = xlab, y = ylab) +
          coord_cartesian(xlim = xlim, ylim = ylim) +
          scale_x_continuous(breaks = seq(0, xlim[2], 50), expand = c(0.035, 0.035)) +
          scale_y_continuous(breaks = seq(ylim[1], ylim[2], 1.5), expand = c(0, 0)) +
          theme_minimal() +
          theme(
               panel.background = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.major.y = element_line(color = "gray", size = 0.2),
               panel.grid.minor = element_blank(),
               axis.line = element_line(color = "black", size = 0.2),
               axis.text = element_text(size = 17, color = "black"),
               axis.title = element_text(size = 17, color = "black"),
               axis.ticks = element_line(color = "black", size = 0.5),
               axis.ticks.length = unit(0.08, "cm"),
               legend.position = "left",
               legend.text = element_text(size = 17, color = "black")
          )
}

# Function for LM + Diagnostic Plots
run_lm <- function(formula, data, show_visualizations = TRUE, visreg_args = NULL) {
     model <- lm(formula, data = data)
     
     # Print summary and CI
     print(summary(model))
     print(confint(model))
     
     if(show_visualizations) {
          # Diagnostic plots
          par(mfrow = c(1,4))
          plot(model, which = 1, add.smooth = FALSE)
          hist(residuals(model), main = "Histogram of Residuals", xlab = "Residuals", freq = FALSE)
          qqnorm(residuals(model), main = "Q-Q Plot of Residuals")
          qqline(residuals(model), lwd = 1)
          plot(model, which = 3, add.smooth = FALSE)
          par(mfrow = c(1,1))
          
          # Optional 2D visreg plot
          if(!is.null(visreg_args)) {
               do.call(visreg2d, c(list(model), visreg_args))
          }
     }
     return(model)
}

#### 1. TRIAD ####

#### 1.1. Load data ####
sTREM2_TRIAD <- read_excel("Input_Tables/sTREM2_SecondarySample_TRIAD.xlsx")

# Apply factor levels
factor_levels_TRIAD <- list(
     sex = c("F", "M"),
     apoee4_status = c("Noncarrier", "Carrier"),
     DX2 = c("CU", "CI"),
     MA_positivity = c("MA-", "MA+")
)
sTREM2_TRIAD <- set_factors(sTREM2_TRIAD, factor_levels_TRIAD)

#### 1.2. Time differences between biomarkers #### 
# List of time-difference columns
time_diff_cols <- c("Correct_Time_difference_CSF_AbPET",
                    "CorrectTime_difference_plasma_AbPET",
                    "Correct_Time_difference_CSF_Plasma")

# Compute median and IQR for all columns
time_diff_summary <- sTREM2_TRIAD %>%
     select(all_of(time_diff_cols)) %>%
     pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
     group_by(variable) %>%
     summarise(
          median = median(value, na.rm = TRUE),
          IQR = IQR(value, na.rm = TRUE),
          .groups = "drop"
     )
time_diff_summary

#### 1.3. Table ####
sTREM2_SecondarySample_TRIAD_filtered <- sTREM2_TRIAD %>%
     filter(DX2 == "CU", age_at_mri >= 65, age_at_mri < 76)

make_table1(
     vars = c("sTREM2", "MA_positivity", "age_at_mri"),
     group = "DX2",
     data = sTREM2_SecondarySample_TRIAD_filtered,
     mean_sd = FALSE,
     overall = TRUE
)

#### 1.4. Stratification ####
sTREM2_TRIAD_by_MA <- sTREM2_TRIAD %>% split(.$MA_positivity)

#### 1.5. Scatterplot ####
scatter_plot(
     data = sTREM2_TRIAD,
     group_var = "MA_positivity",
     x = "Centiloid",
     y = "plasmaGFAP_pgmL_normalized_z"
)

#### 1.6. Regressions ####

##### Model01 ####
Model01 <- run_lm(formula = plasmaGFAP_pgmL_normalized_z ~ Centiloid + age_at_mri + sex + DX2,
                  data = sTREM2_TRIAD_by_MA[["MA-"]],
                  show_visualizations = FALSE)

##### Model02 ####
Model02 <- run_lm(formula = plasmaGFAP_pgmL_normalized_z ~ Centiloid + age_at_mri + sex + DX2,
                  data = sTREM2_TRIAD_by_MA[["MA+"]],
                  show_visualizations = FALSE)

##### Model03 ####
Model03 <- run_lm(formula = plasmaGFAP_pgmL_normalized_z ~ Centiloid*MA_positivity + age_at_mri + sex + DX2,
                  data = sTREM2_TRIAD,
                  show_visualizations = FALSE)

#### 2. WRAP ####

#### 2.1. Load data ####
sTREM2_WRAP <- read_excel("Input_Tables/sTREM2_SecondarySample_WRAP.xlsx")

# Apply factor levels
factor_levels_WRAP <- list(
     Gender = c("F", "M"),
     APOE4_status = c("0", "1"),
     DX2 = c("CU", "CI"),
     AB_positivity = c("0", "1"),
     Batch = c("NTK1", "NTK2"),
     MA_positivity = c("MA-", "MA+")
)
sTREM2_WRAP <- set_factors(sTREM2_WRAP, factor_levels_WRAP)

#### 2.2. Time differences between biomarkers #### 
time_diff_cols <- c("Correct_Difference_age_abPET_sTREM2",
                    "Correct_Difference_age_GFAP_abPET",
                    "Correct_Difference_age_GFAP_sTREM2")

# Compute median and IQR for all columns
time_diff_summary <- sTREM2_WRAP %>%
     select(all_of(time_diff_cols)) %>%
     pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
     group_by(variable) %>%
     summarise(
          median = median(value, na.rm = TRUE),
          IQR = IQR(value, na.rm = TRUE),
          .groups = "drop"
     )
time_diff_summary

#### 2.3. Batch adjustment ####
Model_adj <- lm(sTREM2 ~ Batch, data = sTREM2_WRAP)
sTREM2_WRAP$sTREM2_adj <- resid(Model_adj) + mean(sTREM2_WRAP$sTREM2, na.rm = TRUE)

#### 2.4. Percentile Computation ####
sTREM2_WRAP_filtered <- sTREM2_WRAP %>%
     filter(DX2 == "CU", AgeAtVisit >= 65, AgeAtVisit < 76)

triad_ecdf <- ecdf(sTREM2_SecondarySample_TRIAD_filtered$sTREM2)
triad_ecdf(2524.014) ## 0.423

quantile(sTREM2_WRAP_filtered$sTREM2_adj, probs = 0.423, type = 7)

#### 2.5. Stratification ####
sTREM2_WRAP_by_MA <- sTREM2_WRAP %>% split(.$MA_positivity)

#### 2.6. Scatterplot ####
scatter_plot(sTREM2_WRAP, group_var = "MA_positivity", x = "centiloid", y = "PlasmaGFAP_pgmL_z", xlim = c(-15, 110))

#### 2.7. Regressions ####

##### Model04 ####
Model04 <- run_lm(formula = PlasmaGFAP_pgmL_z ~ centiloid + AgeAtVisit + Gender + DX2,
                  data = sTREM2_WRAP_by_MA[["MA-"]],
                  show_visualizations = FALSE)

##### Model05 ####
Model05 <- run_lm(formula = PlasmaGFAP_pgmL_z ~ centiloid + AgeAtVisit + Gender + DX2,
                  data = sTREM2_WRAP_by_MA[["MA+"]],
                  show_visualizations = FALSE)

##### Model06 ####
Model06 <- run_lm(formula = PlasmaGFAP_pgmL_z ~ centiloid*MA_positivity + AgeAtVisit + Gender + DX2,
                  data = sTREM2_WRAP,
                  show_visualizations = FALSE)
