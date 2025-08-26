#### Load Libraries ####
library(readxl)
library(table1)
library(dplyr)
library(ggplot2)
library(visreg)
library(broom)
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

# Function for LM + Residual Plots
run_lm <- function(formula, data, show_visualizations = TRUE, visreg_args = NULL) {
     model <- lm(formula, data = data)
     
     # Print summary and CI
     print(summary(model))
     print(confint(model))
     
     if(show_visualizations) {
          # Residual plots
          par(mfrow = c(1,3))
          hist(residuals(model), main = "Histogram of Residuals", xlab = "Residuals", freq = FALSE)
          qqnorm(residuals(model), main = "Q-Q Plot of Residuals")
          qqline(residuals(model), lwd = 1)
          plot(model, which = 3)
          par(mfrow = c(1,1))
          
          # Optional 2D visreg plot
          if(!is.null(visreg_args)) {
               do.call(visreg2d, c(list(model), visreg_args))
          }
     }
     return(model)
}

# Scatterplot function
scatter_plot <- function(data, group_var, 
                         x = "plasmaGFAPpgmL_normalized", 
                         y = "Ptau217_Conc_pgmL",
                         xlab = "Plasma GFAP (pg/mL)", 
                         ylab = "Plasma p-tau217 (pg/mL)",
                         colors = c("MA-" = "#B28C36", "MA+" = "#2A3345"),
                         xlim = c(0, 650), 
                         ylim = c(-0.05, 0.5)) {
     
     ggplot(data, aes_string(x = x, y = y, color = group_var, fill = group_var)) +
          geom_point(size = 4) +
          geom_smooth(method = "lm", se = TRUE) +
          scale_color_manual(values = colors) +
          scale_fill_manual(values = colors) +
          labs(x = xlab, y = ylab) +
          coord_cartesian(xlim = xlim, ylim = ylim) +
          scale_x_continuous(breaks = seq(0, xlim[2], 150), expand = c(0.035, 0.035)) +
          scale_y_continuous(breaks = seq(0, ylim[2], 0.1), expand = c(0, 0)) +
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

make_sem_model <- function(outcome, predictor) {
     paste0(
          outcome, " ~ ", predictor, " + Ptau217_Conc_pgmL + plasmaGFAPpgmL_normalized + ",
          "NeoctxAZD_SUVR + age_at_mri + sex + years_of_education\n",
          
          predictor, " ~ Ptau217_Conc_pgmL + plasmaGFAPpgmL_normalized + ",
          "NeoctxAZD_SUVR + age_at_mri + sex\n",
          
          "Ptau217_Conc_pgmL ~ plasmaGFAPpgmL_normalized + NeoctxAZD_SUVR + age_at_mri + sex\n",
          
          "plasmaGFAPpgmL_normalized ~ NeoctxAZD_SUVR + age_at_mri + sex"
     )
}

run_sem <- function(model_string, data, group_label = NULL) {
     fit <- sem(model_string, data = data)
     
     cat("\n==============================================\n")
     if (!is.null(group_label)) {
          cat("SEM Results for:", group_label, "\n")
     }
     cat("==============================================\n")
     
     # Full summary
     print(summary(fit, standardized = TRUE, fit.measures = TRUE, rsq = TRUE, ci = TRUE))
     
     # Key fit indices
     print(fitMeasures(fit, c("cfi", "rmsea", "srmr", "tli")))
     
     # Standardized paths
     print(standardizedSolution(fit, type = "std.all"))
     
     return(fit)
}

#### 1. TSPO_PrimarySample_ptau217 ####

###### Load data ####
TSPO_PrimarySample_ptau217 <- read_excel("Input_Tables/TSPO_PrimarySample_ptau217.xlsx")

# Apply factor levels
factor_levels <- list(
     sex = c("F", "M"),
     DX2 = c("CU", "CI"),
     apoee4_status = c("Noncarrier", "Carrier"),
     MA_positivity = c("MA-", "MA+")
)
TSPO_PrimarySample_ptau217 <- set_factors(TSPO_PrimarySample_ptau217, factor_levels)

##### Table ####
make_table1(vars = c("Ptau217_Conc_pgmL", "plasmaGFAPpgmL_normalized", "PBR_PCC", "age_at_mri", "sex"),
            group = "DX2", data = TSPO_PrimarySample_ptau217,
            mean_sd = FALSE, overall = TRUE)

#### 2. Time differences between biomarkers #### 
# List of time-difference columns
time_diff_cols <- c("Correct_Time_diff_abPET_plasmaPtau217",
                    "Correct_Time_diff_tspoPET_plasmaPtau217",
                    "Correct_Time_diff_plasmaGFAP_plasmaPtau217",
                    "Correct_Time_diff_plasmaPtau217_tauPET",
                    "Correct_Time_diff_plasmaPtau217_MMSE")

# Compute median and IQR for all columns
time_diff_summary <- TSPO_PrimarySample_ptau217 %>%
     select(all_of(time_diff_cols)) %>%
     pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
     group_by(variable) %>%
     summarise(
          median = median(value, na.rm = TRUE),
          IQR = IQR(value, na.rm = TRUE),
          .groups = "drop"
     )
time_diff_summary

#### 3. Stratification #### 
TSPO_PrimarySample_ptau217_by_MA <- TSPO_PrimarySample_ptau217 %>%
     split(.$MA_positivity)

#### 4. Scatterplot ####
scatter_plot(TSPO_PrimarySample_ptau217, group_var = "MA_positivity")

#### 5. Regressions ####

##### Model01 ####
Model01 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ plasmaGFAPpgmL_normalized_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_ptau217_by_MA[["MA-"]])

##### Model02 ####
Model02 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ plasmaGFAPpgmL_normalized_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_ptau217_by_MA[["MA+"]])

##### Model03 ####
Model03 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ plasmaGFAPpgmL_normalized_z*MA_positivity + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_ptau217)

##### Model04 ####
Model04 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ plasmaGFAPpgmL_normalized_z*PBR_PCC_z + age_at_mri + sex + DX2, 
                  data = TSPO_PrimarySample_ptau217,
                  visreg_args = list(xvar = "plasmaGFAPpgmL_normalized_z",
                                     yvar = "PBR_PCC_z", nn = 350,
                                     plot.type = "image",
                                     col = colorRampPalette(c("#B28C36", "#DFCBA3", "#F3F5F8", 
                                                              "#C3CFE1","#94A4C0", "#5F7498", 
                                                              "#3B485F", "#2A3345", "#070A0D"))(200)))

##### 5.1 FDR adjustment ####
# Define models and coefficient names
models <- list(Model01, Model02, Model03, Model04)
coef_names <- c(
     "plasmaGFAPpgmL_normalized_z",
     "plasmaGFAPpgmL_normalized_z",
     "plasmaGFAPpgmL_normalized_z:MA_positivityMA+",
     "plasmaGFAPpgmL_normalized_z:PBR_PCC_z"
)

# Extract p-values
p_values_02 <- sapply(seq_along(models), function(i) {
     coef(summary(models[[i]]))[coef_names[i], "Pr(>|t|)"]
})

# FDR adjustment
adjusted_p_values_02 <- p.adjust(p_values_02, method = "fdr")
adjusted_p_values_02

##### 5.2 ANOVA ####

##### Model05 ####
Model05 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ plasmaGFAPpgmL_normalized_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_ptau217)

##### Model06 ####
Model06 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ PBR_PCC_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_ptau217)

##### Model07 ####
Model07 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ plasmaGFAPpgmL_normalized_z + PBR_PCC_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_ptau217)

##### Model08 ####
Model08 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ plasmaGFAPpgmL_normalized_z * PBR_PCC_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_ptau217)

anova(Model08, Model05)
anova(Model08, Model06)
anova(Model08, Model07)

#### 6. Regressions ####
library(lavaan)
library(semPlot)

##### Run SEMs ####
SEM_model <- make_sem_model("MMSE", "temporal_meta_ROI_infCG")

##### Model09 ####
Model09 <- run_sem(model_string = SEM_model,
                 data = TSPO_PrimarySample_ptau217_by_MA[["MA-"]],
                 group_label = "MA-")

##### Model10 ####
Model10 <- run_sem(model_string = SEM_model,
                 data = TSPO_PrimarySample_ptau217_by_MA[["MA+"]],
                 group_label = "MA+")

#### 7. Sensitivity/Exploratory Analyses ####
# List of (outcome, predictor) combos
sem_specs <- list(
     list(outcome = "MMSE", predictor = "Braak_12_SUVR_infCG"),
     list(outcome = "MMSE", predictor = "Braak_34_SUVR_infCG"),
     list(outcome = "MMSE", predictor = "Braak_56_SUVR_infCG"),
     list(outcome = "Memory_composite_score", predictor = "temporal_meta_ROI_infCG")
)

# Drop missing values for Memory_composite_score beforehand
TSPO_PrimarySample_ptau217_by_MA[["MA-"]] <- TSPO_PrimarySample_ptau217_by_MA[["MA-"]] %>%
     filter(!is.na(Memory_composite_score))
TSPO_PrimarySample_ptau217_by_MA[["MA+"]] <- TSPO_PrimarySample_ptau217_by_MA[["MA+"]] %>%
     filter(!is.na(Memory_composite_score))

# Run all SEMs
fits <- list()
for (spec in sem_specs) {
     outcome <- spec$outcome
     predictor <- spec$predictor
     model_string <- make_sem_model(outcome, predictor)
     
     # MA-
     fit_neg <- run_sem(model_string, 
                        data = TSPO_PrimarySample_ptau217_by_MA[["MA-"]], 
                        group_label = paste0("MA- | ", outcome, " ~ ", predictor))
     
     # MA+
     fit_pos <- run_sem(model_string, 
                        data = TSPO_PrimarySample_ptau217_by_MA[["MA+"]], 
                        group_label = paste0("MA+ | ", outcome, " ~ ", predictor))
     
     # Store
     fits[[paste0("MA-", outcome, predictor)]] <- fit_neg
     fits[[paste0("MA+", outcome, predictor)]] <- fit_pos
}

##### Model11 ####
Model11 <- fits[[1]]

##### Model12 ####
Model12 <- fits[[2]]

##### Model13 ####
Model13 <- fits[[3]]

##### Model14 ####
Model14 <- fits[[4]]

##### Model15 ####
Model15 <- fits[[5]]

##### Model16 ####
Model16 <- fits[[6]]

##### Model17 ####
Model17 <- fits[[7]]

##### Model18 ####
Model18 <- fits[[8]]

##### 8.3 Including Outliers ####

###### Stratification #### 
TSPO_PrimarySample_ptau217_Outliers <- read_excel(
     "Input_Tables/TSPO_PrimarySample_ptau217_Outliers.xlsx"
) %>%
     mutate(
          sex = factor(sex, levels = c("F", "M")),
          DX2 = factor(DX2, levels = c("CU", "CI")),
          MA_positivity = factor(MA_positivity, levels = c("MA-", "MA+"))
     )
TSPO_PrimarySample_ptau217_Outliers_MA <- TSPO_PrimarySample_ptau217_Outliers %>%
     split(.$MA_positivity)

###### Model19 ####
Model19 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ plasmaGFAPpgmL_normalized_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_ptau217_Outliers_MA[["MA-"]],
                  show_visualizations = FALSE)

###### Model20 ####
Model20 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ plasmaGFAPpgmL_normalized_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_ptau217_Outliers_MA[["MA+"]],
                  show_visualizations = FALSE)

###### Model21 ####
Model21 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ plasmaGFAPpgmL_normalized_z * MA_positivity + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_ptau217_Outliers,
                  show_visualizations = FALSE)

###### Model26 ####
Model26 <- run_lm(formula = Ptau217_Conc_pgmL_z ~ plasmaGFAPpgmL_normalized_z * PBR_PCC_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_ptau217_Outliers,
                  show_visualizations = FALSE)
