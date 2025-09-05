#### Load Library ####
library(ggplot2)
library(grid)
library(readxl)
library(dplyr)
library(table1)
library(tidyr)
library(visreg)
library(broom)

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
     }
     
     # Optional 2D visreg plot
     if(!is.null(visreg_args)) {
          do.call(visreg2d, c(list(model), visreg_args))
     }
     
     return(model)
}

# Function to run regression per region
run_region_regressions <- function(data, regions, outcome_y, interaction_var = NULL) {
     results <- lapply(regions, function(region) {
          # Build formula
          if (is.null(interaction_var)) {
               formula <- as.formula(paste(outcome_y, "~", region, "+ age_at_mri + sex + DX2"))
          } else {
               formula <- as.formula(paste(outcome_y, "~", paste0(region, "*", interaction_var), "+ age_at_mri + sex + DX2"))
          }
          
          model <- lm(formula, data = data)
          model_tidy <- tidy(model)
          
          # Extract rows of interest
          if (is.null(interaction_var)) {
               model_tidy <- model_tidy[model_tidy$term == region, c("term", "estimate", "p.value","statistic")]
          } else {
               # Only keep interaction term(s)
               interaction_pattern <- paste0(region, ":")
               model_tidy <- model_tidy[grep(interaction_pattern, model_tidy$term), c("term", "estimate", "p.value","statistic")]
          }
          
          cbind(Region = region, model_tidy)
     })
     
     results_df <- do.call(rbind, results)
     
     # FDR adjustment
     results_df$FDR_P <- p.adjust(results_df$p.value, method = "fdr")
     
     # Rename
     colnames(results_df)[colnames(results_df) == "estimate"] <- "Coefficient"
     colnames(results_df)[colnames(results_df) == "p.value"] <- "P_Value"
     colnames(results_df)[colnames(results_df) == "statistic"] <- "tvalue"
     
     return(results_df)
}

# Scatterplot function
scatter_plot <- function(data, group_var, 
                         x = "NeoctxAZD_SUVR", 
                         y = "plasmaGFAPpgmL_normalized",
                         xlab = "Aβ-PET SUVR", 
                         ylab = "Plasma GFAP (pg/mL)",
                         colors = c("MA-" = "#B28C36", "MA+" = "#2A3345"),
                         xlim = c(1, 3.5), 
                         ylim = c(0, 750)) {
     
     ggplot(data, aes_string(x = x, y = y, color = group_var, fill = group_var)) +
          geom_point(size = 5) +
          geom_smooth(method = "lm", se = TRUE) +
          scale_color_manual(values = colors) +
          scale_fill_manual(values = colors) +
          labs(x = xlab, y = ylab) +
          coord_cartesian(xlim = xlim, ylim = ylim) +
          scale_x_continuous(breaks = seq(0, xlim[2], 0.5), expand = c(0.035, 0.035)) +
          scale_y_continuous(breaks = seq(0, ylim[2], 150), expand = c(0, 0)) +
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

#### 1. TSPO_PrimarySample #### 

###### Load data ####
TSPO_PrimarySample <- read_excel("Input_Tables/TSPO_PrimarySample.xlsx")

# Apply factor levels
factor_levels <- list(
     sex = c("F", "M"),
     DX2 = c("CU", "CI"),
     apoee4_status = c("Noncarrier", "Carrier"),
     AB_threshold = c("0", "1"),
     MA_positivity = c("MA-", "MA+"),
     MA_positivity_DKT_difference_composite = c("MA-", "MA+"),
     MA_positivity_DKT_T_composite = c("MA-", "MA+")
)
TSPO_PrimarySample <- set_factors(TSPO_PrimarySample, factor_levels)

# Define variables for table1
demographic_vars <- c(
     "age_at_mri", "sex", "years_of_education", "apoee4_status", 
     "plasmaGFAPpgmL_normalized", "NeoctxAZD_SUVR", 
     "temporal_meta_ROI_infCG", "PBR_PCC", "MMSE"
)

###### Table by DX #####
make_table1(vars = c(demographic_vars, "MA_positivity"), group = "DX2", 
            data = TSPO_PrimarySample,
            mean_sd = FALSE, overall = TRUE)

###### Table by MA_positivity ##### 
make_table1(vars = c(demographic_vars, "DX2"), group = "MA_positivity", 
            data = TSPO_PrimarySample,
            mean_sd = FALSE, overall = TRUE)

#### 2. Time differences between biomarkers #### 
# List of time-difference columns
time_diff_cols <- c(
     "Correct_Time_diff_abPET_tspoPET",
     "Correct_Time_diff_abPET_plasmaGFAP",
     "Correct_Time_diff_abPET_tauPET",
     "Correct_Time_diff_abPET_MMSE",
     "Correct_Time_diff_tspoPET_plasmaGFAP",
     "Correct_Time_diff_tspoPET_tauPET",
     "Correct_Time_diff_tspoPET_MMSE",
     "Correct_Time_diff_plasmaGFAP_tauPET",
     "Correct_Time_diff_plasmaGFAP_MMSE",
     "Correct_Time_diff_tauPET_MMSE"
)

# Compute median and IQR for all columns
time_diff_summary <- TSPO_PrimarySample %>%
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
TSPO_by_MA <- TSPO_PrimarySample %>%
     split(.$MA_positivity)

#### 4. Scatter Plot #### 
scatter_plot(
     data = TSPO_PrimarySample,
     group_var = "MA_positivity"
)

#### 5. "Overall" Regression Analyses ####

##### Model01 ####
Model01 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + sex + DX2,
                  data = TSPO_by_MA[["MA-"]],
                  show_visualizations = TRUE)

##### Model02 ####
Model02 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + sex + DX2,
                  data = TSPO_by_MA[["MA+"]],
                  show_visualizations = TRUE)

##### Model03 ####
Model03 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z*MA_positivity + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample,
                  show_visualizations = TRUE)

##### Model04 ####
Model04 <- run_lm(
     plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z*PBR_PCC_z + age_at_mri + sex + DX2,
     data = TSPO_PrimarySample, 
     show_visualizations = TRUE,
     visreg_args = list(
          x = "NeoctxAZD_SUVR_z",
          y = "PBR_PCC_z",
          nn = 350,
          plot.type = "image",
          col = colorRampPalette(c("#B28C36", "#DFCBA3", "#F3F5F8",
                                   "#C3CFE1","#94A4C0", "#5F7498",
                                   "#3B485F", "#2A3345", "#070A0D"))(200)
     )
)

##### FDR-adjusted P-values####

# Define models and coefficient names
models <- list(Model01, Model02, Model03, Model04)
coef_names <- c(
     "NeoctxAZD_SUVR_z",
     "NeoctxAZD_SUVR_z",
     "NeoctxAZD_SUVR_z:MA_positivityMA+",
     "NeoctxAZD_SUVR_z:PBR_PCC_z"
)

# Extract p-values
p_values_01 <- sapply(seq_along(models), function(i) {
     coef(summary(models[[i]]))[coef_names[i], "Pr(>|t|)"]
})

# FDR adjustment
adjusted_p_values_01 <- p.adjust(p_values_01, method = "fdr")
adjusted_p_values_01

#### 6. ANOVA comparisons ####

##### Model05 ####
Model05 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample,
                  show_visualizations = FALSE)

##### Model06 ####
Model06 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ PBR_PCC_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample,
                  show_visualizations = FALSE)

##### Model07 ####
Model07 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + PBR_PCC_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample,
                  show_visualizations = FALSE)

##### Model08 ####
Model08 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z*PBR_PCC_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample,
                  show_visualizations = FALSE)

# ANOVA comparisons
anova(Model08, Model05)
anova(Model08, Model06)
anova(Model08, Model07)

#### 7. Region-wise regressions: Plasma GFAP ~ Aß PET (Binary MA) ####
NAV_brain_regions <- c(
     "NAV_HC_z", "NAV_ENTH_z", "NAV_AG_z", "NAV_THAL_z", "NAV_CAUD_z", "NAV_PUTA_z", "NAV_PALL_z",
     "NAV_ACUM_z", "NAV_cACC_z", "NAV_cMFC_z", "NAV_CUN_z", "NAV_FG_z", "NAV_IP_z", "NAV_infTC_z",
     "NAV_isthCC_z", "NAV_OC_z", "NAV_latOFC_z", "NAV_LING_z", "NAV_medOFC_z", "NAV_midTC_z",
     "NAV_paraHC_z", "NAV_paraCent_z", "NAV_parsOPC_z", "NAV_parsORB_z", "NAV_parsTRI_z",
     "NAV_periCalc_z", "NAV_postCent_z", "NAV_PCC_z", "NAV_preCent_z", "NAV_PCUN_z", "NAV_rACC_z",
     "NAV_rmidFC_z", "NAV_SF_z", "NAV_SP_z", "NAV_supTC_z", "NAV_supraMG_z", "NAV_transTC_z", "NAV_insula_z"
)

##### MA negative ####
results_table_NAVregional_MAneg <- run_region_regressions(data = TSPO_by_MA[["MA-"]], 
                                                          regions = NAV_brain_regions,
                                                          outcome_y = "plasmaGFAPpgmL_normalized_z")
print(results_table_NAVregional_MAneg)

##### MA positive ####
results_table_NAVregional_MApos <- run_region_regressions(data = TSPO_by_MA[["MA+"]], 
                                                          regions = NAV_brain_regions,
                                                          outcome_y = "plasmaGFAPpgmL_normalized_z")
print(results_table_NAVregional_MApos)

#### 8. Region-wise regressions: Plasma GFAP ~ global Aß PET x TSPO PET #### 
PBR_brain_regions <- c(
     "PBR_HC_z", "PBR_ENTH_z", "PBR_AG_z", "PBR_THAL_z", "PBR_CAUD_z", "PBR_PUTA_z", "PBR_PALL_z",
     "PBR_ACUM_z", "PBR_cACC_z", "PBR_cMFC_z", "PBR_CUN_z", "PBR_FG_z", "PBR_IP_z", "PBR_infTC_z",
     "PBR_isthCC_z", "PBR_OC_z", "PBR_latOFC_z", "PBR_LING_z", "PBR_medOFC_z", "PBR_midTC_z",
     "PBR_paraHC_z", "PBR_paraCent_z", "PBR_parsOPC_z", "PBR_parsORB_z", "PBR_parsTRI_z",
     "PBR_periCalc_z", "PBR_postCent_z", "PBR_PCC_z", "PBR_preCent_z", "PBR_PCUN_z", "PBR_rACC_z",
     "PBR_rmidFC_z", "PBR_SF_z", "PBR_SP_z", "PBR_supTC_z", "PBR_supraMG_z", "PBR_transTC_z", "PBR_insula_z"
)
results_table_regionalPBR <- run_region_regressions(data = TSPO_PrimarySample,
                                                    regions = PBR_brain_regions,
                                                    interaction_var = "NeoctxAZD_SUVR_z",
                                                    outcome_y = "plasmaGFAPpgmL_normalized_z")
print(results_table_regionalPBR)

#### 9. Region-wise regressions: Tau PET ~ plasma GFAP (Binary MA) #### 
outcome_brain_regions_MK <- c(
     "MK_HC_z", "MK_ENTH_z", "MK_AG_z", "MK_THAL_z", "MK_CAUD_z", "MK_PUTA_z", "MK_PALL_z",
     "MK_ACUM_z", "MK_cACC_z", "MK_cMFC_z", "MK_CUN_z", "MK_FG_z", "MK_IP_z", "MK_infTC_z",
     "MK_isthCC_z", "MK_OC_z", "MK_latOFC_z", "MK_LING_z", "MK_medOFC_z", "MK_midTC_z",
     "MK_paraHC_z", "MK_paraCent_z", "MK_parsOPC_z", "MK_parsORB_z", "MK_parsTRI_z",
     "MK_periCalc_z", "MK_postCent_z", "MK_PCC_z", "MK_preCent_z", "MK_PCUN_z", "MK_rACC_z",
     "MK_rmidFC_z", "MK_SF_z", "MK_SP_z", "MK_supTC_z", "MK_supraMG_z", "MK_transTC_z", "MK_insula_z"
)

##### MA negative ####
results_table_MAneg_MKregional <- run_region_regressions(data = TSPO_by_MA[["MA-"]],
                                                         regions = outcome_brain_regions_MK, 
                                                         outcome_y = "plasmaGFAPpgmL_normalized_z")
print(results_table_MAneg_MKregional)

##### MA positive ####
results_table_MApos_MKregional <- run_region_regressions(data = TSPO_by_MA[["MA+"]],
                                                         regions = outcome_brain_regions_MK, 
                                                         outcome_y = "plasmaGFAPpgmL_normalized_z")
print(results_table_MApos_MKregional)

#### 10. Sensitivity/Exploratory Analyses #### 

##### Other TSPO summary measures ####

###### Stratification #### 

# Difference composite
TSPO_PrimarySample_MA_DKT_difference_composite <-  TSPO_PrimarySample %>%
     split(.$MA_positivity_DKT_difference_composite)

# T composite
TSPO_PrimarySample_MA_DKT_T_composite <-  TSPO_PrimarySample %>%
     split(.$MA_positivity_DKT_T_composite)

###### Scatter - Models for Difference Composite  ####
scatter02 <- scatter_plot(
     data = TSPO_PrimarySample,
     group_var = "MA_positivity_DKT_difference_composite"
)
print(scatter02)

###### Model09 ####
Model09  <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + sex + DX2,
                   data = TSPO_PrimarySample_MA_DKT_difference_composite[["MA-"]],
                   show_visualizations = FALSE)

###### Model10 ####
Model10  <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + sex + DX2,
                   data = TSPO_PrimarySample_MA_DKT_difference_composite[["MA+"]],
                   show_visualizations = FALSE)

###### Model11 ####
Model11  <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z * MA_positivity_DKT_difference_composite +
                        age_at_mri + sex + DX2,
                   data = TSPO_PrimarySample,
                   show_visualizations = FALSE)

###### Model12 ####
Model12  <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z * PBR_DKT_difference_composite_z +
                        age_at_mri + sex + DX2,
                   data = TSPO_PrimarySample,
                   show_visualizations = FALSE,
                   visreg_args = list(
                        xvar = "NeoctxAZD_SUVR_z",
                        yvar = "PBR_DKT_difference_composite_z",
                        nn = 350, plot.type = "image",
                        col = colorRampPalette(
                             c("#B28C36","#DFCBA3","#F3F5F8","#C3CFE1",
                               "#94A4C0","#5F7498","#3B485F","#2A3345","#070A0D")
                        )(200)
                   )
)

###### Scatter - Models for T Composite  ####
scatter03 <- scatter_plot(
     data = TSPO_PrimarySample,
     group_var = "MA_positivity_DKT_T_composite"
)
print(scatter03)

###### Model13 ####
Model13 <- run_lm(formula =  plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_MA_DKT_T_composite[["MA-"]],
                  show_visualizations = FALSE)

###### Model14 ####
Model14 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_MA_DKT_T_composite[["MA+"]],
                  show_visualizations = FALSE)

###### Model15 ####
Model15 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z * MA_positivity_DKT_T_composite +
                       age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample,
                  show_visualizations = FALSE)

###### Model16 ####
Model16 <- run_lm(plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z * PBR_DKT_T_composite_z +
                       age_at_mri + sex + DX2,
                  TSPO_PrimarySample,
                  show_visualizations = FALSE,
                  visreg_args = list(
                       xvar = "NeoctxAZD_SUVR_z",
                       yvar = "PBR_DKT_T_composite_z",
                       nn = 350, plot.type = "image",
                       col = colorRampPalette(
                            c("#B28C36","#DFCBA3","#F3F5F8","#C3CFE1",
                              "#94A4C0","#5F7498","#3B485F","#2A3345","#070A0D")
                       )(200)
                  )
)

##### Sex Effects ####

###### Model17 ####
Model17 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z * sex + age_at_mri + DX2,
                  data = TSPO_by_MA[["MA-"]],
                  show_visualizations = FALSE)

###### Model18 ####
Model18 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z * sex + age_at_mri + DX2,
                  data = TSPO_by_MA[["MA+"]],
                  show_visualizations = FALSE)

###### Stratified by Sex within MA groups ######
TSPO_PrimarySample_MAneg_SexSplit <-  TSPO_by_MA[["MA-"]] %>%
     split(.$sex)
TSPO_PrimarySample_MApos_SexSplit <-  TSPO_by_MA[["MA+"]] %>%
     split(.$sex)

###### Model19 ####
Model19 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + DX2,
                  data = TSPO_PrimarySample_MAneg_SexSplit[["F"]],
                  show_visualizations = FALSE)

###### Model20 ####
Model20 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + DX2,
                  data = TSPO_PrimarySample_MAneg_SexSplit[["M"]],
                  show_visualizations = FALSE)

###### Model21 ####
Model21 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + DX2,
                  data = TSPO_PrimarySample_MApos_SexSplit[["F"]],
                  show_visualizations = FALSE)

###### Model22 ####
Model22 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + DX2,
                  data = TSPO_PrimarySample_MApos_SexSplit[["M"]],
                  show_visualizations = FALSE)

##### Including Outliers ####

###### Stratification #### 
TSPO_PrimarySample_Outliers <- read_excel(
     "Input_Tables/TSPO_PrimarySample_Outliers.xlsx"
) %>%
     mutate(
          sex = factor(sex, levels = c("F", "M")),
          DX2 = factor(DX2, levels = c("CU", "CI")),
          MA_positivity = factor(MA_positivity, levels = c("MA-", "MA+"))
     )
TSPO_PrimarySample_Outliers_MA <- TSPO_PrimarySample_Outliers %>%
     split(.$MA_positivity)

###### Model23 ####
Model23 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_Outliers_MA[["MA-"]],
                  show_visualizations = FALSE)

###### Model24 ####
Model24 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_Outliers_MA[["MA+"]],
                  show_visualizations = FALSE)

###### Model25 ####
Model25 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z * MA_positivity + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_Outliers,
                  show_visualizations = FALSE)

###### Model26 ####
Model26 <- run_lm(formula = plasmaGFAPpgmL_normalized_z ~ NeoctxAZD_SUVR_z * PBR_PCC_z + age_at_mri + sex + DX2,
                  data = TSPO_PrimarySample_Outliers,
                  show_visualizations = FALSE)
