#### Load Library ####
library(readxl)
library(table1)
library(dplyr)
library(ggplot2)
library(grid)

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

# Run regressions programmatically
run_multiple_regressions <- function(outcomes, base_formula, data) {
     lapply(outcomes, function(outcome) {
          model <- lm(as.formula(paste(outcome, base_formula)), data = data)
          summary(model)
     }) |> setNames(outcomes)
}

# Reusable plot function
plot_grouped_bar <- function(data, var, ylab_text, ylim_range) {
     summary_stats <- data %>%
          group_by(Clinical_Group_CU_AD) %>%
          summarise(Mean = mean(.data[[var]], na.rm = TRUE),
                    SD = sd(.data[[var]], na.rm = TRUE), .groups = 'drop')
     
     ggplot(summary_stats, aes(x = Clinical_Group_CU_AD, y = Mean)) + 
          geom_col(fill = "gray", width = 0.7) +
          geom_errorbar(aes(ymin = Mean, ymax = Mean + SD), width = 0.2) +
          geom_point(data = data, aes(y = .data[[var]], color = Clinical_Group_CU_AD),
                     position = position_jitter(width = 0.2), size = 5.5) +
          theme(panel.grid.major.y = element_line(colour = "gray", size = 0.2),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black", size = 0.2),
                legend.position = "left") +
          ylab(ylab_text) + 
          xlab("Clinical Diagnosis") +
          coord_cartesian(ylim = ylim_range) +
          scale_y_continuous(breaks = seq(ylim_range[1], ylim_range[2], 0.2), expand = c(0,0))
}

#### 1. TSPO_signature_CUneg_ADpos ####

###### Load data ####
TSPO_signature_CUneg_ADpos <- read_excel("Input_Tables/TSPO_signature_CUneg_ADpos.xlsx")

# Apply factor levels
factor_levels <- list(
     sex = c("F", "M"),
     Clinical_Group_CU_AD = c("CU", "AD"),
     apoee4_status = c("Noncarrier", "Carrier"),
     AB_threshold = c("0", "1")
)
TSPO_signature_CUneg_ADpos <- set_factors(TSPO_signature_CUneg_ADpos, factor_levels)

##### Table ####

demographic_vars <- c("age_at_mri", "sex", "AB_threshold")
make_table1(vars = demographic_vars,
            group = "Clinical_Group_CU_AD",
            data = TSPO_signature_CUneg_ADpos,
            mean_sd = FALSE)

#### 2. TSPO regional patterns ####

# Collect PBR
continuous_vars <- c("PBR_HC", "PBR_ENTH", "PBR_AG", "PBR_THAL", "PBR_CAUD",
                     "PBR_PUTA", "PBR_PALL", "PBR_ACUM", "PBR_cACC", "PBR_cMFC",
                     "PBR_CUN", "PBR_FG", "PBR_IP", "PBR_infTC", "PBR_isthCC",
                     "PBR_OC", "PBR_latOFC", "PBR_LING", "PBR_medOFC", "PBR_midTC",
                     "PBR_paraHC", "PBR_paraCent", "PBR_parsOPC", "PBR_parsORB",
                     "PBR_parsTRI", "PBR_periCalc", "PBR_postCent", "PBR_PCC",
                     "PBR_preCent", "PBR_PCUN", "PBR_rACC", "PBR_rmidFC", "PBR_SF",
                     "PBR_SP", "PBR_supTC", "PBR_supraMG", "PBR_transTC", "PBR_insula")

# Table for mean (SD) per region
make_table1(vars = continuous_vars, group = "Clinical_Group_CU_AD", data = TSPO_signature_CUneg_ADpos, mean_sd = TRUE, overall = FALSE)

# Difference table between CU A- & AD A+
difference_table <- TSPO_signature_CUneg_ADpos %>%
     group_by(Clinical_Group_CU_AD) %>%
     summarise(across(all_of(continuous_vars), ~mean(.x, na.rm = TRUE), .names = "{.col}")) %>%
     tidyr::pivot_longer(-Clinical_Group_CU_AD, names_to = "Variable", values_to = "Mean") %>%
     tidyr::pivot_wider(names_from = Clinical_Group_CU_AD, values_from = Mean) %>%
     mutate(difference = AD - CU)
print(difference_table)

# Calculate the 90th percentile
print(quantile(difference_table$difference, probs = 0.9))

#### 3. Regressions ####

results_PBR <- run_multiple_regressions(
     continuous_vars, "~ Clinical_Group_CU_AD + age_at_mri + sex", TSPO_signature_CUneg_ADpos
)
print(results_PBR)

#### 4. Bar plots ####
TSPO_signature_CUneg_ADpos$Clinical_Group_CU_AD <- factor(ifelse(test = TSPO_signature_CUneg_ADpos$Clinical_Group_CU_AD == "CU", 
                                                                 yes = "CU AB-", no = "AD AB+"), levels = c("CU AB-", "AD AB+"))

plot_grouped_bar(data = TSPO_signature_CUneg_ADpos,
                 var = "PBR_PCC",
                 ylab_text = "TSPO PET SUVR (posterior cingulate cortex)",
                 ylim_range = c(0.8, 1.8))

plot_grouped_bar(data = TSPO_signature_CUneg_ADpos,
                 var = "PBR_DKT_difference_composite", 
                 ylab_text = "TSPO PET SUVR (difference-derived composite brain region)",
                 ylim_range = c(0.8, 1.8))

plot_grouped_bar(data = TSPO_signature_CUneg_ADpos,
                 var = "PBR_DKT_T_composite", 
                 ylab_text = "TSPO PET SUVR (T-derived composite brain region)",
                 ylim_range = c(0.8, 1.6))

#### 5. Group Comparisons ####

group_comp_results <- lapply(c("PBR_PCC", "PBR_DKT_difference_composite", "PBR_DKT_T_composite"), function(v) {
     aov(as.formula(paste(v, "~ Clinical_Group_CU_AD + age_at_mri + sex")), data = TSPO_signature_CUneg_ADpos) |> summary()
}) |> setNames(c("PCC", "Difference_composite", "T_composite"))

group_comp_results
