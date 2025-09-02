#### Load Libraries ####
library(readxl)
library(table1)
library(dplyr)
library(ggplot2)
library(visreg)
library(stats)
library(graphics)
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

# Scatterplot function
scatter_plot <- function(data, group_var, 
                         x = "centiloid", 
                         y = "plasmaGFAP_z",
                         xlab = "Aβ-PET centiloid", 
                         ylab = "Plasma GFAP (z-scored)",
                         colors = c("MA-" = "#B28C36", "MA+" = "#2A3345"),
                         xlim = c(-15, 155), 
                         ylim = c(-3, 6)) {
     
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

#### 1. sTREM2_SecondarySample_Combined ####

###### Load data ####
sTREM2_SecondarySample_Combined <- read_excel("Input_Tables/sTREM2_SecondarySample_Combined.xlsx")

# Apply factor levels
factor_levels <- list(
     cohort = c("TRIAD", "WRAP"),
     sex = c("F", "M"),
     APOE4_status = c("Noncarrier", "Carrier"),
     DX2 = c("CU", "CI"),
     MA_positivity = c("MA-", "MA+"),
     AB_positivity = c("0", "1")
)
sTREM2_SecondarySample_Combined <- set_factors(sTREM2_SecondarySample_Combined, factor_levels)
# Obs: variables "sTREM2_adj" (WRAP) and "sTREM2" (TRIAD) were both named as "sTREM2" for merging

##### Table ####
make_table1(vars = c("age", "sex", "education_years", "APOE4_status", 
                     "plasmaGFAP_z", "centiloid", "sTREM2_z", "DX2", "MMSE"),
            group = "MA_positivity", data = sTREM2_SecondarySample_Combined,
            mean_sd = FALSE, overall = TRUE)

#### 3. Stratification ####
sTREM2_SecondarySample_Combined_by_MA <- sTREM2_SecondarySample_Combined %>%
     split(.$MA_positivity)

#### 4. Scatterplot ####
scatter_plot(sTREM2_SecondarySample_Combined, group_var = "MA_positivity")

#### 5. Regressions ####

##### Model01 ####
Model01 <- run_lm(plasmaGFAP_z ~ centiloid + age + sex + DX2 + cohort,
                  data = sTREM2_SecondarySample_Combined_by_MA[["MA-"]],
                  show_visualizations = TRUE)

##### Model02 ####
Model02 <- run_lm(plasmaGFAP_z ~ centiloid + age + sex + DX2 + cohort,
                  data = sTREM2_SecondarySample_Combined_by_MA[["MA+"]],
                  show_visualizations = TRUE)

##### Model03 ####
Model03 <- run_lm(plasmaGFAP_z ~ centiloid*MA_positivity + age + sex + DX2 + cohort,
                  data = sTREM2_SecondarySample_Combined,
                  show_visualizations = TRUE)

##### Model04 ####
Model04 <- run_lm(plasmaGFAP_z ~ centiloid*sTREM2_z + age + sex + DX2 + cohort, 
                  data = sTREM2_SecondarySample_Combined,
                  show_visualizations = TRUE,
                  visreg_args = list(xvar = "centiloid",
                                     yvar = "sTREM2_z", nn = 350,
                                     plot.type = "image",
                                     col = colorRampPalette(c("#B28C36", "#DFCBA3", "#F3F5F8", 
                                                              "#C3CFE1","#94A4C0", "#5F7498", 
                                                              "#3B485F", "#2A3345", "#070A0D"))(200)))

##### 5.1 ANOVA ####

##### Model05 ####
Model05 <- run_lm(plasmaGFAP_z ~ centiloid + age + sex + DX2 + cohort,
                  data = sTREM2_SecondarySample_Combined,
                  show_visualizations = FALSE)

##### Model06 ####
Model06 <- run_lm(plasmaGFAP_z ~ sTREM2_z + age + sex + DX2 + cohort,
                  data = sTREM2_SecondarySample_Combined,
                  show_visualizations = FALSE)

##### Model07 ####
Model07 <- run_lm(plasmaGFAP_z ~ centiloid + sTREM2_z + age + sex + DX2 + cohort,
                  data = sTREM2_SecondarySample_Combined,
                  show_visualizations = FALSE)

##### Model08 ####
Model08 <- run_lm(plasmaGFAP_z ~ centiloid*sTREM2_z + age + sex + DX2 + cohort,
                  data = sTREM2_SecondarySample_Combined,
                  show_visualizations = FALSE)

anova(Model08, Model05)
anova(Model08, Model06)
anova(Model08, Model07)

#### 6. Sensitivity/Exploratory Analyses ####

# Stratify by MA and sex
sTREM2_MAneg_sex <- sTREM2_SecondarySample_Combined_by_MA[["MA-"]] %>%
     split(.$sex)
sTREM2_MApos_sex <- sTREM2_SecondarySample_Combined_by_MA[["MA+"]] %>%
     split(.$sex)

##### Model09 ####
Model09 <- run_lm(plasmaGFAP_z ~ centiloid*sex + age + DX2 + cohort,
                  data = sTREM2_SecondarySample_Combined_by_MA[["MA-"]],
                  show_visualizations = FALSE)

##### Model10 ####
Model10 <- run_lm(plasmaGFAP_z ~ centiloid + age + DX2 + cohort,
                  data = sTREM2_MAneg_sex[["F"]],
                  show_visualizations = FALSE)

##### Model11 ####
Model11 <- run_lm(plasmaGFAP_z ~ centiloid + age + DX2 + cohort,
                  data = sTREM2_MAneg_sex[["M"]],
                  show_visualizations = FALSE)

##### Model12 ####
Model12 <- run_lm(plasmaGFAP_z ~ centiloid*sex + age + DX2 + cohort,
                  data = sTREM2_SecondarySample_Combined_by_MA[["MA+"]],
                  show_visualizations = FALSE)

##### Model13 ####
Model13 <- run_lm(plasmaGFAP_z ~ centiloid + age + DX2 + cohort,
                  data = sTREM2_MApos_sex[["F"]],
                  show_visualizations = FALSE)

##### Model14 ####
Model14 <- run_lm(plasmaGFAP_z ~ centiloid + age + DX2 + cohort,
                  data = sTREM2_MApos_sex[["M"]],
                  show_visualizations = FALSE)
