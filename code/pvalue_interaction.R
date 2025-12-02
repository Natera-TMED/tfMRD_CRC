# ---- Original Script ----
library(survival)
#install.packages("survminer")
library(survminer)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("ggfortify")
library(ggfortify)
#install.packages(c("lubridate", "ggsurvfit", "gtsummary", "tidycmprsk"))
#install.packages("coxphf")
library(gtsummary)
library(tidycmprsk)
library(lubridate)
library(dplyr)
library(coxphf)
library(readxl)

rm(list=ls()) # clears all environments plots etc so starts the script clean
setwd("~/Documents/Clinical Data Mining/CRC/Galaxy/tfMRD/")
getwd()

mrd_outcome <- read_excel("tfMRD_Signatera_pval_interaction.xlsx")

# === FILTER TO ELIGIBLE PATIENTS ONLY ===
mrd_outcome <- mrd_outcome %>%
  filter(Eligible == TRUE)

km_mrd <- with(mrd_outcome[which(mrd_outcome$`Latitude_Signatera`!= "NA"),],
               Surv(mrd_outcome$DFS.months.Landmark, mrd_outcome$DFS.Event))
km_mrd_fit <- survfit(Surv(DFS.months.Landmark, DFS.Event) ~ Latitude_Signatera, data = mrd_outcome)
surv_pvalue(km_mrd_fit) #p=0.47

ggsurvplot(km_mrd_fit,
           conf.int = FALSE,
           pval = FALSE,
           risk.table = TRUE,
           size = 1,
           palette = c("#253494", "#238443", "orange", "purple"),
           xlab = "Months from Landmark: MRD timepoint",
           break.time.by = 3,
           ylab = "DFS Probability", 
           legend = "none",
           legend.title = "MRD timepoint\nLatitude/Signatera",
           legend.labs = c("Latitude (-)",
                           "Latitude (+)",
                           "Signatera (-)",
                           "Signaterra (+)"))
summary(km_mrd_fit)$table

# === TOGGLE LANDMARK ON/OFF ===
landmark_on <- FALSE  # Set to FALSE to use full cohort (no landmarking)

if (landmark_on) {
  time_var <- "DFS.month.Landmark"
  analysis_data <- mrd_outcome %>%
    filter(Latitude_Signatera != "NA" & DFS.month.Landmark.2mo > 0)
} else {
  time_var <- "DFS.month.updated"
  analysis_data <- mrd_outcome %>%
    filter(Latitude_Signatera != "NA")
}

surv_formula <- as.formula(paste0("Surv(", time_var, ", DFS.Event) ~ methylationclassification_Adjuvant"))
km_fit <- survfit(surv_formula, data = analysis_data)

cat("\n--- PFS/OS Probabilities at 12 and 24 Months ---\n")
print(summary(km_fit, times = c(24,30,36)))

cat("\n--- Median PFS/OS for Each Group ---\n")
print(summary(km_fit))

#Calculate %PFS/OS at 12 mo, 24 mo
summary(km_mrd_fit, times= c(24,30,36))
#Calculate Median PFS and median OS for both groups
survfit(Surv(time = mrd_outcome$DFS.months.Landmark, event = mrd_outcome$DFS.Event)~Latitude_Signatera1, data = mrd_outcome)







# ---- Cox Model with Latitude x Signatera Interaction and Forest Plot ----

# Create clean dataset including only necessary variables
mrd_outcome_clean <- mrd_outcome %>%
  filter(!is.na(DFS.months.Landmark), !is.na(DFS.Event)) %>%
  transmute(
    DFS.months.Landmark = DFS.months.Landmark,
    DFS.Event = DFS.Event,
    Latitude = factor(methylationclassification),
    Signatera = factor(Signatera)
  )

# Check level counts
cat("\n--- Level Check ---\n")
print(table(mrd_outcome_clean$Latitude))
print(table(mrd_outcome_clean$Signatera))

# Check if both variables have at least 2 levels
if (nlevels(mrd_outcome_clean$Latitude) >= 2 & nlevels(mrd_outcome_clean$Signatera) >= 2) {
  
  # Set reference levels explicitly
  mrd_outcome_clean$Latitude <- relevel(mrd_outcome_clean$Latitude, ref = "0")
  mrd_outcome_clean$Signatera <- relevel(mrd_outcome_clean$Signatera, ref = "0")
  
  # Fit Cox model *using this exact dataset*
  cox_interaction <- coxph(Surv(DFS.months.Landmark, DFS.Event) ~ Latitude * Signatera,
                           data = mrd_outcome_clean)
  
  # Print model summary
  cat("\n--- Cox Model with Interaction ---\n")
  print(summary(cox_interaction))
  
  # Extract and display interaction p-value
  interaction_p <- summary(cox_interaction)$coefficients["Latitude1:Signatera1", "Pr(>|z|)"]
  cat(paste0("\nInteraction p-value: ", signif(interaction_p, 3), "\n"))
  
  # Generate forest plot using the SAME data used in model fitting
  ggforest(cox_interaction, data = model.frame(cox_interaction))
  
  
} else {
  cat("⚠️ Cannot compute interaction: one of the factors has only one level.\n")
}

