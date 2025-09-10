#Demographics Table
{```{r}
rm(list=ls())
setwd("~/Documents/Clinical Data Mining/CRC/Galaxy/tfMRD/")
circ_data <- read_excel("Demographics.xlsx")
circ_data <- circ_data[circ_data$Eligible=="TRUE",]

circ_data_subset <- circ_data %>%
  select(
    Age,
    Gender,
    ECOG,
    PrimSite,
    pT,
    pN,
    Stage,
    ACT,
    BRAF.V600E,
    RAS,
    MSI,
    DFS,
    OS.months)    %>%
  mutate(
    Age = as.numeric(Age),
    Gender = factor(Gender, levels = c("Male", "Female")),
    ECOG = factor(ECOG, levels = c(0, 1)),
    PrimSite = factor(PrimSite, levels = c("Right-sided colon", "Left-sided colon", "Rectum")),
    pT = factor(pT, levels = c("T1-T2", "T3-T4")),
    pN = factor(pN, levels = c("N0", "N1-N2")),
    Stage = factor(Stage, levels = c("I","II", "III", "IV")),
    ACT = factor(ACT, levels = c("TRUE", "FALSE"), labels = c("Adjuvant Chemotherapy", "Observation")),
    BRAF.V600E = factor(BRAF.V600E, levels = c("Negative", "Positive"), labels = c("BRAF wt", "BRAF V600E")),
    RAS = factor(RAS, levels = c("Negative", "Positive"), labels = c("RAS wt", "RAS mut")),
    MSI = factor(MSI, levels = c("Negative", "Positive"), labels = c("MSS", "MSI-H")),
    DFS = as.numeric(OS.months),
    OS.months = as.numeric(OS.months))
table1 <- circ_data_subset %>%
  tbl_summary(
    statistic = list(
      all_continuous() ~ "{median} ({min} - {max})",
      all_categorical() ~ "{n} ({p}%)")) %>%
  bold_labels()
table1
fit1 <- as_flex_table(
  table1,
  include = everything(),
  return_calls = FALSE)
fit1
save_as_docx(fit1, path= "~/Downloads/table1.docx")
```}
#Figure 2 & 3 KM
{rm(list = ls())
  # Figure 2 and 3 KM plot- adjust filter for surveillance group or MRD group

# Set working directory
setwd("~/Documents/Clinical Data Mining/CRC/Galaxy/tfMRD/")
getwd()

# Load necessary libraries
library(readxl)
library(survival)
library(survminer)
library(dplyr)

# ---- Section: Summary Statistics ----
# Read data
mrd_outcome <- read_excel("Full_cohort_surveillance_KM_pt_level.xlsx")

# === FILTER TO ELIGIBLE PATIENTS ONLY ===
mrd_outcome <- mrd_outcome %>%
  filter(Eligible == TRUE)

# === TOGGLE LANDMARK ON/OFF ===
landmark_on <- FALSE   # Set to FALSE to disable landmark analysis
landmark_weeks <- 10
landmark_months <- landmark_weeks / 4.35  # Convert to months


# Subset data: exclude rows with NA in methylationclassification
mrd_filtered <- mrd_outcome %>%
  filter(methylationclassification != "NA")

# Create DFS.landmarked.timepoint if using landmark
if (landmark_on) {
  mrd_landmark <- mrd_filtered %>%
    filter(DFS.months.updated > landmark_months) %>%
    mutate(DFS.landmarked.timepoint = DFS.months.updated - landmark_months)
  
  # Fit KM curve using landmarked data
  km_mrd_fit <- survfit(
    Surv(DFS.landmarked.timepoint, DFS.Event.updated) ~ methylationclassification,
    data = mrd_landmark
  )
  
  x_label <- paste0("Months from Landmark (", landmark_weeks, " weeks post-surgery)")
  
} else {
  # Standard full-cohort KM fit
  km_mrd_fit <- survfit(
    Surv(DFS.months.Landmark, DFS.Event.updated) ~ methylationclassification,
    data = mrd_filtered
  )
  
  x_label <- "Months from Landmark: 10 weeks post-surgery"
}

# Optional p-value display
surv_pvalue(km_mrd_fit)

# KM Plot
ggsurvplot(km_mrd_fit,
           conf.int = FALSE,
           pval = FALSE,
           risk.table = TRUE,
           size = 1,
           palette = c("#253494", "#238443"),
           xlab = x_label,
           break.time.by = 6,
           ylab = "DFS Probability",
           legend = "none",
           legend.title = " ",
           legend.labs = c("ctDNA (-)", "ctDNA (+)")) + ggtitle("Number at risk")

# View survival summary
summary(km_mrd_fit)$table

# ---- Section: Summary Statistics ----
if (landmark_on) {
  time_var <- "DFS.landmarked.timepoint"
  analysis_data <- mrd_landmark
} else {
  time_var <- "DFS.months.Landmark"
  analysis_data <- mrd_filtered
}

# Dynamic survival formula
surv_formula <- as.formula(paste0("Surv(", time_var, ", DFS.Event.updated) ~ methylationclassification"))
km_fit <- survfit(surv_formula, data = analysis_data)

# PFS/OS probabilities at 12 and 24 months
cat("\n--- PFS/OS Probabilities at 12 and 24 Months ---\n")
print(summary(km_fit, times = c(12, 24,36)))

# Median PFS/OS
cat("\n--- Median Survival Table ---\n")
library(survminer)

cat("\n--- Median Survival Times by Group ---\n")
print(surv_median(km_fit))
summary(km_mrd_fit, times= c(24,30,36))
summary(km_fit)$table


# ---- Section: Hazard Ratio (Cox model) ----
if (landmark_on) {
  mrd_subset <- mrd_landmark %>%
    filter(methylationclassification %in% c("0", "1"))
  
  surv_object <- Surv(time = mrd_subset$DFS.landmarked.timepoint, event = mrd_subset$DFS.Event.updated)
  
} else {
  mrd_subset <- mrd_filtered %>%
    filter(methylationclassification %in% c("0", "1"))
  
  surv_object <- Surv(time = mrd_subset$DFS.months.Landmark, event = mrd_subset$DFS.Event.updated)
}

# Relevel factor
mrd_subset$methylationclassification <- factor(mrd_subset$methylationclassification)
mrd_subset$methylationclassification <- relevel(mrd_subset$methylationclassification, ref = "0")

# Fit Cox model
cox_fit <- coxph(surv_object ~ methylationclassification, data = mrd_subset)

# View summary
summary(cox_fit)

# Forest plot of hazard ratio
ggforest(cox_fit, data = mrd_subset)
}
#Figure 2 & 3 Multivariate
{ #Figure 2 & 3 Mulivariate- adjust eligbility filter for each group
  library(dplyr)
  library(survival)
  library(survminer)
  library(coxphf)
  library(readxl)
  
  # best to have each factor as 2 variables but can use more (ie 0/1 or high/low etc)
  
  rm(list=ls()) # clears all environments plots etc so starts the script clean
  
  setwd("~/Documents/Clinical Data Mining/CRC/Galaxy/tfMRD/")
  multi_data <- read_excel("MRD_ACT_Benefit_KM.xlsx")
  multi_data <- subset(multi_data, !is.na(methylationclassification2))
  multi_data <- subset(multi_data, Eligible == TRUE)
  multi_data <- subset(multi_data, methylationclassification2 %in% c(0, 1))
  
  multi_datadf <- as.data.frame(multi_data)
  
  multi_datadf$Gender <- factor(multi_datadf$Gender, levels = c("Female", "Male"), labels = c("Female", "Male"))
  multi_datadf$Age <- factor(multi_datadf$Age.binary, levels = c("0", "1"), labels = c("<70", "≥70"))
  multi_datadf$Tumor_Location <- factor(multi_datadf$`Tumor.Location`, levels = c("Right-sided colon", "Left-sided colon", "Rectum"), labels = c("Right-sided colon", "Left-sided colon", "Rectum"))
  #multi_datadf$Stage <- relevel(factor(multi_datadf$Stage.group, levels = c("0", "1"), labels = c("II", "III")), ref = "II")
  multi_datadf$Stage <- factor(multi_datadf$Stage.group)
  multi_datadf$Stage <- factor(multi_datadf$Stage, levels = c("0", "1"), labels = c("II", "III"))
  multi_datadf$Stage <- relevel(multi_datadf$Stage, ref = "II")
  
  multi_datadf$MRD_Status <- factor(multi_datadf$methylationclassification2, levels = c("0", "1"), labels = c("ctDNA (-), obs", "ctDNA (-), ACT"))
  
  surv_object<-Surv(time = multi_datadf$DFS.month.Landmark, event = multi_datadf$DFS.Event) 
  cox_fit <- coxph(surv_object ~MRD_Status + Gender + Tumor_Location + Age + Stage, data=multi_datadf) 
  # Use this line instead of event only occurs in one of the factors 
  ggforest(cox_fit, data = multi_datadf, main = "Multivariate Regression Model for DFS", refLabel = "Reference Group")
  test.ph <- cox.zph(cox_fit)
  
}
#Figure 4
{##Figure 4##

#install.packages("survival")
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
setwd("~/Documents/Clinical\ Data\ Mining/CRC/Galaxy/tfMRD/")
getwd()


mrd_outcome <- read_excel("MRD_ACT_Benefit_KM.xlsx")

# === FILTER TO ELIGIBLE PATIENTS ONLY ===
# Make sure 'Eligible' column exists and is TRUE/FALSE or logical
mrd_outcome <- mrd_outcome %>%
  filter(Eligible == TRUE)

km_mrd <- with(mrd_outcome[which(mrd_outcome$methylationclassification_Adjuvant!= "NA"),], Surv(mrd_outcome$DFS.month.Landmark.2mo, mrd_outcome$DFS.Event))
km_mrd_fit <- survfit(Surv(DFS.month.Landmark.2mo, DFS.Event) ~ methylationclassification_Adjuvant, data = mrd_outcome)
surv_pvalue(km_mrd_fit) #p=0.47
#Plot version 2

ggsurvplot(km_mrd_fit,
           conf.int = FALSE,
           pval = FALSE,
           #fun = "cumhaz",
           risk.table = TRUE,
           size = 1,
           #linetype = "strata",
           palette = c("#8BC400","#46AC46", "#008BCE","#004677"
                       
           ),
           xlab = "Months from Landmark: 2 month post-surgery",
           break.time.by = 3,
           # xlim = c(0, 84),
           ylab = "DFS Probability", 
           legend = "none",
           legend.title = "MRD status",
           legend.labs = c("tfMRD (-), obs","tfMRD (-), ACT", "tfMRD (+), ACT", "tfMRD (+), obs"))
#surv.median.line = c("hv"))
summary(km_mrd_fit)$table



#color schemes
#253494 is blue
#238443 is green

# ---- Section: Load Data ----
# === TOGGLE LANDMARK ON/OFF ===
landmark_on <- TRUE  # Set to FALSE to use full cohort (no landmarking)

# === Define time variable and dataset dynamically ===
if (landmark_on) {
  time_var <- "DFS.month.Landmark"
  analysis_data <- mrd_outcome %>%
    filter(methylationclassification_Adjuvant != "NA" & DFS.month.Landmark.2mo > 0)  # Add check for landmarked patients
} else {
  time_var <- "DFS.month.updated"
  analysis_data <- mrd_outcome %>%
    filter(methylationclassification_Adjuvant != "NA")
}

# === Create Surv object dynamically and fit KM model ===
surv_formula <- as.formula(paste0("Surv(", time_var, ", DFS.Event) ~ methylationclassification_Adjuvant"))
km_fit <- survfit(surv_formula, data = analysis_data)

# === Calculate %PFS/OS at 12 and 24 months ===
cat("\n--- PFS/OS Probabilities at 12 and 24 Months ---\n")
print(summary(km_fit, times = c(24,30,36)))

# === Calculate Median PFS/OS for each group ===
cat("\n--- Median PFS/OS for Each Group ---\n")
print(summary(km_fit))
summary(km_fit)$table

#---adjusted HR calculation---

library(dplyr)
library(survival)
library(survminer)
library(coxphf)
library(readxl)


rm(list=ls()) # clears all environments plots etc so starts the script clean

setwd("~/Documents/Clinical Data Mining/CRC/Galaxy/tfMRD/")
multi_data <- read_excel("MRD_ACT_Benefit_KM.xlsx")
multi_data <- subset(multi_data, !is.na(methylationclassification2))
multi_data <- subset(multi_data, Eligible == TRUE)
multi_data <- subset(multi_data, methylationclassification2 %in% c(0, 1))

multi_datadf <- as.data.frame(multi_data)

multi_datadf$Gender <- factor(multi_datadf$Gender, levels = c("Female", "Male"), labels = c("Female", "Male"))
multi_datadf$Age <- factor(multi_datadf$Age.binary, levels = c("0", "1"), labels = c("<70", "≥70"))
multi_datadf$Tumor_Location <- factor(multi_datadf$`Tumor.Location`, levels = c("Right-sided colon", "Left-sided colon", "Rectum"), labels = c("Right-sided colon", "Left-sided colon", "Rectum"))

multi_datadf$MRD_Status <- factor(multi_datadf$methylationclassification2, levels = c("0", "1"), labels = c("ctDNA (-), obs", "ctDNA (-), ACT"))
# modify to fit adjusted HR grouping as needed

surv_object<-Surv(time = multi_datadf$DFS.month.Landmark, event = multi_datadf$DFS.Event) 
cox_fit <- coxph(surv_object ~MRD_Status + Gender + Tumor_Location + Age, data=multi_datadf) 

ggforest(cox_fit, data = multi_datadf, main = "Multivariate Regression Model for DFS", refLabel = "Reference Group")
test.ph <- cox.zph(cox_fit)

}