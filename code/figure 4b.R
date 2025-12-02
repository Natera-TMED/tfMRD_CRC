library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(survminer)
library(survival)

# load the data
df<- readRDS("~/Desktop/project/tfMRD/act_data.rds")
surv_fit <- survfit(Surv(time, status) ~ category_merged, data = df)

ggsurv <- ggsurvplot(
  fit = surv_fit,
  data = df,
  conf.int = FALSE,
  pval = FALSE, 
  risk.table = TRUE,
  break.time.by = 6,
  xlab = "Months (DFS from ACT end)",
  ylab = "Disease-free survival probability",
  legend.title = "Category"
)
cox_model <- coxph(Surv(time, status) ~ category_merged, data = df)
sm <- summary(cox_model)

hr_tbl <- tibble(
  term = rownames(sm$coefficients),
  HR   = sm$coefficients[, "exp(coef)"],
  LCL  = sm$conf.int[, "lower .95"],
  UCL  = sm$conf.int[, "upper .95"],
  P    = sm$coefficients[, "Pr(>|z|)"]
) %>%
  mutate(
    cmp = gsub("^category", "", term),
    HR  = ifelse(is.finite(HR), round(HR, 2), NA_real_),
    LCL = ifelse(is.finite(LCL), round(LCL, 2), NA_real_),
    UCL = ifelse(is.finite(UCL), round(UCL, 2), NA_real_),
    P   = format.pval(P, digits = 3)
  )

ref_level <- levels(df$category_merged)[1]
hr_block <- paste0(
  "Cox HRs vs ",
  paste0(hr_tbl$cmp, ": ", hr_tbl$HR, " (", hr_tbl$LCL, "–",hr_tbl$UCL, "), p=", hr_tbl$P), collapse = "\n"
)


x_pos <- max(1, 0.05 * max(df$time, na.rm = TRUE))
ggsurv$plot <- ggsurv$plot +
  annotate("text", x = x_pos, y = 0.2, hjust = 0,
           label = hr_block, size = 3.8)

ggsurv