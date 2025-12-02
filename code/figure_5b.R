# load library
library(dplyr)
library(ggplot2)
# Figure 5a
plot_df <- leadtime_df %>%
  select(BundlingID, first_positive_month, relapse_month = months_from_surgery, PtsID) %>%
  pivot_longer(
    cols = c(first_positive_month, relapse_month),
    names_to = "Timepoint_Type",
    values_to = "Months"
  ) %>%
  mutate(
    Timepoint_Type = recode(
      Timepoint_Type,
      first_positive_month = "First Positive",
      relapse_month = "Relapse"
    )
  )

# Generate the plot
ggplot(plot_df,
       aes(
         x = Months,
         y = factor(PtsID),
         color = Timepoint_Type,
         shape = Timepoint_Type
       )) +
  geom_point(size = 3) +
  geom_line(aes(group = BundlingID), color = "gray50") +
  labs(
    title = "Time from Surgery: First Positive vs. Relapse",
    x = "Months from Surgery",
    y = "Patient ID",
    color = "Timepoint",
    shape = "Timepoint"
  ) +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  theme_minimal()
