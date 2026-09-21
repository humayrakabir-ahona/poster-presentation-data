###############################################################################
# Kaplan-Meier plot: Ferroptosis-associated molecular state vs overall survival
# TCGA-LUAD, median-score stratification
#
# Requires: survival, survminer, ggplot2
#   install.packages(c("survival", "survminer", "ggplot2"))
###############################################################################

library(survival)
library(survminer)
library(ggplot2)

## ---------------------------------------------------------------------------
## 1. DATA
## ---------------------------------------------------------------------------
# Replace this block with your real data. You need a data.frame with:
#   time   : overall survival time, in days (numeric)
#   status : event indicator (1 = death/event, 0 = censored)
#   group  : factor with levels "High", "Low" (ferroptosis-associated score,
#            split at the median score)
#
# e.g. df <- read.csv("your_clinical_data.csv")
#      df$group <- factor(ifelse(df$score > median(df$score), "High", "Low"),
#                          levels = c("High", "Low"))

set.seed(42)
n_high <- 252
n_low  <- 251

# Simulated example data only — calibrated so HR ~ 0.80 and log-rank P ~ 0.13,
# matching the reported statistics. DELETE this simulation block once you
# plug in real data above.
sim_group <- function(n, group_label, rate, max_followup = 7300) {
  event_time <- rexp(n, rate = rate)
  censor_time <- runif(n, min = 200, max = max_followup)
  time <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)
  data.frame(time = time, status = status, group = group_label)
}

df <- rbind(
  sim_group(n_high, "High", rate = 0.00035),
  sim_group(n_low,  "Low",  rate = 0.00035 / 0.80)
)
df$group <- factor(df$group, levels = c("High", "Low"))

## ---------------------------------------------------------------------------
## 2. FIT
## ---------------------------------------------------------------------------
fit <- survfit(Surv(time, status) ~ group, data = df)

cox <- coxph(Surv(time, status) ~ group, data = df)
cox_summary <- summary(cox)

hr     <- cox_summary$conf.int[1, "exp(coef)"]
hr_lo  <- cox_summary$conf.int[1, "lower .95"]
hr_hi  <- cox_summary$conf.int[1, "upper .95"]
logrank_p <- survdiff(Surv(time, status) ~ group, data = df)$pvalue

# Labels shown in the plot legend
n_labels <- table(df$group)
group_labels <- c(
  High = paste0("High ferroptosis-associated score (n=", n_labels["High"], ")"),
  Low  = paste0("Low ferroptosis-associated score (n=",  n_labels["Low"],  ")")
)

annotation_text <- paste0(
  "Median-score stratification\n",
  "Log-rank P = ", formatC(logrank_p, digits = 3, format = "f"), "\n",
  "Cox HR = ", formatC(hr, digits = 2, format = "f"),
  " (95% CI ", formatC(hr_lo, digits = 2, format = "f"),
  "\u2013", formatC(hr_hi, digits = 2, format = "f"), ")"
)

## ---------------------------------------------------------------------------
## 3. PLOT
## ---------------------------------------------------------------------------
palette_colors <- c(High = "#1F4E79", Low = "#B03A2E")  # navy / brick red

p <- ggsurvplot(
  fit,
  data              = df,
  palette           = unname(palette_colors[levels(df$group)]),
  legend            = c(0.30, 0.14),          # inside plot, lower-left
  legend.title      = "",
  legend.labs       = group_labels,
  censor.shape       = "+",
  censor.size        = 4,
  size                = 1,
  xlab               = "Overall survival time (days)",
  ylab               = "Overall survival probability",
  break.x.by         = 1000,
  ylim               = c(0, 1),
  ggtheme            = theme_minimal(base_size = 16),
  font.x             = c(20, "plain", "black"),
  font.y             = c(20, "plain", "black"),
  font.tickslab      = c(16, "plain", "black"),
  font.legend        = c(16, "plain", "black")
)

p$plot <- p$plot +
  labs(title = "Ferroptosis-associated molecular state and overall survival \u2014 TCGA-LUAD") +
  theme(
    plot.title       = element_text(face = "bold", size = 22, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88"),
    axis.line        = element_line(color = "black"),
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(t = 15, r = 20, b = 40, l = 15)
  ) +
  annotate(
    "label",
    x = max(df$time) * 0.62, y = 0.97,
    label = annotation_text,
    hjust = 0, vjust = 1,
    size = 5, fill = "white", color = "black",
    label.size = 0.4, label.r = unit(0.1, "lines")
  ) +
  labs(caption = "Exploratory public-dataset association \u2014 not evidence of therapeutic efficacy") +
  theme(
    plot.caption = element_text(face = "bold", size = 14, hjust = 0.5,
                                 margin = margin(t = 15))
  )

## ---------------------------------------------------------------------------
## 4. SAVE
## ---------------------------------------------------------------------------
ggsave("Figure_TCGA_LUAD_ferroptosis_score_KM.png",
       plot = p$plot, width = 13, height = 8, dpi = 300, bg = "white")

print(p)
