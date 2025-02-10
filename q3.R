library(readxl)
library(data.table)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)

options(width = 600)
plot_theme <- theme(
  axis.line = element_line(color = "black"),
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
)

plot_dc_over_time_by_sbj_per_marker <- function(plot_data) {
  marker <- unique(plot_data$marker)
  if (length(marker) > 1) {
    stop("plot_dc_over_time_by_sbj_per_marker function works on one marker.")
  }
  # Use this plot to assess random effect from subjects
  m <- ggplot(plot_data, aes(
    x=sample_time, y=log(digital_count), color=trt_group, group=sbj_id)
  ) +
  geom_point(size=2) +
  geom_smooth(method="glm", linetype="dashed", se=FALSE) +
  facet_wrap(~ sbj_id, nrow=3, scales="free_y") +
  theme_bw() + plot_theme
  # dump to the current wkdir.
  out <- file.path(plot_dir, paste("dc_over_time", marker, "pdf", sep="."))
  ggsave(out, m, height=10, width=18)
}

fit_per_marker <- function(marker_data) {
  # make sure to fit model per marker
  if (length(unique(marker_data$marker)) > 1) {
    stop(
      "fit_per_marker function fits a model per marker. Found more than one in the give data."
    )
  }
  
  # fit liner mixed effect model with fixed effects sample_time and trt_group,
  # and subject as random effect.
  # we can't model the data with random effect of slope (change of digital
  # count over time) as we do not have enough data to do so.
  # From the dc_over_time plots for each marker, it appears that different
  # subjects have not the same slope.
  fit <- lmer(
    digital_count ~ sample_time + trt_group + (1|sbj_id),
    data=marker_data
  )
  # check normality
  normality_test <- shapiro.test(residuals(fit))
  if (normality_test$p.value <= 0.05) {
      # use log transformation if original normality test fails.
      marker_data[, digital_count := log(digital_count + 1)]
      fit <- lmer(
        digital_count ~ sample_time + trt_group + (1|sbj_id),
        data=marker_data
      )
      normality_test <- shapiro.test(residuals(fit))
      # I kill the program here for manual intervention to check on alternative
      # way to transform the data. Or perhaps we could try non-linear models.
      if (normality_test$p.value <= 0.05) {
        print("Normality assumption still violates after trying log transformation.")
        stop("Review the distribution of digital_count and figure out alternative transformation method.")
      }
  }
  # we want to see if there is a difference in change of digital count over
  # time among trt_group. Adding interaction between sample_time and trt_group.
  # although the question we are asking is not about this.
  fit_interaction <- lmer(
    digital_count ~ sample_time * trt_group + (1|sbj_id),
    data=marker_data
  )
  # compare model with and without interaction term.
  # return the model with interaction term if it fits better.
  model_comp <- anova(fit, fit_interaction)
  model_comp_pval <- model_comp["fit_interaction", "Pr(>Chisq)"]
  # for the given data, model with interaction does not fit better.
  # there is no significant difference in change of digital count over time
  # among treatment groups.
  if (model_comp_pval <= 0.05) {
    list(
      model=fit_interaction, normality_pval = normality_test$p.value
    )
  } else {
    list(
      model=fit, normality_pval = normality_test$p.value
    )
  }
}

# function to collect test results given model
test_per_marker <- function(model) {
  anova_res <- anova(model)
  # get p value for testing effects of treatment group.
  trt_group_effect_pvalue <- anova_res["trt_group", "Pr(>F)"]
  
  # compare avg difference between pairs of sample time given model.
  em_time <- emmeans(model, ~ sample_time)
  # names(em_time_pairs_dt)
  # contrast, estimate, SE, df, t.ratio, p.value
  em_time_pairs_dt <- as.data.table(pairs(em_time, reverse=TRUE))
  # return trt_group_effect_pvalue and table of diff between pair of sample time.
  list(
    trt_group_effect_pvalue = trt_group_effect_pvalue,
    time_pair_diff = em_time_pairs_dt[, c("contrast", "p.value")]
  )
}

# Start here
plot_dir <- "./q3_plots"
if (! dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# assume the data file is under the same current folder as the Rscript
data <- "./data.xlsx"
if (! file.exists(data)) {
  stop("Please make sure data.xlsx file is under the same folder as this script")
}
dt <- setDT(read_excel(data, sheet="data"))
# make sample_time and trt_group factors
dt[, ":="(
  sample_time = factor(
  sample_time,
  levels=c("DAY1", "DAY8", "DAY15", "DAY22", "DAY29")
  ),
  trt_group = factor(trt_group, levels=c("TA", "TB", "TC"))
)]

# visualize digital counts over treatment group for each marker
# look for effect of trt_group on digital counts for each marker
m <- ggplot(dt, aes(x=trt_group, y=log(digital_count), color=marker)) +
  geom_boxplot() +
  facet_wrap(~ marker, nrow=1) +
  labs(x="Treatment group", y="Digital count in log scale") +
  theme_bw() + plot_theme
ggsave(file.path(plot_dir, "dc_over_trt_group.pdf"), m, height=10, width=14)

# visualize digital counts over treatment group for each marker
m <- ggplot(dt, aes(x=sample_time, y=log(digital_count), color=marker)) +
  geom_boxplot() +
  facet_wrap(~ marker, nrow=1) +
  labs(x="Sample time", y="Digital count in log scale") +
  theme_bw() + plot_theme
ggsave(file.path(plot_dir, "dc_over_time.box.pdf"), m, height=10, width=14)

# holds for the final result
res_dt <- data.table()
markers <- unique(dt$marker)
# fit and get test result for each marker
for (i in seq_len(length(markers))) {
  data <- dt[marker==markers[i]]
  # visualize dc over time per subject per marker
  plot_dc_over_time_by_sbj_per_marker(data)
  # fit model
  fit_res <- fit_per_marker(data)
  # get test result from model
  test_res <- test_per_marker(fit_res$model)
  # we want to test difference between day22 and day8
  # this is hard-coded but can be easily extent to general comparisons
  day22_day8_diff_pvalue <- test_res$time_pair_diff[contrast=="DAY22 - DAY8"][["p.value"]]
  # collect and process results
  marker_res_dt <- data.table(
    "marker"=markers[i],
    "trt_group_effect_pvalue"=test_res$trt_group_effect_pvalue,
    "day22_day8_diff_pvalue"=day22_day8_diff_pvalue
  )
  marker_res_dt[, ":="(
    significance_trt_group_effect = ifelse(trt_group_effect_pvalue <= 0.05, "yes", "no"),
    significance_day22_day8_diff = ifelse(day22_day8_diff_pvalue <= 0.05, "yes", "no")
  )]
  res_dt <- rbind(res_dt, marker_res_dt)
}
# print out the final result
# names(res_dt)
# 1.marker
# 2. trt_group_effect_pvalue
# 3. day22_day8_diff_pvalue
# 4. significance_trt_group_effect: yes or no
# 5. significance_day22_day8_diff: yes or no
print(res_dt)
fwrite(res_dt, "q3.res.tsv", sep="\t", row.names=FALSE, quote=FALSE)


