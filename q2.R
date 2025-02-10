library(data.table)
library(ggplot2)
library(readxl)

check_normality <- function(v) {
  # require a minimum of 3 samples.
  if(length(v) < 3) {
    return(FALSE)
  }
  # Given the limited number of measurement for each marker under
  # each treatment group, the test for normality can be not reliable.
  test_result <- shapiro.test(v)
  return(test_result$p.value > 0.05)
}

perform_test <- function(grp) {
  t1_dc <- log(grp$DAY1)
  t2_dc <- log(grp$DAY8)
  diffs <- t2_dc - t1_dc
  # Given the limited number of data, we can simply use
  # wilcoxon signed rank test.
  # To use wilcoxon test regardless of normality tests, set
  # is_normal <- FALSE
  is_normal <- check_normality(diffs)
  
  result <- list()
  if(is_normal) {
    # run parametric paired t test when normality assumption holds.
    t_test <- t.test(t2_dc, t1_dc, paired = TRUE)
    result$test_type <- "Paired t-test"
    result$pval <- t_test$p.value
  } else {
    # run non-parametric wilcoxon signed rank test when assumption is not met.
    wilcox_test <- wilcox.test(t2_dc, t1_dc, paired = TRUE)
    result$test_type <- "Wilcoxon signed-rank test"
    result$pval <- wilcox_test$p.value
  }
  
  return(result)
}

data <- "./data.xlsx"
if (! file.exists(data)) {
  stop("Please make sure data.xlsx file is under the same folder as this script")
}
dt <- setDT(read_excel(data, sheet="data"))
dt[, sample_time := factor(
  sample_time, levels=c("DAY1", "DAY8", "DAY15", "DAY22", "DAY29")
)]

# Reshape the data: long to wide
# marker, trt_group, subj_id, DAY1, DAY8 
sample_time_pairs <- c("DAY1", "DAY8")
paired_dt <- dcast(
  dt[sample_time %in% sample_time_pairs],
  trt_group + marker + sbj_id ~ sample_time,
  value.var = "digital_count"
)
paired_dt <- paired_dt[!is.na(DAY1) & !is.na(DAY8)]
# Run test per marker and trt_group combo
q2_res_dt <- paired_dt[,
  perform_test(.SD),
  by=list(marker, trt_group),
  .SDcols=c(sample_time_pairs, "sbj_id")
]
# correct for multi-testing
q2_res_dt[, pval_adj := p.adjust(pval, method="BH")]
q2_res_dt[, significant_diff := ifelse(pval_adj <= 0.05, "yes", "no")]
fwrite(q2_res_dt, "q2.res.tsv", sep="\t", row.names=FALSE, quote=FALSE)

