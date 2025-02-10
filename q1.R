library(data.table)
library(ggplot2)
library(readxl)

plot_theme <- theme(
  axis.line = element_line(color = "black"),
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
)

data <- "./data.xlsx"
if (! file.exists(data)) {
  stop("Please make sure data.xlsx file is under the same folder as this script")
}
dt <- setDT(read_excel(data, sheet="data"))
# compute sample mean and standard deviation
dt[, ":="(
  avg_dc = mean(log(digital_count), na.rm=TRUE),
  sd_dc = sd(log(digital_count), na.rm=TRUE)
), by=list(trt_group, sample_time, marker)]
dt[, sample_time := factor(
  sample_time, levels=c("DAY1", "DAY8", "DAY15", "DAY22", "DAY29")
)]
grp_cols <- c("trt_group", "sample_time", "marker")
cols <- c(grp_cols, "avg_dc", "sd_dc")
# names(cols)
# 1. trt_group
# 2. sample_time
# 3. marker
# 4. avg_dc
# 5. sd_dc
dt <- unique(dt, by=grp_cols)[, ..cols]
fwrite(dt, "q1.res.tsv", sep="\t", row.names=FALSE, quote=FALSE)

m <- ggplot(dt, aes(x=sample_time, y=avg_dc, color=trt_group)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin = avg_dc - sd_dc, ymax = avg_dc + sd_dc), width = 0.2) +
  facet_wrap(~ marker + trt_group, ncol = 3) +
  labs(
    x="Sample time",
    y="Mean digital count in log scale plus/minus 1 unit sd",
    color="Treatment group"
  ) +
  theme_bw() + plot_theme
ggsave("q1.pdf", m, height=12, width=15)
