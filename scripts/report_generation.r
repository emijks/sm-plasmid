library(ggplot2)
library(dplyr)
library(knitr)

args <- commandArgs(trailingOnly = TRUE)
seq <- args[1]
ref <- args[2]
tmp <- args[3]


setwd(tmp)

flagstat <- paste(readLines(paste0(seq, "_flagstat.txt")), collapse = "\n")

depth_data <- read.table(file.path(paste0(seq, "_depth.txt")),
  header = FALSE,
  col.names = c("reference", "position", "depth"),
  sep = "\t")

stats <- depth_data %>%
  group_by(reference) %>%
  summarise(
    mean_depth = mean(depth),
    min_depth = min(depth),
    max_depth = max(depth),
    no_coverage = sum(depth == 0) 
  )

mean_depth <- mean(depth_data$depth)  

p_depth <- ggplot(depth_data, aes(x = position, y = depth)) +
  geom_line(color = "#00b7ff", size = 0.3) + 
  geom_hline(yintercept = mean_depth, linetype = "dashed", color = "red") +
  annotate("text", x = max(depth_data$position) * 0.95, y = mean_depth * 1.05,
           label = paste0("Mean: ", round(mean_depth, 2)), 
           color = "#740202", size = 3, hjust = 1) +
  labs(x = "Position (bp)", y = "Depth") + 
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"))

ggsave(paste0(seq, "_depth_plot.png"), p_depth, width = 8, height = 4)


report <- paste0(
  "---\n",
  "title: 'Report for sample ", seq, "'\n",
  "output: html_document\n",
  "---\n\n",
  "## Sample Information\n\n",
  "- SEQ: ", seq, "\n",
  "- REF: ", ref, "\n\n",
  "## Alignment Statistics\n\n",
  "```\n", flagstat, "\n```\n\n",
  "## Coverage Statistics\n\n",
  "```{r, echo=FALSE}\n",
  "library(knitr)\n",
  "kable(stats, caption = 'Coverage Statistics')\n",
  "```\n\n",
  "## Plots\n\n",
  "### Coverage Depth\n\n",
  "<img src='", seq, "_depth_plot.png' width='800' height='400'>\n\n",
  "### IGV Visualization\n\n",
  "<img src='", seq, ".igv.png' width='800' height='3500'>\n\n"
)


# report <- paste0(
#   "---\n",
#   "title: 'Report for sample ", seq, "'\n",
#   "output: html_document\n",
#   "---\n\n",
#   "## Sample Information\n\n",
#   "- SEQ: ", seq, "\n",
#   "- REF: ", ref, "\n\n",
#   "## Alignment Statistics\n\n",
#   "```\n", flagstat, "\n```\n\n",
#   "## Coverage Statistics\n\n",
#   "```{r, echo=FALSE}\n",
#   "library(knitr)\n",
#   "kable(stats, caption = 'Coverage Statistics')\n",
#   "```\n\n",
#   "## Plots\n\n",
#   "### Coverage Depth\n\n",
#   "<img src='", seq, "_depth_plot.png' width='800' height='400'>\n\n"
# )

writeLines(report, paste0(seq, "_report.Rmd"))

# html
rmarkdown::render(paste0(seq, "_report.Rmd"), output_file = paste0(seq, "_report.html"))
file.remove(paste0(seq, "_report.Rmd"))
cat("Report: ", paste0(seq, "_report.html"), "\n")