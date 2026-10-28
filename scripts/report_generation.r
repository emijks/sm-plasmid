#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(knitr)
  library(rmarkdown)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: report_generation.r <seq> <ref> <tmp> [fastqc_r1_zip] [fastqc_r2_zip]")
}

seq <- args[1]
ref <- args[2]
tmp <- args[3]
fastqc_r1_zip <- if (length(args) >= 4) args[4] else NA_character_
fastqc_r2_zip <- if (length(args) >= 5) args[5] else NA_character_

setwd(tmp)

# -----------------------------
# Display names
# -----------------------------
strip_sample_name <- function(x) {
  x <- basename(x)
  x <- sub("\\.(fastq|fq)(\\.gz)?$", "", x, ignore.case = TRUE)
  x <- sub("\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("_L[0-9]+_R[12]$", "", x, ignore.case = TRUE)
  x <- sub("_R[12]$", "", x, ignore.case = TRUE)
  x
}

strip_ref_name <- function(x) {
  x <- basename(x)
  x <- sub("\\.(fa|fasta|fna)$", "", x, ignore.case = TRUE)
  x
}

seq_display <- strip_sample_name(seq)
ref_display <- strip_ref_name(ref)

# -----------------------------
# Formatting helpers
# -----------------------------
fmt_int <- function(x) {
  ifelse(is.na(x), "NA", formatC(round(x), format = "d", big.mark = ","))
}

fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits, big.mark = ","))
}

fmt_pct <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", paste0(formatC(x, format = "f", digits = digits), "%"))
}

read_or_empty <- function(path) {
  if (file.exists(path)) readLines(path, warn = FALSE) else character(0)
}

extract_first_number <- function(lines, pattern) {
  hit <- grep(pattern, lines, value = TRUE)
  if (length(hit) == 0) return(NA_real_)
  as.numeric(sub("^\\s*([0-9]+).*", "\\1", hit[1]))
}

extract_percent <- function(lines, pattern) {
  hit <- grep(pattern, lines, value = TRUE)
  if (length(hit) == 0) return(NA_real_)
  pct <- regmatches(hit[1], regexpr("[0-9]+\\.?[0-9]*(?=%)", hit[1], perl = TRUE))
  if (length(pct) == 0) return(NA_real_)
  as.numeric(pct)
}

# -----------------------------
# FastQC parsing from *_fastqc.zip
# -----------------------------
fastqc_zip_exists <- function(zip_file) {
  !is.na(zip_file) && nzchar(zip_file) && file.exists(zip_file)
}

zip_member <- function(zip_file, pattern) {
  listing <- unzip(zip_file, list = TRUE)$Name
  hit <- listing[grepl(pattern, listing)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_fastqc_summary <- function(zip_file, read_label) {
  if (!fastqc_zip_exists(zip_file)) {
    return(data.frame(Read = read_label, Module = "FastQC zip not found", Status = "NA"))
  }
  summary_member <- zip_member(zip_file, "summary.txt$")
  if (is.na(summary_member)) {
    return(data.frame(Read = read_label, Module = "summary.txt not found", Status = "NA"))
  }
  con <- unz(zip_file, summary_member)
  on.exit(close(con), add = TRUE)
  txt <- readLines(con, warn = FALSE)
  if (length(txt) == 0) {
    return(data.frame(Read = read_label, Module = "Empty FastQC summary", Status = "NA"))
  }
  parts <- strsplit(txt, "\t")
  data.frame(
    Read = read_label,
    Module = vapply(parts, function(x) ifelse(length(x) >= 2, x[2], NA_character_), character(1)),
    Status = vapply(parts, function(x) ifelse(length(x) >= 1, x[1], NA_character_), character(1)),
    stringsAsFactors = FALSE
  )
}

extract_fastqc_basic_stats <- function(zip_file, read_label) {
  out_empty <- data.frame(
    Read = character(0),
    Metric = character(0),
    Value = character(0),
    stringsAsFactors = FALSE
  )

  if (!fastqc_zip_exists(zip_file)) return(out_empty)

  data_member <- zip_member(zip_file, "fastqc_data.txt$")
  if (is.na(data_member)) return(out_empty)

  con <- unz(zip_file, data_member)
  on.exit(close(con), add = TRUE)
  txt <- readLines(con, warn = FALSE)

  start <- grep("^>>Basic Statistics", txt)
  end_all <- grep("^>>END_MODULE", txt)

  if (length(start) == 0 || length(end_all) == 0) return(out_empty)

  end <- end_all[end_all > start[1]][1]
  if (is.na(end) || end <= start[1] + 1) return(out_empty)

  block <- txt[(start[1] + 1):(end - 1)]
  block <- block[!grepl("^#", block)]
  block <- block[nzchar(block)]

  if (length(block) == 0) return(out_empty)

  parts <- strsplit(block, "\t", fixed = FALSE)
  metric <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
  value <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))

  keep <- !is.na(metric) & nzchar(metric)
  if (!any(keep)) return(out_empty)

  data.frame(
    Read = rep(read_label, sum(keep)),
    Metric = metric[keep],
    Value = value[keep],
    stringsAsFactors = FALSE
  )
}


extract_fastqc_plots <- function(zip_file, read_label, out_dir) {
  out_empty <- data.frame(
    Read = character(0),
    Plot = character(0),
    File = character(0),
    stringsAsFactors = FALSE
  )

  if (!fastqc_zip_exists(zip_file)) return(out_empty)

  listing <- unzip(zip_file, list = TRUE)$Name
  plot_members <- listing[grepl("/Images/.*\\.png$", listing, ignore.case = TRUE)]

  if (length(plot_members) == 0) return(out_empty)

  # имя образца в каталоге: графики разных образцов не должны пересекаться
  plot_dir <- file.path(out_dir, paste0(seq, "_fastqc_plots_", read_label))
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  rows <- lapply(plot_members, function(member) {
    original_name <- basename(member)
    clean_name <- paste0(read_label, "_", original_name)
    target <- file.path(plot_dir, clean_name)

    unzip(zip_file, files = member, exdir = plot_dir, junkpaths = TRUE, overwrite = TRUE)
    extracted <- file.path(plot_dir, original_name)

    if (file.exists(extracted) && normalizePath(extracted, mustWork = FALSE) != normalizePath(target, mustWork = FALSE)) {
      file.rename(extracted, target)
    }

    plot_label <- sub("\\.png$", "", original_name, ignore.case = TRUE)
    plot_label <- gsub("_", " ", plot_label)
    plot_label <- tools::toTitleCase(plot_label)

    data.frame(
      Read = read_label,
      Plot = plot_label,
      File = target,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(rows)
}

make_fastqc_plot_html <- function(plot_df) {
  if (nrow(plot_df) == 0) {
    return("<p><em>FastQC plots were not available.</em></p>\n")
  }

  html <- c()
  for (read_label in unique(plot_df$Read)) {
    html <- c(html, paste0("<h3>", read_label, " FastQC plots</h3>"))
    html <- c(html, "<div class='fastqc-grid'>")

    sub_df <- plot_df[plot_df$Read == read_label, , drop = FALSE]
    for (i in seq_len(nrow(sub_df))) {
      html <- c(
        html,
        paste0(
          "<div class='fastqc-card'>",
          "<div class='fastqc-title'>", sub_df$Plot[i], "</div>",
          "<img src='", sub_df$File[i], "'>",
          "</div>"
        )
      )
    }

    html <- c(html, "</div>")
  }

  paste(html, collapse = "\n")
}

fastqc_summary <- bind_rows(
  read_fastqc_summary(fastqc_r1_zip, "R1"),
  read_fastqc_summary(fastqc_r2_zip, "R2")
) %>%
  select(Read, Module, Status)

fastqc_basic <- bind_rows(
  extract_fastqc_basic_stats(fastqc_r1_zip, "R1"),
  extract_fastqc_basic_stats(fastqc_r2_zip, "R2")
)

fastqc_plots <- bind_rows(
  extract_fastqc_plots(fastqc_r1_zip, "R1", tmp),
  extract_fastqc_plots(fastqc_r2_zip, "R2", tmp)
)

fastqc_plot_html <- make_fastqc_plot_html(fastqc_plots)

# -----------------------------
# Input files generated by Snakemake rule report
# -----------------------------
flagstat_file <- paste0(seq, "_flagstat.txt")
depth_file <- paste0(seq, "_depth.txt")
depth_png <- paste0(seq, "_depth_plot.png")
igv_png <- paste0(seq, ".igv.png")

flagstat_lines <- read_or_empty(flagstat_file)
flagstat_text <- paste(flagstat_lines, collapse = "\n")

if (!file.exists(depth_file)) {
  stop("Depth file not found: ", depth_file)
}

depth_data <- read.table(
  depth_file,
  header = FALSE,
  col.names = c("reference", "position", "depth"),
  sep = "\t"
)

# -----------------------------
# Coverage metrics and E-score
# -----------------------------
coverage <- depth_data$depth
mean_depth <- mean(coverage)
median_depth <- median(coverage)
min_depth <- min(coverage)
max_depth <- max(coverage)
sd_depth <- sd(coverage)
cv_depth <- ifelse(mean_depth == 0, NA_real_, sd_depth / mean_depth)
zero_bases <- sum(coverage == 0)
breadth_1x <- mean(coverage >= 1) * 100
breadth_10x <- mean(coverage >= 10) * 100
breadth_50pct_mean <- mean(coverage >= 0.5 * mean_depth) * 100

e_score <- ifelse(mean_depth == 0, 0, mean(pmin(coverage / mean_depth, 1)))

# -----------------------------
# Flagstat metrics
# -----------------------------
total_reads <- extract_first_number(flagstat_lines, "in total")
primary_mapped <- extract_first_number(flagstat_lines, "primary mapped")
mapped_reads <- extract_first_number(flagstat_lines, " mapped \\(")
mapped_pct <- extract_percent(flagstat_lines, " mapped \\(")
properly_paired <- extract_first_number(flagstat_lines, "properly paired")
properly_paired_pct <- extract_percent(flagstat_lines, "properly paired")
duplicates <- extract_first_number(flagstat_lines, "duplicates")

# -----------------------------
# One non-duplicated metrics table
# -----------------------------
metrics_table <- data.frame(
  Metric = c(
    "Total reads",
    "Mapped reads",
    "Primary mapped reads",
    "Properly paired reads",
    "Duplicate reads",
    "Mean depth",
    "Median depth",
    "Minimum depth",
    "Maximum depth",
    "Coverage SD",
    "Coverage CV",
    "Zero-covered bases",
    "Breadth >=1x",
    "Breadth >=10x",
    "Breadth >=50% mean",
    "Coverage E-score"
  ),
  Value = c(
    fmt_int(total_reads),
    ifelse(is.na(mapped_pct), fmt_int(mapped_reads), paste0(fmt_int(mapped_reads), " (", fmt_pct(mapped_pct), ")")),
    fmt_int(primary_mapped),
    ifelse(is.na(properly_paired_pct), fmt_int(properly_paired), paste0(fmt_int(properly_paired), " (", fmt_pct(properly_paired_pct), ")")),
    fmt_int(duplicates),
    paste0(fmt_num(mean_depth, 1), "x"),
    paste0(fmt_num(median_depth, 1), "x"),
    paste0(fmt_int(min_depth), "x"),
    paste0(fmt_int(max_depth), "x"),
    paste0(fmt_num(sd_depth, 1), "x"),
    fmt_num(cv_depth, 3),
    fmt_int(zero_bases),
    fmt_pct(breadth_1x, 2),
    fmt_pct(breadth_10x, 2),
    fmt_pct(breadth_50pct_mean, 2),
    fmt_num(e_score, 3)
  ),
  Description = c(
    "All reads reported by samtools flagstat.",
    "Reads aligned to the reference.",
    "Primary mapped alignments only.",
    "Read pairs in expected orientation and distance.",
    "Reads marked as PCR or optical duplicates.",
    "Average per-base coverage.",
    "Median per-base coverage.",
    "Lowest per-base coverage.",
    "Highest per-base coverage.",
    "Standard deviation of per-base coverage.",
    "Coverage variability; lower is more uniform.",
    "Reference bases with depth 0.",
    "Reference bases covered by at least one read.",
    "Reference bases covered by at least 10 reads.",
    "Bases with coverage at least half of the mean.",
    "Coverage uniformity; 1.000 is perfectly even."
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# Coverage plot: compact, readable
# -----------------------------
plot_label <- paste0("Mean: ", fmt_num(mean_depth, 1), "x\nE-score: ", fmt_num(e_score, 3))

p_depth <- ggplot(depth_data, aes(x = position, y = depth)) +
  geom_line(color = "#2F6F9F", linewidth = 0.35) +
  geom_hline(yintercept = mean_depth, linetype = "dashed", color = "#B03A2E", linewidth = 1.15) +
  annotate(
    "label",
    x = max(depth_data$position, na.rm = TRUE) * 0.985,
    y = max(depth_data$depth, na.rm = TRUE) * 0.95,
    label = plot_label,
    hjust = 1,
    vjust = 1,
    size = 3.4,
    label.size = 0.25,
    fill = "white"
  ) +
  labs(x = "Position (bp)", y = "Depth") +
  theme_minimal(base_size = 12, base_family = "Arial") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#E5E7EB", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "#E5E7EB", linewidth = 0.3),
    axis.title = element_text(face = "bold", color = "#1F2933"),
    axis.text = element_text(color = "#374151")
  )

ggsave(depth_png, p_depth, width = 8.2, height = 2.9, dpi = 170, bg = "white")

# -----------------------------
# CSS + R Markdown report
# -----------------------------
css <- paste0(
  "<style>\n",
  "body { font-family: Arial, Helvetica, sans-serif; font-size: 14px; color: #1f2933; background: #ffffff; }\n",
  ".main-container { max-width: 1120px !important; }\n",
  "h1, h2, h3 { font-family: Arial, Helvetica, sans-serif; font-weight: 650; color: #102A43; }\n",
  "h1.title { font-size: 28px; margin-bottom: 8px; }\n",
  ".subtitle { font-size: 18px; margin: 2px 0 26px 0; color: #334E68; line-height: 1.45; }\n",
  ".subtitle strong { font-size: 22px; color: #102A43; }\n",
  ".section-note { color: #52606D; font-size: 13px; margin-top: -4px; margin-bottom: 12px; }\n",
  "table { font-family: Arial, Helvetica, sans-serif; font-size: 13px; width: 100%; }\n",
  "th { background-color: #EEF2F7; color: #102A43; font-weight: 650; }\n",
  "td, th { padding: 7px 9px !important; vertical-align: top; }\n",
  "caption { color: #52606D; font-size: 13px; padding-bottom: 6px; }\n",
  "pre, code { font-family: 'Menlo', 'Consolas', monospace; font-size: 12px; }\n",
  "pre { background-color: #F8FAFC; border: 1px solid #E5E7EB; border-radius: 6px; }\n",
  "#TOC { border: 1px solid #D9E2EC; border-radius: 12px; padding: 10px 8px; background: #F8FAFC; }\n",
  "#TOC:before { content: 'Contents'; display: block; font-weight: 700; color: #102A43; margin: 4px 8px 10px 8px; font-size: 15px; }\n",
  "#TOC ul { margin-left: 6px; }\n",
  "#TOC li { margin-bottom: 6px; line-height: 1.35; }\n",
  ".tocify-item { font-size: 13px; }\n",
  ".coverage-img { width: 72%; max-width: 820px; min-width: 560px; border: 1px solid #E5E7EB; border-radius: 8px; background: white; display: block; margin: 8px auto 20px auto; }\n",
  ".igv-img { width: 100%; max-width: 1120px; border: 1px solid #E5E7EB; border-radius: 8px; background: white; display: block; margin: 8px auto 20px auto; }\n",
  ".fastqc-status-PASS { color: #0B6B3A; font-weight: 700; }\n",
  ".fastqc-status-WARN { color: #9A6700; font-weight: 700; }\n",
  ".fastqc-status-FAIL { color: #B42318; font-weight: 700; }\n",
  ".fastqc-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); gap: 18px; margin: 12px 0 24px 0; }\n",
  ".fastqc-card { border: 1px solid #E5E7EB; border-radius: 10px; padding: 10px; background: #FFFFFF; }\n",
  ".fastqc-card img { width: 100%; height: auto; display: block; }\n",
  ".fastqc-title { font-weight: 650; color: #102A43; margin-bottom: 8px; font-size: 13px; }\n",
  "details { margin-top: 8px; margin-bottom: 20px; }\n",
  "summary { cursor: pointer; font-weight: 650; color: #334E68; margin-bottom: 8px; }\n",
  "hr { border: none; border-top: 1px solid #E5E7EB; margin: 24px 0; }\n",
  "</style>\n"
)

# Add CSS class to FastQC status values after kable rendering via HTML span.
fastqc_summary_display <- fastqc_summary %>%
  mutate(Status = paste0("<span class='fastqc-status-", Status, "'>", Status, "</span>"))

rmd_file <- paste0(seq, "_report.Rmd")
html_file <- paste0(seq, "_report.html")

report <- paste0(
  "---\n",
  "title: 'Plasmid alignment report'\n",
  "output:\n",
  "  html_document:\n",
  "    toc: true\n",
  "    toc_depth: 2\n",
  "    toc_float:\n",
  "      collapsed: false\n",
  "      smooth_scroll: true\n",
  "    number_sections: false\n",
  "    theme: cosmo\n",
  "    highlight: tango\n",
  "    df_print: paged\n",
  "---\n\n",
  css,
  "<div class='subtitle'>Sample: <strong>", seq_display, "</strong><br>Reference: <strong>", ref_display, "</strong></div>\n\n",

  "## Summary metrics\n\n",
  "```{r summary-metrics, echo=FALSE, message=FALSE, warning=FALSE}\n",
  "knitr::kable(metrics_table, caption = 'Main alignment and coverage metrics')\n",
  "```\n\n",

  "## Coverage profile\n\n",
  "<p class='section-note'>Per-base depth across the reference. The red dashed line marks the mean depth.</p>\n\n",
  "<img class='coverage-img' src='", depth_png, "'>\n\n",

  "## IGV alignment snapshot\n\n",
  "<p class='section-note'>Compact IGV snapshot of the BAM and VCF tracks.</p>\n\n",
  if (file.exists(igv_png)) {
    paste0("<img class='igv-img' src='", igv_png, "'>\n\n")
  } else {
    "<p><em>IGV snapshot not found.</em></p>\n\n"
  },

  "## FastQC metrics\n\n",
  "<p class='section-note'>FastQC summaries are embedded from the FastQC zip outputs, so the report remains readable after download.</p>\n\n",
  "### Basic statistics\n\n",
  "```{r fastqc-basic, echo=FALSE, message=FALSE, warning=FALSE}\n",
  "if (nrow(fastqc_basic) > 0) {\n",
  "  knitr::kable(fastqc_basic, caption = 'FastQC basic statistics')\n",
  "} else {\n",
  "  cat('<p><em>FastQC basic statistics were not available.</em></p>')\n",
  "}\n",
  "```\n\n",
  "### FastQC plots\n\n",
  fastqc_plot_html, "\n\n",

  "### Module status\n\n",
  "```{r fastqc-summary, echo=FALSE, message=FALSE, warning=FALSE, results='asis'}\n",
  "knitr::kable(fastqc_summary_display, caption = 'FastQC module status', escape = FALSE)\n",
  "```\n\n",

  "## Appendix\n\n",
  "<details>\n",
  "<summary>Complete samtools flagstat output</summary>\n\n",
  "```\n", flagstat_text, "\n```\n\n",
  "</details>\n\n"
)

writeLines(report, rmd_file, useBytes = TRUE)

rmarkdown::render(
  rmd_file,
  output_file = html_file,
  quiet = TRUE,
  envir = environment(),
  encoding = "UTF-8"
)

file.remove(rmd_file)
cat("Report: ", html_file, "\n")
