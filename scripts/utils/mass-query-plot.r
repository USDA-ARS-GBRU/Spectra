#!/usr/bin/env Rscript
#' Mass Query Plotting Script
#' Plots coverage of extreme k-mers along genomic sequences.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(optparse)
  library(readr)
  library(scales)
})

# --- Constants and Themes ---

BASE_THEME <- theme_bw() + theme(
  axis.title.x = element_text(size = 5),
  axis.title.y = element_text(size = 5),
  legend.key.size = unit(5, "pt"),
  legend.spacing = unit(-3, "pt"),
  axis.text.x = element_text(angle = 0, vjust = 0.5, size = 5),
  axis.text.y = element_text(size = 5),
  line = element_blank(),
  axis.ticks = element_line(),
  legend.title = element_blank()
)

# --- Helper Functions ---

create_coverage_plot <- function(values, legend = FALSE, axes = TRUE, y_max = NULL) {
  values <- values %>%
    group_by(Sequence, Start, End) %>%
    arrange(desc(Bin), .by_group = TRUE) %>%
    mutate(ymax = cumsum(Count), ymin = ymax - Count) %>%
    ungroup()

  p <- ggplot(values, aes(x = (Start + End) / 2, ymin = ymin, ymax = ymax, fill = Bin)) +
    geom_ribbon()
  
  p <- p + scale_x_continuous(
      limits = c(min(values$Start), max(values$End)),
      n.breaks = 10,
      expand = c(0, 0),
      labels = scales::scientific
    ) +
    xlab("Window Position (nucleotide)") +
    ylab("Counts") +
    BASE_THEME +
    theme(plot.margin = margin(t = 2.5, l = 2.5, b = 2.5, r = 2.5),
          legend.position = "bottom")

  if (!is.null(y_max)) {
    p <- p + scale_y_continuous(limits = c(0, y_max), expand = c(0, 0))
  } else {
    p <- p + scale_y_continuous(expand = c(0, 0))
  }

  if (!axes) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(t = 0, l = -2.7, b = -2.7, r = 0)
    )
  }

  if (!legend) {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}

# --- Main Logic ---

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input coverage tsv file", dest = "input_filename"),
  make_option(c("-w", "--window-size"), type = "character", default = NULL, help = "Subset range 'N,M'", dest = "window_range"),
  make_option(c("-s", "--sequence"), type = "character", default = NULL, help = "Subset sequences 'A,B,C' or regex", dest = "sequences"),
  make_option(c("-e", "--regex"), action = "store_true", default = FALSE, help = "Use regex for subsetting", dest = "use_regex"),
  make_option(c("-o", "--output-prefix"), type = "character", default = NULL, help = "Output filename prefix", dest = "output_prefix"),
  make_option(c("-f", "--output-format"), type = "character", default = "png", help = "Output format (png, svg)", dest = "output_format"),
  make_option(c("-r", "--resolution"), type = "numeric", default = 300, help = "DPI resolution", dest = "resolution"),
  make_option(c("-l", "--show-legend"), action = "store_true", default = FALSE, help = "Display legend", dest = "show_legend"),
  make_option(c("-x", "--scale"), type = "numeric", default = 1, help = "X-axis scale", dest = "x_scale"),
  make_option(c("-a", "--axes"), action = "store_false", default = TRUE, help = "Display axes text", dest = "show_axes"),
  make_option(c("-k", "--keep-scale"), action = "store_true", default = FALSE, help = "Keep absolute scale", dest = "keep_scale"),
  make_option(c("-y", "--y-max"), type = "numeric", default = NULL, help = "Fixed Y-axis limit", dest = "y_max"),
  make_option(c("-u", "--uniform-y"), action = "store_true", default = FALSE, help = "Uniform Y across all plots", dest = "use_uniform_y")
)

parser <- OptionParser(usage = "%prog -i coverage.tsv [options]", option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$input_filename)) {
  cat("Error: No input TSV specified. Use -h for help.\n")
  quit(status = 1)
}

if (is.null(opt$output_prefix)) {
  opt$output_prefix <- tools::file_path_sans_ext(basename(opt$input_filename))
}
if (opt$output_format == "svg") suppressPackageStartupMessages(library(svglite))

# Load data
values <- readr::read_tsv(opt$input_filename, show_col_types = FALSE)

# Filtering
if (!is.null(opt$sequences)) {
  if (opt$use_regex) {
    values <- values %>% filter(grepl(opt$sequences, Sequence))
  } else {
    values <- values %>% filter(Sequence %in% unlist(strsplit(opt$sequences, ",")))
  }
}
if (!is.null(opt$window_range)) {
  coords <- as.numeric(unlist(strsplit(opt$window_range, ",")))
  values <- values %>% filter(Start >= coords[1], End <= coords[2])
}

if (nrow(values) == 0) {
  cat("Error: No data remains after filtering.\n")
  quit(status = 1)
}

# Y-axis scaling
global_y_max <- opt$y_max
if (opt$use_uniform_y) {
  global_y_max <- values %>%
    group_by(Sequence, Start, End) %>%
    summarise(total = sum(Count), .groups = "drop") %>%
    pull(total) %>%
    max()
}

seq_names <- unique(values$Sequence)
for (seq in seq_names) {
  temp_values <- values %>% filter(Sequence == seq)
  
  if (opt$keep_scale) {
    temp_range <- (max(temp_values$End) - min(temp_values$Start) + 1) / (1e6 * opt$x_scale)
    temp_length <- temp_range + (if (opt$show_legend) 2 else 0.5)
    if (!opt$show_axes) temp_length <- temp_length - 0.5
  } else {
    temp_length <- 10
  }
  
  # Split into low and high bins (typically pct001-pct050 and pct051-pct100)
  values_low <- temp_values %>% filter(as.integer(sub("pct", "", Bin)) <= 50)
  values_high <- temp_values %>% filter(as.integer(sub("pct", "", Bin)) > 50)

  if (nrow(values_low) > 0) {
    p_low <- create_coverage_plot(values_low, legend = opt$show_legend, axes = opt$show_axes, y_max = global_y_max)
    out_low <- paste0(opt$output_prefix, "_", seq, "_low.", opt$output_format)
    ggsave(filename = out_low, device = opt$output_format, width = temp_length, height = 3,
           units = "in", dpi = opt$resolution, limitsize = FALSE)
  }

  if (nrow(values_high) > 0) {
    p_high <- create_coverage_plot(values_high, legend = opt$show_legend, axes = opt$show_axes, y_max = global_y_max)
    out_high <- paste0(opt$output_prefix, "_", seq, "_high.", opt$output_format)
    ggsave(filename = out_high, device = opt$output_format, width = temp_length, height = 3,
           units = "in", dpi = opt$resolution, limitsize = FALSE)
  }
}
