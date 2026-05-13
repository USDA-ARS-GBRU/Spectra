#!/usr/bin/env Rscript
#' Spectra Visualization Script
#' Generates linear or circular plots of k-mer distributions along genomic sequences.

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

#' Build color palette for triplets
palette_builder <- function(triplet, palette = "base") {
  bases <- c("A", "C", "G", "T")
  colors <- switch(
    palette,
    "base"      = c("C6", "6C", "3C", "10"),
    "base2"     = c("BF", "6F", "46", "1F"),
    "base3"     = c("CF", "6F", "56", "2F"),
    "dual"      = c("96", "3C", "3C", "96"),
    "acontrast" = c("CC", "33", "33", "3F")
  )

  if (is.null(colors)) colors <- c("C6", "6C", "3C", "10") # Default to base

  # Construct color hex
  b1 <- substr(triplet, 1, 1)
  b2 <- substr(triplet, 2, 2)
  b3 <- substr(triplet, 3, 3)

  color <- paste0("#",
                  colors[match(b1, bases)],
                  colors[match(b2, bases)],
                  colors[match(b3, bases)])
  return(color)
}

#' Standard linear spectra plot
create_spectra_plot <- function(values, triplet_colors, legend = FALSE, facet = FALSE,
                               frequencies = FALSE, ylims = TRUE, scale = 1, axes = TRUE) {
  if (!frequencies) {
    values <- values %>% mutate(value = value / (End - Start + 1))
  }

  values <- values %>%
    group_by(Library, Sequence, Start, End) %>%
    arrange(desc(name), .by_group = TRUE) %>%
    mutate(ymax = cumsum(value), ymin = ymax - value) %>%
    ungroup()

  p <- ggplot(values, aes(x = (Start + End) / 2, ymin = ymin, ymax = ymax, fill = name)) +
    geom_ribbon()

  if (ylims) {
    p <- p + scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  } else {
    p <- p + scale_y_continuous(expand = c(0, 0))
  }

  xrange <- max(values$End) - min(values$Start) + 1
  if (log10(xrange) > log10(scale * 1e6) + 1) {
    scale <- scale * 10
  }

  p <- p + scale_fill_manual(values = triplet_colors) +
    scale_x_continuous(
      limits = c(min(values$Start), max(values$End)),
      n.breaks = 10,
      expand = c(0, 0),
      labels = scales::scientific
    ) +
    xlab("Window Position (nucleotide)") +
    ylab(if (frequencies) "Frequency" else "Proportion") +
    BASE_THEME +
    theme(plot.margin = margin(t = 2.5, l = 2.5, b = 2.5, r = 2.5))

  if (facet) {
    p <- p + facet_grid(rows = vars(Library))
  }
  if (!legend) {
    p <- p + theme(legend.position = "none")
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
  return(p)
}

#' Circular genome-wide plot
create_circular_plot <- function(values, triplet_colors, legend = FALSE, frequencies = FALSE,
                               ylims = TRUE, limit = 0) {
  suppressPackageStartupMessages(library(egg))

  # Definable sequence spacer size
  spacer_length <- 1000000
  spacer_values <- rbind(
    values %>% filter(End == min(End), Sequence == values$Sequence[1]) %>%
      mutate(Start = 1, End = 1000, value = 0, Sequence = "Spacer"),
    values %>% filter(End == min(End), Sequence == values$Sequence[1]) %>%
      mutate(Start = 1001, End = spacer_length - 1000, value = 0, Sequence = "Spacer"),
    values %>% filter(End == min(End), Sequence == values$Sequence[1]) %>%
      mutate(Start = spacer_length - 999, End = spacer_length, value = 0, Sequence = "Spacer")
  )

  # Modify sequence positions to be absolute
  seq_names <- unique(values$Sequence)
  new_values <- values %>% filter(Sequence == seq_names[1])
  current_max <- max(new_values$End)

  if (length(seq_names) > 1) {
    for (index in seq_names[2:length(seq_names)]) {
      new_values <- rbind(
        new_values,
        spacer_values %>% mutate(Start = Start + current_max, End = End + current_max)
      )
      current_max <- current_max + spacer_length
      current_values <- values %>% filter(Sequence == index)
      new_values <- rbind(
        new_values,
        current_values %>% mutate(Start = Start + current_max, End = End + current_max)
      )
      current_max <- max(new_values$End) + spacer_length
    }
  }

  if (frequencies) {
    p <- ggplot() + geom_area(data = new_values, aes(fill = name, x = (Start + End) / 2, y = value),
                             stat = "identity", position = "stack")
  } else {
    p <- ggplot() + geom_area(data = new_values, aes(fill = name, x = (Start + End) / 2, y = value / (End - Start + 1)),
                             stat = "identity", position = "stack")
  }

  if (ylims) {
    p <- p + scale_y_continuous(limits = c(0, 1), expand = c(0.5, 0.5))
  } else {
    p <- p + scale_y_continuous(expand = c(0.5, 0.5))
  }

  p <- p + scale_fill_manual(values = triplet_colors)

  plot_size <- if (limit > 0) limit + ((length(seq_names) - 1) * spacer_length) else current_max

  p <- p + scale_x_continuous(
    limits = c(1, plot_size),
    n.breaks = 30,
    expand = c(0, 0),
    labels = scales::scientific
  ) +
    xlab("Window Position (nucleotide)") +
    ylab("Proportion") +
    BASE_THEME +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(t = -270, l = -270, b = -270, r = -270, unit = "pt")
    )

  if (!legend) {
    p <- p + theme(legend.position = "none")
  }

  p <- p + coord_polar(start = 0)
  return(p)
}

# --- Main Logic ---

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "Input triplet tsv file", dest = "input_filename"),
  make_option(c("-w", "--window-size"), type = "character", default = NULL, help = "Subset range 'N,M'", dest = "window_range"),
  make_option(c("-s", "--sequence"), type = "character", default = NULL, help = "Subset sequences 'A,B,C' or regex", dest = "sequences"),
  make_option(c("-n", "--libraries"), type = "character", default = NULL, help = "Subset libraries 'A,B,C' or regex", dest = "libraries"),
  make_option(c("-e", "--regex"), action = "store_true", default = FALSE, help = "Use regex for subsetting", dest = "use_regex"),
  make_option(c("--graphlength"), type = "numeric", default = 0, help = "Max plot length for circular", dest = "graph_length"),
  make_option(c("--gff-file"), type = "character", default = NULL, help = "Overlapping GFF annotations", dest = "gff_file"),
  make_option(c("--gff-tracks"), type = "character", default = NULL, help = "GFF types 'A,B,C'", dest = "gff_tracks"),
  make_option(c("--trf-file"), type = "character", default = NULL, help = "Overlapping TRF-tsv annotations", dest = "trf_file"),
  make_option(c("-o", "--output-prefix"), type = "character", default = NULL, help = "Output filename prefix", dest = "output_prefix"),
  make_option(c("-f", "--output-format"), type = "character", default = "png", help = "Output image format (png, svg, etc)", dest = "output_format"),
  make_option(c("-r", "--resolution"), type = "numeric", default = 300, help = "DPI resolution", dest = "resolution"),
  make_option(c("-q", "--freq"), action = "store_true", default = FALSE, help = "Data is already frequencies", dest = "is_frequencies"),
  make_option(c("-l", "--show-legend"), action = "store_true", default = FALSE, help = "Display color legend", dest = "show_legend"),
  make_option(c("-y", "--ylims"), action = "store_false", default = TRUE, help = "Limit y-axis between 0,1", dest = "use_ylims"),
  make_option(c("-x", "--scale"), type = "numeric", default = 1, help = "X-axis scale (Mb per inch)", dest = "x_scale"),
  make_option(c("-a", "--axes"), action = "store_false", default = TRUE, help = "Display axes text", dest = "show_axes"),
  make_option(c("-k", "--keep-scale"), action = "store_true", default = FALSE, help = "Keep absolute scale", dest = "keep_scale"),
  make_option(c("--palette"), type = "character", default = "base", help = "Palette: base, dual", dest = "palette_type"),
  make_option(c("-c", "--circular"), action = "store_true", default = FALSE, help = "Circular plot", dest = "is_circular"),
  make_option(c("-t", "--transparent"), action = "store_true", default = FALSE, help = "Transparent background", dest = "is_transparent"),
  make_option(c("--ngaps"), type = "character", default = NULL, help = "GFF of N-gap coordinates", dest = "ngaps_file")
)

parser <- OptionParser(usage = "%prog -i triplet.tsv [options]", option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$input_filename)) {
  cat("Error: No input TSV specified. Use -h for help.\n")
  quit(status = 1)
}

# Determine script directory for resource files
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_all[grep(file_arg, args_all)])
if (length(script_path) == 0) script_path <- "."
script_dir <- dirname(script_path)

# Prepare output naming
if (is.null(opt$output_prefix)) {
  opt$output_prefix <- tools::file_path_sans_ext(basename(opt$input_filename))
}
if (opt$output_format == "svg") suppressPackageStartupMessages(library(svglite))

# Load data
values <- readr::read_tsv(opt$input_filename, show_col_types = FALSE)

# Filtering
if (!is.null(opt$libraries)) {
  if (opt$use_regex) {
    values <- values %>% filter(grepl(opt$libraries, Library))
  } else {
    values <- values %>% filter(Library %in% unlist(strsplit(opt$libraries, ",")))
  }
}
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

lib_names <- unique(values$Library)
seq_names <- unique(values$Sequence)

# Load auxiliary tracks
load_gff <- function(path, range_str) {
  if (is.null(path)) return(NULL)
  suppressPackageStartupMessages(library(ape))
  g <- read.gff(path)
  if (!is.null(range_str)) {
    crd <- as.numeric(unlist(strsplit(range_str, ",")))
    g <- g %>% filter(start >= crd[1], end <= crd[2])
  }
  return(g)
}

gff_data <- load_gff(opt$gff_file, opt$window_range)
if (!is.null(opt$gff_tracks) && !is.null(gff_data)) {
  gff_data <- gff_data %>% filter(type %in% unlist(strsplit(opt$gff_tracks, ",")))
}
ngaps_data <- load_gff(opt$ngaps_file, opt$window_range)

trf_data <- NULL
if (!is.null(opt$trf_file)) {
  trf_data <- readr::read_tsv(opt$trf_file, show_col_types = FALSE)
  if (!is.null(opt$window_range)) {
    crd <- as.numeric(unlist(strsplit(opt$window_range, ",")))
    trf_data <- trf_data %>% filter(Start >= crd[1], End <= crd[2])
  }
}

# Pivot and color setup
values <- values %>% tidyr::pivot_longer(cols = starts_with(c("A", "C", "G", "T")))
palette_csv <- file.path(script_dir, "includes", "paletteMatrix_base.csv")

if (file.exists(palette_csv)) {
  palette_order <- read.csv(palette_csv, header = TRUE)
  triplet_names <- unique(values$name)
  palette_names <- c()
  for (col in 1:ncol(palette_order)) {
    for (row in 1:nrow(palette_order)) {
      if (palette_order[row, col] %in% triplet_names) {
        palette_names <- c(palette_names, as.character(palette_order[row, col]))
      }
    }
  }
} else {
  palette_names <- unique(values$name)
}
triplet_colors <- sapply(palette_names, palette_builder, palette = opt$palette_type)

# --- Plotting Execution ---

if (opt$is_circular) {
  out_name <- paste0(opt$output_prefix, ".", opt$output_format)
  p <- create_circular_plot(values, triplet_colors, legend = opt$show_legend,
                          frequencies = opt$is_frequencies, ylims = opt$use_ylims,
                          limit = opt$graph_length)
  if (opt$is_transparent) {
    p <- p + theme(panel.background = element_rect(fill = "transparent"),
                 plot.background = element_rect(fill = "transparent", color = NA))
  }
  ggsave(filename = out_name, device = opt$output_format, width = 10, height = 10,
         units = "in", dpi = opt$resolution, limitsize = FALSE, bg = "transparent")
} else {
  for (seq in seq_names) {
    out_name <- paste0(opt$output_prefix, "_", seq, ".", opt$output_format)
    temp_values <- values %>% filter(Sequence == seq)

    # Scaling logic
    if (opt$keep_scale) {
      temp_range <- (max(temp_values$End) - min(temp_values$Start) + 1) / (1e6 * opt$x_scale)
      temp_length <- temp_range + (if (opt$show_legend) 2 else 0.5)
      if (!opt$show_axes) temp_length <- temp_length - 0.5
    } else {
      temp_length <- 10
    }

    height_factor <- if (length(lib_names) > 1) length(unique(temp_values$Library)) else 1

    p <- create_spectra_plot(temp_values, triplet_colors, legend = opt$show_legend,
                           facet = (length(lib_names) > 1), frequencies = opt$is_frequencies,
                           ylims = opt$use_ylims, scale = opt$x_scale, axes = opt$show_axes)

    if (!is.null(ngaps_data)) {
      seq_gaps <- ngaps_data %>% filter(seqid == seq)
      if (nrow(seq_gaps) > 0) {
        p <- p + geom_rect(data = seq_gaps, aes(xmin = start - 10000, xmax = end + 10000,
                                             ymin = 0, ymax = 1),
                          inherit.aes = FALSE, fill = "black")
      }
    }

    if (!is.null(trf_data)) {
      seq_trf <- trf_data %>% filter(Sequence == seq)
      if (nrow(seq_trf) > 0) {
        p <- p + geom_line(data = seq_trf, aes(x = (End + Start) / 2, y = (Proportion - 1.12) / 4),
                          color = "black", size = 0.25) +
          scale_y_continuous(limits = c(-0.28, 1), expand = c(0, 0))
      }
    }

    if (opt$is_transparent) {
      p <- p + theme(panel.background = element_rect(fill = "transparent"),
                   plot.background = element_rect(fill = "transparent", color = NA))
    }

    ggsave(filename = out_name, device = opt$output_format, width = temp_length,
           height = 1 + height_factor * 2, units = "in", dpi = opt$resolution, limitsize = FALSE)

    # Optional GFF track plot
    if (!is.null(gff_data)) {
      seq_gff <- gff_data %>% filter(seqid == seq)
      if (nrow(seq_gff) > 0) {
        height_gff <- length(unique(seq_gff$type))
        pg <- ggplot(seq_gff, aes(xmin = start, xmax = end, ymin = 0.1, ymax = 0.2, fill = strand)) +
          geom_rect() +
          scale_y_continuous(limits = c(0.1, 0.2), expand = c(0, 0)) +
          scale_x_continuous(expand = c(0, 0)) +
          theme_bw() + ylab("") + xlab("Position (bp)") +
          theme(legend.position = "none", axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                strip.text.y = element_text(size = 4, colour = "black", angle = 90)) +
          facet_grid(rows = vars(type))

        if (opt$is_transparent) {
          pg <- pg + theme(panel.background = element_rect(fill = "transparent"),
                         plot.background = element_rect(fill = "transparent", color = NA))
        }

        gff_out <- paste0(opt$output_prefix, "_gff_", seq, ".", opt$output_format)
        ggsave(filename = gff_out, device = opt$output_format, width = temp_length + 0.26,
               height = 0.6 + (0.35 * height_gff), units = "in", dpi = opt$resolution, limitsize = FALSE)
      }
    }
  }
}
