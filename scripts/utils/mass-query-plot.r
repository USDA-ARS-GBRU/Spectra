#!/usr/bin/env Rscript
####
##### Plot highest and lowest kmers along a sequence based on mass-query.py output

#Load dependencies
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(optparse)
})

baseTheme = theme_bw() + theme(
  axis.title.x = element_text(size=5),
  axis.title.y = element_text(size=5),
  legend.key.size = unit(5, "pt"),
  legend.spacing = unit(-3, "pt"),
  axis.text.x = element_text(angle=0, vjust=0.5, size=5),
  axis.text.y = element_text(size=5),
  line = element_blank(),
  axis.ticks = element_line(),
  legend.title = element_blank()
)

coveragePlot = function(values, legend=FALSE, range, scale, axes){
  values = values %>% group_by(Sequence, Start, End) %>%
    arrange(desc(Bin), .by_group=TRUE, ) %>%
    mutate(ymax = cumsum(Count), ymin = ymax - Count) %>%
    ungroup()
  p=ggplot(values, aes(x=(Start+End)/2, ymin=ymin, ymax=ymax, fill=Bin))+geom_ribbon()+baseTheme
  
  xrange <- max(values$End) - min(values$Start) + 1
  if(log10(xrange) > log10(scale*1000000)+1){
    scale = scale * 10
  }
  breaks <- xrange %/% (scale*1000000)
  if (breaks < 2) breaks = 2
  p = p + scale_x_continuous(
      limits=c(min(values$Start), max(values$End)),
      n.breaks=breaks,
      expand=c(0,0),
      labels=scales::scientific
    ) +
    xlab("Window Position (nucleotide)") +
    ylab("Counts") +
    baseTheme + theme(plot.margin = margin(t=2.5, l=2.5, b=2.5, r=2.5))
  if (!axes) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(t=0, l=-2.7, b=-2.7, r=0)
    )
  }
  return(p)
}

option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL, help="coverage tsv file' [default %default]", dest="input_filename"),
  make_option(c("-w", "--window-size"), type="character", default=NULL, help="Generate plot only from the closest intervals of \"N,M\". Otherwise generate a plot for all positions [default %default]", dest="window_size"),
  make_option(c("-s", "--sequence"), type="character", default=NULL, help="Generate plot only from the sequences with names in \"A,B,C\" or by regular expression if -e flag specified. Otherwise generate a plot for each sequence id [default %default]", dest="sequences"),
  make_option(c("-e", "--regex"), action="store_true", default=FALSE, help="Uses regex to subset sequence and library names [default %default]", dest="regex"),
  make_option(c("-o","--output-prefix"), type="character", default=NULL, help="Output prefix [default %default]", dest="output_filename"),
  make_option(c("-f","--output-format"), type="character", default="png", help="Output image format [default %default]", dest="output_type"),
  make_option(c("-r","--resolution"), type="numeric", default=300, help="Plotting dpi resolution [default %default]", dest="resolution"),
  make_option(c("-l","--show-legend"), action="store_true", default=FALSE, help="Display triplet color legend [default %default]", dest="legend"),
  make_option(c("-x","--scale"), type="numeric", default=1, help="Scale of x-axis. Plot each n (mb) over 1 inch [default %default]", dest="scale"),
  make_option(c("-a","--axes"), action="store_false", default=TRUE, help="Display axes text [default %default]", dest="axes"),
  make_option(c("-k","--keep-scale"), action="store_true", default=FALSE, help="Incorporate scale [default %default]", dest="keep")
)

options(error=traceback)
parser = OptionParser(usage = "%prog -i coverage.tsv [options]",option_list=option_list)
opt = parse_args(parser)


# Prepare outputfile space. If no output given
if(is.null(opt$output_filename)){
  split_file = strsplit(opt$input_filename,"/")[[1]]
  output_file = c(split_file[length(split_file)],opt$output_type)
}else{
  if(opt$output_type=="svg"){
    suppressPackageStartupMessages(library(svglite))
  }
  output_file = c(opt$output_filename, opt$output_type)
}

# read tsv values or exit with error
if(is.null(opt$input_filename)){
  cat("Error: No tsv specified. See usage 'with spectra-plot.r -h'\n")
  quit()
}else{
  values = readr::read_tsv(opt$input_filename, show_col_types = FALSE)
}

if(!is.null(opt$sequences)){
  if(opt$regex){
    values = values %>% filter(grepl(opt$sequences, Sequence))
  }else{
    values = values %>% filter(Sequence%in%as.vector(opt$sequences))
  }
}

if(!is.null(opt$window_size)){
  coords = strsplit(opt$window_size,",")
  values = values %>% filter(Start >= as.numeric(coords[[1]][1]))
  values = values %>% filter(End <= as.numeric(coords[[1]][2]))
}
height.factor = 1
seq.names = unique(values$Sequence)
for(seq in seq.names){
  temp.values = values %>% filter(Sequence==seq)
  
  # Calculate the necessary size to match the current scale
  if(opt$keep){
    temp.range =  (max(temp.values$End) - min(temp.values$Start) + 1) / (1000000 * opt$scale)
    if(opt$legend){
      temp.length = temp.range + 2
    }else{
      temp.length = temp.range + 0.5
    }
    if(!opt$axes){
      temp.length = temp.length - 0.5
    }
  }else{
    temp.length=10
  }
  
  p_low = coveragePlot(temp.values%>%filter(as.integer(sub("pct", "", Bin)) <= 50), legend=opt$legend, temp.range, opt$scale, opt$axes)
  ggsave(filename=paste0(output_file[1], '_', seq, '_low.', output_file[2]),device=output_file[2], width=temp.length, height=1+height.factor*2, units="in", dpi=opt$resolution, limitsize=F)
  p_high = coveragePlot(temp.values%>%filter(as.integer(sub("pct", "", Bin)) > 50), legend=opt$legend, temp.range, opt$scale, opt$axes)
  ggsave(filename=paste0(output_file[1], '_', seq, '_high.', output_file[2]),device=output_file[2], width=temp.length, height=1+height.factor*2, units="in", dpi=opt$resolution, limitsize=F)
}