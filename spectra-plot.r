#!/usr/bin/env Rscript
####
##### Repeats along (catermerize) sequence for multiple libraries


#load dependencies
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))

# Functions for processing
triplerColor = function(triplet){
  bases = c("A","C","G","T")
  colors = c("DD","9F","60","11")
  #colors = c("CC","99","66","33")
  
  color = paste0("#",colors[which(substr(triplet,1,1) == bases)],colors[which(substr(triplet,2,2) == bases)],colors[which(substr(triplet,3,3) == bases)])
  return(color)
}

#Argument parser
option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL, help="triplet tsv file' [default %default]", dest="input_filename"),
  make_option(c("-w", "--window-size"), type="character", default=NULL, help="if specified, generate plot only from the closest intervals of \"N,M\". Otherwise generate a plot for all positions [default %default]", dest="window_size"),
  make_option(c("-s", "--sequence"), type="character", default=NULL, help="if specified, generate plot only from the sequences with names in \"A,B,C\". Otherwise generate a plot for  each sequence id [default %default]", dest="sequences"),
  make_option(c("-g", "--gff"), type="character", default=NULL, help="if specified, generate plot of overlapping gene annotations. [default %default]", dest="gff"),
  make_option(c("-l","--show-legend"), action="store_true", default=FALSE, help="display triplet color legend [default %default]", dest="legend"),
  make_option(c("-o","--output-file"), type="character", default=NULL, help="output filename [default %default]", dest="output_filename"),
  make_option(c("-c","--collapse"), action="store_true", default=FALSE, help="collapse reverse compliments [default %default]", dest="reverse"),
  make_option(c("-r","--resolution"), type="numeric", default=300, help="plotting dpi resolution [default %default]", dest="resolution"),
  make_option(c("-S","--simplify"), action="store_true", default=FALSE, help="only display the first of each rotated triplet [default %default]", dest="simplify_bases"),
  make_option(c("-f","--freq"), action="store_true", default=FALSE, help="data already transformed to frequencies [default %default]", dest="frequencies")
  
)
options(error=traceback)
parser <- OptionParser(usage = "%prog -i triplet.ysv [options]",option_list=option_list)
opt = parse_args(parser)

# Read triplet tsv
if(is.null(opt$input_filename)){
  cat("No tsv specified. See usage \"with spectra-plot.r -h\"")
  quit()
}else{
  values = read.csv(opt$input_filename,sep="\t",stringsAsFactors = FALSE)
}

if(!is.null(opt$sequences)){
  values = values %>% filter(Sequence%in%as.vector(opt$sequences))
}

if(!is.null(opt$window_size)){
  coords = strsplit(opt$window_size,",")
  values = values %>% filter(Start >= as.numeric(coords[[1]][1]))
  values = values %>% filter(End <= as.numeric(coords[[1]][2]))
}

values = values %>% tidyr::pivot_longer(cols=starts_with(c("A","C","G","T")))

if(opt$simplify_bases){
  values$name = unname(sapply(values$name,simplify))
}

tripletNames=names(table(values$name))
tripletColors=sapply(tripletNames,triplerColor)

size = c(min(values$Start),max(values$End))
height_factor = length(table(values$Sequence))

if(opt$frequencies){
  p = ggplot() +
    geom_area(data=values, aes(fill=name,x=(End+Start)/2,y=(value)), stat="identity",position="stack")
}else{
  p = ggplot() +
    geom_area(data=values, aes(fill=name,x=(End+Start)/2,y=(value)/(End-Start+1)), stat="identity",position="stack")
}
p = p + scale_fill_manual(values=tripletColors) +
  scale_y_continuous(limits=c(0,1),expand = c(0, 0)) +
  scale_x_continuous(limits=size,breaks=waiver(),minor_breaks = waiver(),n.breaks = 10, expand = c(0, 0)) +
  xlab("Window Position (nucleotide)") +
  ylab("Proportion") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size=6),
    axis.title.y = element_text(size=6),
    legend.text = element_text(size=5),
    legend.key.size = unit(5,"pt"),
    legend.spacing = unit(-3,"pt"),
    axis.text.x = element_text(angle=0,vjust=0.5),
    line=element_blank(),
    axis.ticks = element_line(),
    legend.title= element_blank()
  )

if(height_factor>1){
  p = p + facet_grid(rows = vars(values$Sequence))
}

if(!opt$legend){
  p = p + theme(legend.position = "none")
}

if(!is.null(opt$gff) && !is.null(opt$sequences)){
  gff = read.csv(opt$gff, sep='\t',header=FALSE)
  genes = gff[gff$V1 %in% as.vector(opt$sequences),]
  if(is.null(opt$window_size)){
    window = c(min(values$Start),c(max(values$End)))
    
  }
  names(genes) = c("Sequence","Source","Type","Start","End","t1","Direction","t2","Info")
  gff.plot = geom_segment(data=genes, aes(x=Start,xend=End,y=-.025,yend=-.025,color=Direction),size = 3)
  p = p + gff.plot + scale_y_continuous(limits=c(-.05,1), expand = c(0, 0)) 
}

if(is.null(opt$output_filename)){
  split_file = strsplit(opt$input_filename,"/")[[1]]
  output_file = paste0(split_file[length(split_file)],".png")
}else{
  split_file = strsplit(opt$input_filename,".")[[1]]
  output_type = split_file[length(split_file)]
  if(output_type=="svg"){
    suppressPackageStartupMessages(library(svglite))
  }
  output_file = opt$output_filename
}

#ggsave(filename=output_file,device=opt$output_filetype,width=10,height=1+height_factor*2,units="in",dpi=300,limitsize=F)
ggsave(filename=output_file,width=10,height=1+height_factor*2,units="in",dpi=opt$resolution,limitsize=F)
