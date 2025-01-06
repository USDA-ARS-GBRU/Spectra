#!/usr/bin/env Rscript
####
##### Repeats along (catermerize) sequence for multiple libraries

#load dependencies
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))

spectraPlot = function(values, tripletColors, legend=FALSE, facet=FALSE, frequencies=FALSE, ylims=TRUE, range, scale, axes, paletteNames){

	if(frequencies){
		p = ggplot() + geom_area(data=values, aes(fill=name, x=(End+Start)/2, y=value), stat="identity", position="stack")
	}else{
		p = ggplot() + geom_area(data=values, aes(fill=name, x=(End+Start)/2,y=value/(End-Start+1)), stat="identity", position="stack")
	}
	if(ylims){
        p = p + scale_y_continuous(limits=c(0,1), expand=c(0, 0))
	}else{
		p = p + scale_y_continuous(expand=c(0, 0))
	}
	####
    range = max(values$End)-min(values$Start)+1
	####

	breaks = range %/% (scale*1000000)
	# automatically lower the scale if sequence is too small (less than 1 scalar)
	if(breaks < 2){
	    breaks=2
	}
	# reminder: change n.breaks=10 to n.breaks=breaks
	#
	p = p + scale_fill_manual(values=tripletColors) +
	    scale_x_continuous(
            limits=c(min(values$Start), max(values$End)),
            breaks=waiver(),
            minor_breaks=waiver(),
            n.breaks=10,
            expand=c(0, 0),
            labels=scales::scientific
		) +
		xlab("Window Position (nucleotide)") +
		ylab("Proportion") +
		theme_bw() +
		theme(
		    axis.title.x = element_text(size=5),
			axis.title.y = element_text(size=5),
			legend.text = element_text(size=5),
			legend.key.size = unit(5, "pt"),
			legend.spacing = unit(-3, "pt"),
			axis.text.x = element_text(angle=0, vjust=0.5, size=5),
			axis.text.y = element_text(size=5),
			line = element_blank(),
			axis.ticks = element_line(),
			legend.title = element_blank(),
			plot.margin = margin(t=2.5, l=2.5, b=2.5, r=2.5)
		)

	if(facet){
		p = p + facet_grid(rows = vars(values$Library))
	}
	if(!legend){
		p = p + theme(legend.position = "none")
	}
	if(!axes){
		p = p +
		theme(
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

circularPlot = function(values, tripletColors, legend=FALSE, frequencies=FALSE, ylims=TRUE, axes=FALSE, paletteNames, limit=0){
	suppressPackageStartupMessages(library(egg))

	#definable sequence spacer size
	spacerLength = 1000000
	spacerValues = rbind(
		values %>% filter(End==min(values$End),Sequence==values$Sequence[1]) %>% mutate(Start=1, End=1000, value=0, Sequence='Spacer'),
		values %>% filter(End==min(values$End),Sequence==values$Sequence[1]) %>% mutate(Start=1001, End=spacerLength-1000, value=0, Sequence='Spacer'),
		values %>% filter(End==min(values$End),Sequence==values$Sequence[1]) %>% mutate(Start=spacerLength-999, End=spacerLength, value=0, Sequence='Spacer')
	)
	
	# modify sequence positions to be absolute
	sequences_names=names(table(values$Sequence))
	newValues = values %>% filter(Sequence==sequences_names[1])
	current_max = max(newValues$End)
	
	for(index in sequences_names[2:length(sequences_names)]){
		newValues = rbind(
			newValues,
			spacerValues %>% mutate(Start=Start+current_max,End=End+current_max)
		)
		current_max = current_max + spacerLength
		current_values = values %>% filter(Sequence==index)
		newValues = rbind(
			newValues,
			current_values %>% mutate(Start=Start+current_max, End=End+current_max)
		)
		current_max = max(newValues$End)
	}
	
	if(frequencies){
		p = ggplot() + geom_area(data=newValues, aes(fill=name, x=(Start+End)/2, y=value), stat="identity", position="stack")
	}else{
		p = ggplot() + geom_area(data=newValues, aes(fill=name, x=(Start+End)/2,y=value/(End-Start+1)), stat="identity", position="stack")
	}
	if(ylims){
        p = p + scale_y_continuous(limits=c(0,1), expand=c(.5, .5))
	}else{
		p = p + scale_y_continuous(expand=c(0, 1))
	}
	p = p + scale_fill_manual(values=tripletColors)
	if(limit>0){
		plot_size = limit + ((length(sequences_names)-1) * spacerLength)
	}else{
		plot_size = current_max
	}
	p = p + scale_x_continuous(
		limits=c(1, plot_size),
		breaks=waiver(),
		minor_breaks=waiver(),
		n.breaks=30,
		expand=c(0, 0),
		labels=scales::scientific
	)+
	xlab("Window Position (nucleotide)") +
	ylab("Proportion") +
	theme_bw() +
	theme(
		axis.title.x = element_text(size=5),
		axis.title.y = element_text(size=5),
		legend.text = element_text(size=5),
		legend.key.size = unit(5, "pt"),
		legend.spacing = unit(-3, "pt"),
		axis.text.x = element_text(angle=0, vjust=0.5, size=5),
		axis.text.y = element_text(size=5),
		line = element_blank(),
		axis.ticks = element_line(),
		legend.title = element_blank(),
		plot.margin = margin(t=0, l=0, b=0, r=0)
	)# +
	#facet_grid(cols = vars(values$Sequence))
	
	if(!legend){
		p = p + theme(legend.position = "none")
	}
	if(!axes){
		p = p +
		theme(
			axis.text.x = element_blank(),
			axis.text.y = element_blank(),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.ticks = element_blank(),
			plot.margin = margin(t=0, l=-2.7, b=-2.7, r=0)
		)
	}
	p = p + coord_polar(start = 0)
	return(p)
}

trfPlot = function(){}

paletteBuilder = function(triplet,palette=opt$palette){
	bases = c("A","C","G","T")
	colors = switch(
		palette,
		'base'= c("C6","6C","3C","10"),
		'base2'= c("BF","6F","46","1F"),
		'base3'= c("CF","6F","56","2F"),
		'dual'= c("96","3C","3C","96"),
		'acontrast'= c("CC","33","33","3F")
	)

	color = paste0("#",colors[which(substr(triplet,1,1) == bases)],colors[which(substr(triplet,2,2) == bases)],colors[which(substr(triplet,3,3) == bases)])
	return(color)
}

#Argument parser
option_list <- list(
	make_option(c("-i","--input"), type="character", default=NULL, help="triplet tsv file' [default %default]", dest="input_filename"),
	make_option(c("-w", "--window-size"), type="character", default=NULL, help="Generate plot only from the closest intervals of \"N,M\". Otherwise generate a plot for all positions [default %default]", dest="window_size"),
	make_option(c("-s", "--sequence"), type="character", default=NULL, help="Generate plot only from the sequences with names in \"A,B,C\" or by regular expression if -e flag specified. Otherwise generate a plot for each sequence id [default %default]", dest="sequences"),
	make_option(c("-n", "--libraries"), type="character", default=NULL, help="Generate plot only from the libraries with names in \"A,B,C\" or by regular expression if -e flag specified. Otherwise generate a plot for each library id [default %default]", dest="libraries"),
	make_option(c("-e", "--regex"), action="store_true", default=FALSE, help="Uses regex to subset sequence and library names [default %default]", dest="regex"),
	make_option(c("-u", "--graphlength"), type="numeric", default=0, help="Designate the length of sequence to plot (for partial graphs in circular plot). 0 = No limit [default %default]", dest="length"),
	make_option(c("-g", "--gff-file"), type="character", default=NULL, help="Generate plot of overlapping gene annotations from supplied gff. [default %default]", dest="gffFile"),
	make_option(c("-t", "--gff-tracks"), type="character", default=NULL, help="Curate which gff types to use in types \"A,B,C\". [default %default]", dest="gffTracks"),
	make_option(c("-z", "--trf-file"), type="character", default=NULL, help="Generate plot of overlapping trf annotations from supplied trf-tsv. [default %default]", dest="trfFile"),
	make_option(c("-l","--show-legend"), action="store_true", default=FALSE, help="Display triplet color legend [default %default]", dest="legend"),
	make_option(c("-o","--output-file"), type="character", default=NULL, help="Output filename [default %default]", dest="output_filename"),
	make_option(c("-r","--resolution"), type="numeric", default=300, help="Plotting dpi resolution [default %default]", dest="resolution"),
	make_option(c("-f","--freq"), action="store_true", default=FALSE, help="Data already transformed to frequencies [default %default]", dest="frequencies"),
	make_option(c("-y","--ylims"), action="store_false", default=TRUE, help="Limit results to y-axes between 0,1 [default %default]", dest="ylims"),
	make_option(c("-x","--scale"), type="numeric", default=1, help="Scale of x-axis. Plot each n (mb) over 1 inch [default %default]", dest="scale"),
	make_option(c("-a","--axes"), action="store_false", default=TRUE, help="Display axes text [default %default]", dest="axes"),
	make_option(c("-k","--keep-scale"), action="store_true", default=FALSE, help="Incorporate scale [default %default]", dest="keep"),
	make_option(c("-p","--palette"), type="character", default='base', help="Spectra color palette. Available palettes: base, dual [default %default]", dest="palette"),
	make_option(c("-c","--circular"), action="store_true", default=FALSE, help="Invoke a circular, full genome plot [default %default]", dest="circular")
)
options(error=traceback)
parser = OptionParser(usage = "%prog -i triplet.tsv [options]",option_list=option_list)
opt = parse_args(parser)

# prepare outputfile space
if(is.null(opt$output_filename)){
	split_file = strsplit(opt$input_filename,"/")[[1]]
	output_file = c(split_file[length(split_file)],"png")
}else{
	split_file = strsplit(opt$output_filename,"\\.")[[1]]
	output_type = split_file[length(split_file)]
	if(output_type=="svg"){
		suppressPackageStartupMessages(library(svglite))
	}
	output_file = c(split_file[1], output_type)
}

# script info for finding utils folder
scriptCommands = commandArgs(trailingOnly = FALSE)
scriptArg = "--file="
scriptLocation = sub(scriptArg, "", scriptCommands[grep(scriptArg, scriptCommands)])

# gff preparation
gff = NULL
if(!is.null(opt$gffFile)){
	suppressPackageStartupMessages(library(ape))
	gff = read.gff(opt$gffFile)
	# filter
	if(!is.null(opt$gffTracks)){
		gff = gff %>% filter(type%in%as.vector(opt$gffTracks))
	}
	if(!is.null(opt$window_size)){
		coords = strsplit(opt$window_size,",")
		gff = gff %>% filter(start >= as.numeric(coords[[1]][1]))
		gff = gff %>% filter(end <= as.numeric(coords[[1]][2]))
	}
}

# trf preparation
trf = NULL
if(!is.null(opt$trfFile)){
	trf = read.csv(opt$trfFile,sep="\t")
		# filter
	if(!is.null(opt$window_size)){
		coords = strsplit(opt$window_size,",")
		trf = trf %>% filter(start >= as.numeric(coords[[1]][1]))
		trf = trf %>% filter(end <= as.numeric(coords[[1]][2]))
	}
}

# read tsv values or exit with error
if(is.null(opt$input_filename)){
  cat("Error: No tsv specified. See usage 'with spectra-plot.r -h'\n")
  quit()
}else{
	values = read.csv(opt$input_filename,sep="\t",stringsAsFactors = FALSE)
}

# filter out wasteful data
if(!is.null(opt$libraries)){
	if(opt$regex){
		values = values %>% filter(grepl(opt$libraries, Library))
	}else{
		values = values %>% filter(Library%in%as.vector(opt$libraries))
	}
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

lib.names = names(table(values$Library))
seq.names = names(table(values$Sequence))

if(length(lib.names)>1 && !is.null(opt$gffFile)){
	print("Warning: the gff track will be omitted because multiple libraries are being plotted.")
	gff = NULL
}

# Pivot table for stacking of columns
values = values %>% tidyr::pivot_longer(cols=starts_with(c("A","C","G","T")))

# Color and Palette building
paletteOrderDF = read.csv(paste0(dirname(scriptLocation),"/includes/paletteMatrix_base.csv"))

tripletNames=names(table(values$name))
paletteNames = c()
for(colNum in 1:ncol(paletteOrderDF)){
	for(rowNum in 1:nrow(paletteOrderDF)){
		if(paletteOrderDF[rowNum,colNum] %in% tripletNames){
			paletteNames = c(paletteNames,paletteOrderDF[rowNum,colNum][1])
		}
	}
}

tripletColors=sapply(paletteNames,paletteBuilder)

if(opt$circular){
	# write multiple plots in a single frame if sequence names are the same
		seq.filename = paste0(output_file[1], '.', output_file[2])

		p = circularPlot(values, tripletColors, legend=opt$legend, frequencies=opt$frequencies, ylims=opt$ylims, opt$axes, paletteNames, opt$length)
		trackOffset=-0.03
		ggsave(filename=seq.filename,device=output_file[2], width=10, height=10, units="in", dpi=opt$resolution, limitsize=F)
}else{
	# write multiple plots in a single frame if sequence names are the same
	for(seq in seq.names){
		seq.filename = paste0(output_file[1], '_', seq, '.', output_file[2])
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

		if(length(lib.names)>1){
			faceted=TRUE
			height.factor = length(names(table(temp.values$Library)))
		}else{
			faceted=FALSE
			height.factor = 1
		}
		p = spectraPlot(temp.values, tripletColors, legend=opt$legend, facet=faceted, frequencies=opt$frequencies, ylims=opt$ylims, temp.range, opt$scale, opt$axes, paletteNames)
		trackOffset=-0.03
		if(!is.null(gff)){
			p = p + geom_segment(data=gff%>%filter(seqid==seq), aes(x=start, xend=end, y=trackOffset, yend=trackOffset, color=strand), size=4) + scale_y_continuous(limits=c(-.06, 1), expand=c(0, 0))
			trackOffset = trackOffset - .03
		}
		if(!is.null(trf)){
			if(nrow(trf%>%filter(Sequence==seq))>0){
				p = p + geom_line(data=trf%>%filter(Sequence==seq), aes(x=(End+Start)/2, y=(Proportion-1.12)/4), color="black", size=0.25) + scale_y_continuous(limits=c(-.28, 1), expand=c(0, 0))
			}
			trackOffset = trackOffset - .03
		}
		ggsave(filename=seq.filename,device=output_file[2], width=temp.length, height=1+height.factor*2, units="in", dpi=opt$resolution, limitsize=F)
	}
}
