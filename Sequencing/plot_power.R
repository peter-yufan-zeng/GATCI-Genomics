### plot_power.R ##################################################################################
# Generate customized power plot (modified from that output by SeqSeq/SMPower).

### LIBRARIES AND SETUP ###########################################################################
library(BoutrosLab.plotting.general);

# set working directory
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/Exomes/config_files/');

# read in config file
cfg <- read.config.file();

### READ DATA  ####################################################################################
# get plot arguments
max_sample_size <- 250;
vertical_line_locations <- 139;

# get power results
power <- read.delim('/u/sprokopec/ATC/Exomes/SNV_callers/somaticsniper_1.0.5.0/Filtered/PowerAnalysis/2017-08-08_Output/power/power_cds.txt');
plot.data <- power[power$sample.size <= max_sample_size,];

p1 <- unique(plot.data$p1);
myColours <- RColorBrewer::brewer.pal(max(3, length(p1)), name = 'Set2')[1 : length(p1)];

### CREATE PLOT ###################################################################################
# go to output directory
setwd(cfg$plot.directory);

# make some legends
legends <- list(
	legend = list(
		title = 'Type I Error',
		labels = '0.05',
		size = 0,
		border = FALSE
		),
	legend = list(
		title = 'BMR',
		labels = scientific.notation(unique(plot.data$p0), digits = 2),
		size = 0,
		border = FALSE
		),
	legend = list(
		title = 'p1',
		labels = paste0(c(0.5,1,2,5,10),'%'),
		colours = myColours
		)
	);

myLegend.grob <- legend.grob(
	legends,
	title.just = 'left'
	);

# make a power plot!
create.scatterplot(
	power ~ sample.size,
	plot.data,
	filename = generate.filename('power', 'cds', 'tiff'),
	groups = plot.data$p1,
	xlab.label = 'Sample Size',
	ylab.label = 'Power',
	abline.v = vertical_line_locations,
	abline.lty = 2,
	abline.lwd = 3,
	abline.col = 'grey50',
	lwd = 3,
	type = 'l',
	col = myColours,
	xat = seq(0, max_sample_size, max_sample_size / 10),
	xlimits = c(0, max_sample_size),
	yat = seq(0, 1, 0.25),
	ylimits = c(0, 1.02),
	xaxis.cex = 0.75,
	yaxis.cex = 0.75,
	xlab.cex = 1,
	ylab.cex = 1,
	xaxis.tck = 0.5,
	yaxis.tck = 0.5,
	width = 11,
	right.padding = 5,
	legend = list(
		right = list(fun = myLegend.grob)
		),
	resolution = 800
	);

### SessionInfo ###################################################################################
save.session.profile(generate.filename('Power','Plot','txt'));
