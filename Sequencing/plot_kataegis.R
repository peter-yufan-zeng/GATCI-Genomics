### plot_kataegis.R ###############################################################################
# Purpose: WGS was performed on ATC tumour and normal samples to identify genes containing 
#	recurrent somatic SNVs in tumour samples.
#	SNVs were called using SomaticSniper and resulting vcfs were evaluated for kataegis events
#	using seqKat (v1.1).
# Initial write: 2016-05-02
# Update: 2017-08-15

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general);

# general parameters
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/WholeGenome/kataegis/');

# general parameters
date <- Sys.Date();

### READ IN DATA FROM RECSNV ######################################################################
# get cytoband locations
cytobands <- read.delim('/u/sprokopec/ATC/other/cytoBandpositions_hg38.txt');
centromeres <- cytobands[cytobands$gieStain == 'acen',];

# get collapsed summary
summaryData <- read.delim('2019-04-11_combined_results.collapsed.summary');

# get significant event data
eventData <- read.delim('2019-04-11_combined_results.full.summary');
tmpData <- eventData[which(eventData$binary.kat == 1 & eventData$significance == 2),];
tmpData <- tmpData[which(tmpData$p.value < 0.05),];
tmpData <- tmpData[which(tmpData$num.snvs >= 10),];
tmpData <- unique(tmpData[which(tmpData$num.tcx > 4),c('sample.name','chromosome')]);
names(tmpData) <- c('sample','chr');

### VISUALIZATION ##################################################################################
# indicate colour scheme
colour.scheme <- list(
	colours = rev(c("goldenrod1", "hotpink1", "olivedrab3", "dodgerblue2", "black", "firebrick2")),
	levels = c("T>A", "T>C", "T>G", "C>A", "C>T", "C>G")
	);

# make the key
base.change.key <- list(
	points = list(
		col = rev(colour.scheme$colours),
		pch = 19,
		cex = 1,
		alpha = 0.75
		),
	text = list(
		lab = c(
			expression("T"%->%"A / A"%->%"T"),
			expression("T"%->%"C / A"%->%"G"),
			expression("T"%->%"G / A"%->%"C"),
			expression("C"%->%"A / G"%->%"T"),
			expression("C"%->%"T / G"%->%"A"),
			expression("C"%->%"G / G"%->%"C")
			),
		cex = 0.7
		),
	space = 'top',
	columns = 6,
	between = 0.7 
	);

# organize data
allData <- list();
for (i in 1:nrow(tmpData)) {

	sample.id <- as.character(tmpData$sample[i]);
	chr <- as.character(tmpData$chr[i]);

	# get significant event data
	these.files <- sort(list.files(path = paste0(sample.id,'/',sample.id), pattern = paste0('finaltable_chr',chr, '_')));
	signif.events <- read.delim(paste0(sample.id,'/',sample.id,'/',these.files[1]));
	signif.events <- signif.events[which(signif.events$significance == 2),];
	signif.events <- signif.events[which(signif.events$num.tcx/signif.events$num.snvs > 0.1),];
	
	signif.events$mid <- apply(
		signif.events[,c('window.start','window.end')],
		1,
		function(i) { i[1] + (i[2] - i[1])/2 }
		);

	allData[[sample.id]][[chr]]$signif <- signif.events;

	# get SNV data
	tmp <- read.delim(paste0(sample.id, '/', sample.id, '_snvs.bed'));
	my.chr <- if (chr == 23) { 'X' } else if (chr == 24) { 'Y' } else { chr }
	tmp <- tmp[which(tmp$CHR == paste0('chr', my.chr)),];

	# determine intermutational distance
	tmp$inter.distance <- NA;

	for (snv in 2:nrow(tmp)) {
		tmp[snv,]$inter.distance <- tmp[snv,]$POS - tmp[snv-1,]$POS;
		}

	# indicate and refactor base change
	tmp$base.change <- paste0(tmp$REF, '>', tmp$ALT);
	tmp$base.change[tmp$base.change == 'A>T'] <- 'T>A';
	tmp$base.change[tmp$base.change == 'A>C'] <- 'T>G';
	tmp$base.change[tmp$base.change == 'A>G'] <- 'T>C';
	tmp$base.change[tmp$base.change == 'G>A'] <- 'C>T';
	tmp$base.change[tmp$base.change == 'G>T'] <- 'C>A';
	tmp$base.change[tmp$base.change == 'G>C'] <- 'C>G';
	tmp$base.change <- factor(
		tmp$base.change,
		levels = c("T>A","T>C","T>G","C>A","C>T","C>G")
		);

	# indicate base change colours
	tmp$colour <- colour.scheme$colours[factor(tmp$base.change, levels = colour.scheme$levels)];

	# format values for plotting
	tmp$x <- tmp$POS/10**6;
	tmp$y <- log10(tmp$inter.distance);

	allData[[sample.id]][[chr]]$data <- tmp;

	}

# make the plots!
setwd('plots');

for (i in 1:nrow(tmpData)) {

	sample.id <- as.character(tmpData$sample[i]);
	chr <- as.character(tmpData$chr[i]);
	my.chr <- if (chr == 23) { 'X' } else if (chr == 24) { 'Y' } else { chr }

	signifEvents <- allData[[sample.id]][[chr]]$signif;
	signifEvents <- signifEvents[which(signifEvents$num.snvs >= 10),];

	myData <- allData[[sample.id]][[chr]]$data;

	textData <- data.frame(
		labels = paste0('p = ', round(signifEvents$p.value, 3)),
		start = signifEvents$window.start/10**6,
		end = signifEvents$window.end/10**6,
		keep = 1,
		x.pos = signifEvents$mid/10**6,
		y.pos = NA
		);

	if (nrow(textData) > 1) {
		for (j in 2:nrow(textData)) {
			if (textData[j-1,]$end > textData[j,]$start) { textData[j,]$keep <- 0; }
			}
		}

	textData <- textData[which(textData$keep == 1),];

	# make the plot!
	create.scatterplot(
		y ~ x,
		myData,
		filename = generate.filename(sample.id, paste0('Chromosome', my.chr, '_kataegis'), 'tiff'),
		resolution = 500,
		width = 8,
		height = 4,
		xlimits = auto.axis(myData$x, log.scaled = FALSE)$at[c(1, length(auto.axis(myData$x, log.scaled = FALSE)$at))],
		xat = auto.axis(myData$x, log.scaled = FALSE)$at,
		ylimits = c(-0.1,8),
		yat = seq(0,8,2),
		col = myData$colour,
		alpha = 0.65,
		cex = 0.75,
		abline.v = c(
			centromeres[which(centromeres[,1] == paste0('chr',my.chr)),][1,2]/10**6,
			centromeres[which(centromeres[,1] == paste0('chr',my.chr)),][2,3]/10**6
			),
		abline.lty = 2,
		add.rectangle = TRUE,
		col.rectangle = 'grey80', #rep_len(c('grey90', 'grey60'), length.out = nrow(signifEvents)),
		alpha.rectangle = 0.5,
		xleft.rectangle = signifEvents$window.start/10**6,
		xright.rectangle = signifEvents$window.end/10**6,
		ytop.rectangle = 10,
		ybottom.rectangle = -1,
		xaxis.tck = c(0.5,0),
		yaxis.tck = c(0.5,0),
		xlab.label = paste0("Chromosome ", my.chr, ": Genomic Position (Mbp)"),
		ylab.label = parse(text = 'bold("Intermutational Distance"~ "(log"["10"]*")")'),
		xlab.cex = 1,
		ylab.cex = 1,
		xaxis.cex = 1,
		yaxis.cex = 1,
		top.padding = 3,
		right.padding = 2,
		use.legacy.settings = TRUE,
		key = base.change.key,
		add.text = TRUE,
		text.labels = textData$labels,
		text.x = textData$x.pos,
		text.y = rep_len(c(1,0.5), length.out = nrow(textData)),
		text.cex = 0.7
		);
	}

### WRITE SESSION INFO TO FILE ####################################################################
save.session.profile(generate.filename('KataegisPlots', 'SessionInfo', 'txt'));
