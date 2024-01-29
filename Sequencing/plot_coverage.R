### plot_coverage.R ###############################################################################
# Whole exome sequencing was performed at the Broad and Baylor. Raw data was aligned to GRCh38
# using BWA-mem.
# Coverage was calculated using BEDTools with the average taken over all target sequences. 

### LIBRARIES AND SETUP ###########################################################################
library(BoutrosLab.plotting.general);

# set working directory
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/Exomes/config_files/');

# read in config file
cfg <- read.config.file();

### READ DATA  ####################################################################################
# get sample information
phenodata <- read.delim(
	'/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/other/phenodata.txt',
	as.is = TRUE
	);
phenodata[phenodata$SampleName == 'ANPT0021M',]$SampleName <- 'ANPT0021PM';
phenodata[phenodata$SampleName == 'ANPT0145PP',]$SampleName <- 'ANPT0145P';

load('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/other/2018-09-11_ATC_clinical_covariates.RData');

# set data directory
setwd(cfg$cov.data.directory);

# get most recent data file
file.names <- list.files();
file.names <- file.names[grep('ANPT',file.names)];

# read in data
coverage.data <- list();

for (i in file.names) {
	
	sample.name <- unlist(strsplit(as.character(i), '_'))[1];
	raw.data <- tryCatch(
		expr = read.table(i),
		error = function(e) {NULL}
		);
	if (is.null(raw.data)) { next; }
	else {
		if ( grepl('R$|RM$', sample.name)) { cov <- cfg$normal.cov; } else { cov <- cfg$tumour.cov; }
		x <- sum(raw.data[raw.data$V2 >= cov,]$V5);
		coverage.data[[sample.name]] <- x*100;
		}
	}

tmp <- do.call(rbind, coverage.data);

data.to.plot <- data.frame(
	'Sample' = rownames(tmp),
	'Percent' = tmp[,1]
	);

### FORMAT DATA ###################################################################################
data.to.plot <- merge(
	data.to.plot,
	phenodata,
	by.x = 'Sample',
	by.y = 'SampleName'
	);

# fill in gaps in phenodata
tmp.clinical <- atc.clinical$covariates;
tmp.clinical$Sex <- as.character(tmp.clinical$Sex);
tmp.clinical$Age <- as.numeric(factor(tmp.clinical$Age, levels = c('<70', '>=70'), labels = c(0,1)));

for (i in 1:nrow(data.to.plot)) {
	smp <- as.character(data.to.plot[i,]$Sample);
	if (is.na(data.to.plot[i,]$Sex)) {
		data.to.plot[i,]$Sex <- tmp.clinical[which(rownames(tmp.clinical) == smp),]$Sex;
		}
	if (is.na(data.to.plot[i,]$Age)) {
		data.to.plot[i,]$Age <- tmp.clinical[which(rownames(tmp.clinical) == smp),]$Age;
		}
	}

# indicate bar colour
data.to.plot$Colour <- 'grey10';
data.to.plot[data.to.plot$Type == 'Reference',]$Colour <- 'grey60';

# order samples by decreasing coverage (keep patient samples together)
tmp <- data.to.plot[which(data.to.plot$Differentiation != 'Reference'),c('Individual','Sample','Percent')];
tmp <- rbind(tmp, data.to.plot[data.to.plot$Individual == 'ANPT0141',c('Individual','Sample','Percent')]);
tmp <- tmp[order(-tmp$Percent),];
tmp <- droplevels(tmp[!duplicated(tmp$Individual),]);

data.to.plot$Individual <- factor(
	data.to.plot$Individual,
	tmp$Individual
	);

data.to.plot <- data.to.plot[order(data.to.plot$Individual, data.to.plot$Differentiation),];

data.to.plot$Sample <- factor(
	data.to.plot$Sample,
	data.to.plot$Sample
	);

### CREATE PLOT ###################################################################################
# make colour list for sequencing site
sex.colours   <- c('lightskyblue', 'lightpink');
sex.covariate <- as.character(data.to.plot$Sex);
sex.covariate[is.na(sex.covariate)] <- 'grey80';
sex.covariate[sex.covariate == 'M'] <- sex.colours[1];
sex.covariate[sex.covariate == 'F'] <- sex.colours[2];

# fill in age colours
age.colours	<- default.colours(2,'spiral.afternoon');
age.covariate	<- as.character(data.to.plot$Age);
age.covariate[is.na(age.covariate)] <- 'grey80';
age.covariate[which(age.covariate == '2')] <- 'grey80';
age.covariate[which(age.covariate == '0')] <- age.colours[2];
age.covariate[which(age.covariate == '1')] <- age.colours[1];

# get group colours
sampleType.colours <- default.colours(5,'spiral.sunrise');

# match colours to sample type
sampleType.covariate <- as.character(data.to.plot$Sample);

for (i in 1:nrow(data.to.plot)) {
	name <- as.character(data.to.plot$Sample[i]);

	if (grepl('R$|RM$', name))      { x <- sampleType.colours[5]; }
	else if (grepl('C$|BC$',name))  { x <- sampleType.colours[4]; }
	else if (grepl('PW$|PI$',name)) { x <- sampleType.colours[3]; }
	else if (grepl('PM',name))	{ x <- sampleType.colours[1]; }
	else                       	{ x <- sampleType.colours[2]; }
	sampleType.covariate[i] <- x;
	}

sample.covariate <- data.frame(
	Sex = sex.covariate,
	Age = age.covariate,
	Type = sampleType.covariate
	);
rownames(sample.covariate) <- data.to.plot$Sample;

heatmap.data <- sample.covariate;
i <- 0;

for (covariate in 1:ncol(sample.covariate)) {

	heatmap.data[,covariate] <- as.numeric(sample.covariate[,covariate]) + i;
	i <- max(heatmap.data[,covariate]);

	}

# get colours
heatmap.colours <- c(); 
for (i in 1:ncol(sample.covariate)) { heatmap.colours <- c(heatmap.colours, levels(sample.covariate[,i])); }

# make a legend
covariate.legend <- list(
	legend = list(
		colours = c(sex.colours), #'darkgrey'),
		labels = c('Male','Female'), #,'N/A'),
		title = 'Sex'
		),
	legend = list(
		colours = c(age.colours, 'grey80'),
		labels = c(expression(''>='70'),expression(''<'70'),'N/A'),
		title = 'Age'
		),
	legend = list(
		colours = sampleType.colours,
		labels = rev(c('Reference','Cell line','DTC','ATC','Metastasis')),
		title = 'Sample Type'
		)
	);

#if (!cfg$type %in% c('WEX','WGS')) { covariate.legend <- covariate.legend[c(1,2,4)]; }

covariate.legends <- legend.grob(
	legends = covariate.legend,
	size = 2,
	title.just = 'left',
	use.legacy.settings = TRUE
	);

# get yaxis label
#if (cfg$type %in% c('WEX','AmpliSeq')) { ylabel <- expression(bold('               Percent target bases covered' >= '50x')); }
#if ('WGS' == cfg$type) { ylabel <- expression(bold('               Percent bases covered' >= '20x')); }
#if ('WGS' == cfg$type) { ylabel <- expression(bold('               Percent bases covered')); }

# make plot
coverage.plot <- create.barplot(
	Percent ~ Sample, 
	data.to.plot, 
	xlab.label = NULL, 
	yaxis.cex = 1.5, 
	ylimits = c(0,100), 
	yat = seq(0,100,20), 
	xaxis.lab = rep('',nrow(data.to.plot)),
	ylab.label = 'Percent target\nbases covered',
	col = data.to.plot$Colour, 
	xaxis.tck = 0,
	yaxis.tck = c(1,0),
	bottom.padding = 5,
	left.padding = -5,
	abline.h = 80,
	abline.lty = 2,
	abline.lwd = 1,
	abline.col = 'red',
	border.lwd = 0.2,
	legend = list(
		inside = list(fun = draw.key, args = list(key = list(
			text = list(lab = paste0(c(cfg$tumour.cov, cfg$normal.cov),'x'), cex = 1.5, font = 2, col = c('grey10','grey50')),
			columns = 2, between = 0.1
			)), x = 0.85, y = 0.98)
		)
	);

# make heatmap
covariate.heatmap <- create.heatmap(
	t(heatmap.data),
	same.as.matrix = TRUE,
	cluster.dimensions = 'none',
	print.colour.key = FALSE,
	colour.scheme = heatmap.colours,
	total.colours = length(heatmap.colours)+1,
	yaxis.lab = NA
	);

create.multipanelplot(
	plot.objects = list(coverage.plot, covariate.heatmap),
	filename = generate.filename(cfg$project.stem, 'GRCh38_Coverage', 'tiff'),
	height = 6,
	width = 12,
	plot.objects.heights = c(4, 1),
	plot.objects.widths = 1,
	layout.width = 1,
	layout.height = 2,
	left.legend.padding = 0,
	right.legend.padding = 0,
	top.legend.padding = 0,
	bottom.legend.padding = 0,
	use.legacy.settings = TRUE,
	legend = list(
		right = list(fun = covariate.legends)
		)
	#resolution = 800
	);


### MAKE THE PLOT!! ###
setwd(cfg$plot.directory);
if (cfg$type %in% c('WEX', 'AmpliSeq')) {

	# combine plots
	create.multiplot(
		plot.objects = list(covariate.heatmap, coverage.plot),
		filename = generate.filename(cfg$project.stem, 'GRCh38_Coverage', 'tiff'),
		height = 6,
		width = 12,
		panel.heights = c(5,1),
		y.relation = 'free',
		yaxis.cex = 1.5,
		xaxis.tck = 0,
		y.spacing = 1,
		ylab.label = expression(bold('               Percent target\nbases covered')),
		ylab.cex = 2.5,
		ylab.padding = 10,
		print.new.legend = TRUE,
		use.legacy.settings = TRUE,
		legend = list(
			right = list(fun = covariate.legends),
			inside = list(fun = draw.key, args = list(key = list(
				text = list(lab = paste0(c(cfg$tumour.cov, cfg$normal.cov),'x'), cex = 1.5, font = 2, col = c('grey10','grey50')),
				columns = 2, between = 0.1
				)), x = 0.85, y = 0.98)
			),
		resolution = 800
		);
	}

if (cfg$type == 'WGS') {

	ylabel <- expression(bold('               Percent bases covered'));

	# combine plots
	create.multiplot(
		plot.objects = list(covariate.heatmap, coverage.plot),
		filename = generate.filename(cfg$project.stem, 'GRCh38_Coverage', 'tiff'),
		height = 6,
		width = 10,
		panel.heights = c(5,1),
		y.relation = 'free',
		yaxis.cex = 1,
		xaxis.tck = 0,
		y.spacing = 0.75,
		ylab.label = ylabel,
		ylab.cex = 1.5,
		ylab.padding = 8,
		print.new.legend = TRUE,
		use.legacy.settings = TRUE,
		legend = list(
			right = list(fun = covariate.legends),
			inside = list(fun = draw.key, args = list(key = list(
				text = list(lab = c('50x','30x'), cex = 1.5, font = 2, col = c('grey10','grey50')),
				columns = 2, between = 0.1
				)), x = 0.01, y = 0.98)
			),
		resolution = 800
		);
	}

write.table(
	data.to.plot,
	file = generate.filename(cfg$project.stem, '50xCoverage','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SessionInfo ###################################################################################
save.session.profile(generate.filename('SessionInfo', 'Coverage','txt'));
