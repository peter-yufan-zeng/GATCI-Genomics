### plot_contamination.R ##########################################################################
# Whole exome sequencing was performed. Raw data was aligned to GRCh38 using BWA-mem.
# GATK was used to realign/recalibrate indels and to call germline variants.
# ContEst was used to estimate cross-sample contamination. 

### LIBRARIES AND SETUP ###########################################################################
library(BoutrosLab.plotting.general);

# set working directory
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/Exomes/config_files/');

# read in config file
cfg <- read.config.file();

### READ DATA  ####################################################################################
# get sample information
phenodata <- read.delim(cfg$phenodata, as.is = TRUE);
phenodata[phenodata$SampleName == 'ANPT0021M',]$SampleName <- 'ANPT0021PM';
phenodata[phenodata$SampleName == 'ANPT0145PP',]$SampleName <- 'ANPT0145P';

load('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/other/2018-09-11_ATC_clinical_covariates.RData');

# filter by data type
this.col <- colnames(phenodata)[grep(cfg$type, colnames(phenodata))];
phenodata <- phenodata[!is.na(phenodata[,this.col]),];

# remove normal-only sample
phenodata <- droplevels(phenodata[!phenodata$SampleName == 'ANPT0141R',]); 

# only want paired samples
paired.samples <- unique(as.character(phenodata[phenodata$Type == 'Reference',]$Individual));

# set data directory
setwd(cfg$gatk.directory);

# create an empty list to hold data
contest.data <- list();

for (sample.id in paired.samples) {

	these.dirs <- list.files(pattern = sample.id);

	for (dir in these.dirs) {

		these.files <- list.files(path = paste0(dir,'/contest/'), pattern = '.tsv', recursive = TRUE);

		for (file in these.files) {

			type	<- if (grepl('per_bam', file)) { 'META' } else { 'LANE' }
			id	<- sub('contest_','',unlist(strsplit(file,'/'))[2]);
			raw.data <- tryCatch(
				expr = read.delim(paste0(dir,'/contest/',file)),
				error = function(e) { NULL }
				);

			if (is.null(raw.data)) {
				raw.data <- tryCatch(
					expr = read.delim(paste0(dir,'/contest/',file), skip = 1),
					error = function(e) { NULL }
					);
				}

			contest.data[[type]][[id]] <- if (!is.null(raw.data)) { raw.data } else { NA }
			}
		}
	}

### FORMAT DATA ###################################################################################
# format meta data
tmp.meta <- do.call(rbind, contest.data[['META']]);

meta.data <- data.frame(
	'Sample' = rownames(tmp.meta),
	'Percent' = tmp.meta[,4],
	'Conf.low' = tmp.meta[,6],
	'Conf.high' = tmp.meta[,7]
	);

meta.data <- merge(
	meta.data,
	phenodata[phenodata$Individual %in% paired.samples,],
	by.x = 'Sample',
	by.y = 'SampleName',
	all.y = TRUE
	);

# order patients by contamination
tmp <- meta.data[which(meta.data$Differentiation != 'Reference'),c('Individual','Sample','Percent')];
tmp <- rbind(tmp, meta.data[meta.data$Individual == 'ANPT0141',c('Individual','Sample','Percent')]);
tmp <- droplevels(tmp[!duplicated(tmp$Individual),]);
tmp <- tmp[order(-tmp$Percent),];

meta.data$Individual <- factor(
	meta.data$Individual,
	tmp$Individual
	);

meta.data <- meta.data[order(meta.data$Individual, meta.data$Differentiation),];

meta.data$Sample <- factor(
	meta.data$Sample,
	meta.data$Sample
	);

# format lane data
tmp.lane <- do.call(rbind, contest.data[['LANE']]);
tmp.lane$Sample <- sapply(rownames(tmp.lane), function(i) { unlist(strsplit(as.character(i), '\\.'))[1] } );

lane.data <- data.frame(
	'Sample' = tmp.lane$Sample,
	'Lane' = tmp.lane[,1],
	'Percent' = tmp.lane[,4],
	'Conf.low' = tmp.lane[,6],
	'Conf.high' = tmp.lane[,7]
	);

lane.data <- merge(
	lane.data,
	phenodata[phenodata$Individual %in% paired.samples,c('SampleName','Individual')],	
	by.x = 'Sample',
	by.y = 'SampleName',
	all.y = TRUE
	);

lane.data$Sample <- factor(
	lane.data$Sample,
	meta.data$Sample
	);

lane.data <- lane.data[order(lane.data$Sample),];

# fill in gaps in phenodata
tmp.clinical <- atc.clinical$covariates;
tmp.clinical$Sex <- as.character(tmp.clinical$Sex);
tmp.clinical$Age <- as.numeric(as.character(factor(tmp.clinical$Age, levels = c('<70', '>=70'), labels = c('0','1'))));

for (i in 1:nrow(meta.data)) {
        smp <- as.character(meta.data[i,]$Sample);
        if (is.na(meta.data[i,]$Sex)) {
                meta.data[i,]$Sex <- tmp.clinical[which(rownames(tmp.clinical) == smp),]$Sex;
                }
        if (is.na(meta.data[i,]$Age)) {
                meta.data[i,]$Age <- tmp.clinical[which(rownames(tmp.clinical) == smp),]$Age;
                }
        }

### CREATE PLOT ###################################################################################
# make colour list for sequencing site
sex.colours   <- c('lightskyblue', 'lightpink');
sex.covariate <- as.character(meta.data$Sex);
sex.covariate[is.na(sex.covariate)] <- 'grey80';
sex.covariate[sex.covariate == 'M'] <- sex.colours[1];
sex.covariate[sex.covariate == 'F'] <- sex.colours[2];

# fill in age colours
age.colours     <- default.colours(2,'spiral.afternoon');
age.covariate   <- as.character(meta.data$Age);
age.covariate[is.na(age.covariate)] <- 'grey80';
age.covariate[which(age.covariate == '0')] <- age.colours[2];
age.covariate[which(age.covariate == '1')] <- age.colours[1];

# get group colours
sampleType.colours <- default.colours(5,'spiral.sunrise');

# match colours to sample type
sampleType.covariate <- as.character(meta.data$Sample);

for (i in 1:nrow(meta.data)) {
	name <- as.character(meta.data$Sample[i]);

	if (grepl('R$|RM$', name))      { x <- sampleType.colours[5]; }
	else if (grepl('C$|BC$',name))  { x <- sampleType.colours[4]; }
	else if (grepl('PW$|PI$',name)) { x <- sampleType.colours[3]; }
	else if (grepl('PM',name))	{ x <- sampleType.colours[1]; }
	else		       	{ x <- sampleType.colours[2]; }
	sampleType.covariate[i] <- x;
	}


# convert covariates to a heatmap
sample.covariate <- data.frame(
	Sex = sex.covariate,
	Age = age.covariate,
	Type = sampleType.covariate
	);
rownames(sample.covariate) <- meta.data$Sample;

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
		colours = c(sex.colours), #, 'darkgrey'),
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

covariate.legends <- legend.grob(
	legends = covariate.legend,
	size = 2,
	title.just = 'left',
	use.legacy.settings = TRUE
	);

# set axes parameters
max.percent <- max(c(meta.data$Percent, lane.data$Percent));

# make meta plot
meta.plot <- create.scatterplot(
	Percent ~ Sample,
	meta.data, 
	xlab.label = NULL, 
	yaxis.cex = 1.5, 
	ylimits = if ((max.percent > 5) & (max.percent < 10)) { c(0,8.5) } else if (max.percent < 5) { c(0,5) }, 
	yat = if ((max.percent > 5) & (max.percent < 10)) { seq(0,8,2) } else if (max.percent < 5) { seq(0,5,1) },
	xaxis.lab = rep('',nrow(meta.data)),
	xaxis.tck = 0,
	yaxis.tck = c(0.5,0),
	cex = 0.5,
	add.line.segments = TRUE,
	line.end = list(meta.data$Conf.high),
	line.start = list(meta.data$Conf.low),
	abline.h = 3,
	abline.col = "red"
	);

# make lane plot
lane.plot <- create.boxplot(
	Percent ~ Sample,
	lane.data,
	xlab.label = NULL, 
	yaxis.cex = 1.5, 
	ylimits = if ((max.percent > 5) & (max.percent < 10)) { c(0,8.5) } else if (max.percent < 5) { c(0,5) }, 
	yat = if ((max.percent > 5) & (max.percent < 10)) { seq(0,8,2) } else if (max.percent < 5) { seq(0,5,1) },
	xaxis.lab = rep('',nrow(meta.data)),
	xaxis.tck = 0,
	yaxis.tck = c(0.5,0),
	symbol.cex = 0.5,
	abline.h = 3,
	abline.col = "red"	
	);

# make heatmap
covariate.heatmap <- create.heatmap(
	t(heatmap.data),
	same.as.matrix = TRUE,
	cluster.dimensions = 'none',
	print.colour.key = FALSE,
	colour.scheme = heatmap.colours,
	total.colours = length(heatmap.colours)+1,
	yaxis.lab = NA,
	xaxis.lab = NA
	);

# combine plots
setwd(cfg$plot.directory);

create.multiplot(
	plot.objects = list(covariate.heatmap, lane.plot, meta.plot),
	filename = generate.filename(cfg$project.stem, 'ContaminationEstimate', 'tiff'),
	height = 6.5,
	width = 12,
	panel.heights = c(4,2,1),
	y.relation = 'free',
	yaxis.cex = 1.5,
	xaxis.cex = 0.5,
	xaxis.tck = 0,
	xaxis.rot = 45,
	xat = which(meta.data$Percent >= 3),
	xaxis.lab = meta.data[which(meta.data$Percent >= 3),]$Sample,
	y.spacing = 0.8,
	ylab.label = 'ContaminationEstimate (%)',
	ylab.cex = 2,
	ylab.padding = 8,
	top.padding = 2,
	print.new.legend = TRUE,
	use.legacy.settings = TRUE,
	legend = list(
		right = list(fun = covariate.legends),
		inside = list(fun = draw.key, args = list(key = list(text = list(lab = parse(text = "bold('Sample level')"), fontface = 'bold', cex = 1.5))), x = 0.8, y = 0.98),
		inside = list(fun = draw.key, args = list(key = list(text = list(lab = parse(text = "bold('Lane level')"), fontface = 'bold', cex = 1.5))), x = 0.8, y = 0.41)
		)
	#resolution = 800
	);

write.table(
	meta.data,
	file = generate.filename(cfg$project.stem, 'ContaminationEstimates_META', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

write.table(
	lane.data,
	file = generate.filename(cfg$project.stem, 'ContaminationEstimates_LANE', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SessionInfo ###################################################################################
save.session.profile(generate.filename('SessionInfo', 'ContEst','txt'));
