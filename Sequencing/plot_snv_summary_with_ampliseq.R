### plot_snv_summary_with_ampliseq.R ###############################################################
# Purpose: Whole-exome sequencing was performed on ATC tumour and normal samples to identify genes 
#	containing recurrent somatic SNVs in tumour samples.
#	SNVs were called using MuTect (or somaticSniper) and resulting vcfs were run through recSNV
#	to generate basechange data, functional mutation data and a recurrent mutation list.
#	AmpliSeq validation was used to filter FPs.
# Initial write: 2016-02-03
# Modified:	 2016-05-25
# Update: 	 2017-10-27

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general);

# general parameters
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/Exomes/plots/');

# general parameters
date <- Sys.Date();
cfg <- read.config.file();

### READ IN DATA FROM RECSNV ######################################################################
# get sample annotation
phenodata <- read.delim('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/other/phenodata.txt');

# read in basechange summary
basechangeSummary <- read.delim(gsub('filtered_pon', 'ampliseq_filtered', cfg$basechange.data));

# read in functional change summary
functionalSummary <- read.delim(gsub('filtered_pon', 'ampliseq_filtered', cfg$functional.data));

# reorder variable levels to match heatmap data
functionalSummary$variable <- factor(
	functionalSummary$variable,
	levels = c('Nonsynonymous', 'Stopgain-stoploss', 'Splicing', 'UTR')
	);
functionalSummary <- functionalSummary[!is.na(functionalSummary$variable),]

# read in recurrence matrix
recurrenceData <- read.delim(gsub('filtered_pon', 'ampliseq_filtered', cfg$full.recurrence));
recurrenceData[recurrenceData == 5] <- NA;
recurrenceData <- recurrenceData[,!grepl('ANPT0173', colnames(recurrenceData))];

# get mutation counts per patient
mutationRates <- read.delim('2017-09-19_ATC_mutation_rates.tsv');
mutationRates <- mutationRates[!grepl('ANPT0173', mutationRates$Sample),];

# load covariate/clinical info
load('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/other/ATC_clinical_covariates.RData');
 
### ORGANIZE DATA #################################################################################
# clean up phenodata
phenodata <- phenodata[!is.na(phenodata$WEX.Site),];
phenodata <- phenodata[!phenodata$Type == 'Reference',];

phenodata$ID <- as.character(phenodata$Individual);
phenodata[phenodata$SampleName == 'ANPT0012PP',]$ID <- 'ANPT0012PP';
phenodata[phenodata$SampleName == 'ANPT0012PW',]$ID <- 'ANPT0012PW';
phenodata[phenodata$SampleName == 'ANPT0015PP',]$ID <- 'ANPT0015P';
phenodata[phenodata$SampleName == 'ANPT0015PW',]$ID <- 'ANPT0015W';
phenodata[phenodata$SampleName == 'ANPT0021P', ]$ID <- 'ANPT0021';
phenodata[phenodata$SampleName == 'ANPT0021M',]$ID <- 'ANPT0021M';
phenodata[phenodata$SampleName == 'ANPT0028PP',]$ID <- 'ANPT0028P';
phenodata[phenodata$SampleName == 'ANPT0028PW',]$ID <- 'ANPT0028W';
phenodata[phenodata$SampleName == 'ANPT0127PW',]$ID <- 'ANPT0127PW';
phenodata[phenodata$SampleName == 'ANPT0134PP',]$ID <- 'ANPT0134PP';
phenodata[phenodata$SampleName == 'ANPT0134PW',]$ID <- 'ANPT0134PW';
phenodata[phenodata$SampleName == 'ANPT0135PP',]$ID <- 'ANPT0135PP';
phenodata[phenodata$SampleName == 'ANPT0135PW',]$ID <- 'ANPT0135PW';
phenodata[phenodata$SampleName == 'ANPT0137PP',]$ID <- 'ANPT0137P';
phenodata[phenodata$SampleName == 'ANPT0137PW',]$ID <- 'ANPT0137W';
phenodata[phenodata$SampleName == 'ANPT0139PP',]$ID <- 'ANPT0139P';
phenodata[phenodata$SampleName == 'ANPT0139PW',]$ID <- 'ANPT0139W';
phenodata[phenodata$SampleName == 'ANPT0143PW',]$ID <- 'ANPT0143PW';
phenodata[phenodata$SampleName == 'ANPT0143PI',]$ID <- 'ANPT0143PI';
phenodata[phenodata$SampleName == 'ANPT0144PP',]$ID <- 'ANPT0144PP';
phenodata[phenodata$SampleName == 'ANPT0144PW',]$ID <- 'ANPT0144PW';
phenodata[phenodata$SampleName == 'ANPT0147PP',]$ID <- 'ANPT0147P';
phenodata[phenodata$SampleName == 'ANPT0147PW',]$ID <- 'ANPT0147W';
phenodata[phenodata$SampleName == 'ANPT0148PP',]$ID <- 'ANPT0148P';
phenodata[phenodata$SampleName == 'ANPT0148PW',]$ID <- 'ANPT0148W';
phenodata[phenodata$SampleName == 'ANPT0155PP',]$ID <- 'ANPT0155P';
phenodata[phenodata$SampleName == 'ANPT0155PW',]$ID <- 'ANPT0155W';
phenodata[phenodata$SampleName == 'ANPT0160PI',]$ID <- 'ANPT0160I';
phenodata[phenodata$SampleName == 'ANPT0161P',]$ID <- 'ANPT0161P';
phenodata[phenodata$SampleName == 'ANPT0161PI',]$ID <- 'ANPT0161PI';
phenodata[phenodata$SampleName == 'ANPT0165PI',]$ID <- 'ANPT0165PI';
phenodata[phenodata$SampleName == 'ANPT0177PP',]$ID <- 'ANPT0177P';
phenodata[phenodata$SampleName == 'ANPT0177PW',]$ID <- 'ANPT0177W';
phenodata[phenodata$SampleName == 'ANPT0179PP',]$ID <- 'ANPT0179P';
phenodata[phenodata$SampleName == 'ANPT0179PW',]$ID <- 'ANPT0179W';
phenodata[phenodata$SampleName == 'ANPT0182PP',]$ID <- 'ANPT0182PP';
phenodata[phenodata$SampleName == 'ANPT0182PW',]$ID <- 'ANPT0182PW';

# remove samples with no SNV data
phenodata <- phenodata[which(phenodata$ID %in% colnames(recurrenceData)),];

# add clinical
phenodata <- merge(phenodata[,c(1,2,14,9:13)], atc.clinical$covariates, by.x = 'SampleName', by.y = 'row.names');

# indicate sample groups
cell.line <- phenodata[phenodata$Type == 'CellLine',]$ID;
non.atc <- phenodata[grep('W$|M$|I$', phenodata$ID),]$ID;
primary <- phenodata[phenodata$Type == 'ATC',]$ID;

# using only paired data, find recurrently altered genes
recurrenceCount <- data.frame(
	Gene  = rownames(recurrenceData),
	Count = apply(recurrenceData[,primary], 1, function(i) { length(i[!is.na(i)]) }),
	ATC.Prop = apply(recurrenceData[,primary], 1, function(i) { length(i[!is.na(i)])/length(i) }),
	Cell.Prop = apply(recurrenceData[,cell.line], 1, function(i) { length(i[!is.na(i)])/length(i) })
	);
recurrenceCount <- recurrenceCount[order(-recurrenceCount$Count),];

# filter out genes primarily mutated in contaminated samples
rec.genes <- as.character(recurrenceCount[recurrenceCount$Count > 11,]$Gene);

recurrenceCount <- recurrenceCount[recurrenceCount$Gene %in% rec.genes,];
recurrenceCount <- recurrenceCount[order(-recurrenceCount$ATC.Prop),];

alt.genes <- as.character(recurrenceCount[grepl('^MUC|TTN', recurrenceCount$Gene),]$Gene);
keep.genes <- as.character(recurrenceCount[!grepl('^MUC|TTN', recurrenceCount$Gene),]$Gene);

# remove non-recurrent genes
recurrenceData  <- recurrenceData[keep.genes,];

# set gene order
recurrenceCount <- recurrenceCount[keep.genes,];
recurrenceCount$Gene <- factor(
	recurrenceCount$Gene,
	levels = rev(recurrenceCount$Gene)
	);

# sort samples by occurrence
tmpData <- t(recurrenceData);
tmpData[!is.na(tmpData)] <- 1;
tmpData <- tmpData[do.call(order,transform(tmpData)),];

# get mutation counts
sampleCount <- data.frame(
	Sample = mutationRates$SampleID,
	Count = log10(mutationRates$MutsPerMb)
	);

sampleCount$Sample <- factor(
	sampleCount$Sample,
	levels = rownames(tmpData)
	);
sampleCount <- sampleCount[order(sampleCount$Sample),];

# ensure same sample ordering
phenodata$ID <- factor(
	phenodata$ID,
	levels = sampleCount$Sample
	);
phenodata <- phenodata[order(phenodata$ID),];

basechangeSummary$Sample <- factor(
	basechangeSummary$Sample,
	levels = sampleCount$Sample
	);
functionalSummary$Sample <- factor(
	functionalSummary$Sample,
	levels = sampleCount$Sample
	);
recurrenceData <- recurrenceData[,as.character(sampleCount$Sample)];

### VISUALIZATION #################################################################################
## make covariates
colour.data <- atc.clinical$covariate.data[as.character(phenodata$SampleName),];
rownames(colour.data) <- phenodata$ID;

# don't show age/event covariates for cell lines
colour.data[cell.line,c('Age','Event')] <- 'white';

# trim covariates
na.counts <- apply(colour.data,2,function(i) { length(i[grepl('grey80',i)])/length(i) } );
keep.covariates <- names(na.counts[which(na.counts < 0.5)]);

colour.data <- colour.data[,keep.covariates];
heatmap.data <- colour.data;
heatmap.colours <- c();

i <- 0;
for (covariate in 1:ncol(heatmap.data)) {
	heatmap.data[,covariate] <- as.numeric(as.factor(colour.data[,covariate])) + i;
	heatmap.colours <- c(heatmap.colours, levels(as.factor(colour.data[,covariate])));
	i <- max(heatmap.data[,covariate]);
	}

# update legend
clinical.legend <- atc.clinical$clinical.legend;
clinical.legend[[1]]$colours <- clinical.legend[[1]]$colours[1:4];
clinical.legend[[1]]$labels <- c('Metastasis','ATC','W/PDTC', 'Cell line');

clinical.legend.trim <- list(
	legend = clinical.legend[[1]],
	legend = clinical.legend[[2]],
	legend = clinical.legend[[3]],
	legend = clinical.legend[[8]]
	);

# make the heatmaps
atc.covariate <- create.heatmap(
	t(heatmap.data[rownames(heatmap.data) %in% c(primary, non.atc),]),
	same.as.matrix = TRUE,
	cluster.dimensions = 'none',
	print.colour.key = FALSE,
	colour.scheme = heatmap.colours,
	total.colours = length(heatmap.colours)+1,
	at = 0:length(heatmap.colours)
	);

symbol.locations <- list(
	  borders = list(
		  list(   
			  xright = 2.05, # bottom
			  xleft = 0.98, # top
			  ybottom = 1, # left
			  ytop = 13, # right
			  size = 4,
			  col = 'black'
			  ),
		  list(   
			  xright = 2.05, # bottom
			  xleft = 0.98, # top
			  ybottom = 1, # left
			  ytop = 13, # right
			  size = 4,
			  col = 'black'
			  )
		  )
	  );

cell.covariate <- create.heatmap(
	t(heatmap.data[rownames(heatmap.data) %in% cell.line,]),
	#filename = 'test.tiff', resolution = 500,
	same.as.matrix = TRUE,
	cluster.dimensions = 'none',
	print.colour.key = FALSE,
	colour.scheme = heatmap.colours,
	total.colours = length(heatmap.colours)+1,
	at = 0:length(heatmap.colours),
	axes.lwd = 0,
	symbols = symbol.locations
	);

# create colour schemes
# basechange
basechangeColours	 <- force.colour.scheme(1:nlevels(basechangeSummary$variable), 'chromosome');
names(basechangeColours) <- levels(basechangeSummary$variable);

# mutational consequence
functionalColours	 <- c('darkseagreen4', 'darkturquoise', 'gold1') #, 'orchid4');
names(functionalColours) <- c('Nonsynonymous', 'Stopgain/Stoploss', 'Splicing') #, 'UTR');

# create legend/colour key
seq.legend <- list(
	legend = list(
		colours = basechangeColours,
		labels = c(
			expression("A"%->%"C / T"%->%"G"),
			expression("A"%->%"G / T"%->%"C"),
			expression("A"%->%"T / T"%->%"A"),
			expression("G"%->%"A / C"%->%"T"),
			expression("G"%->%"C / C"%->%"G"),
			expression("G"%->%"T / C"%->%"A")
			),
		title = 'Basechange'
		),
	legend = list(
		colours = functionalColours,
		labels = names(functionalColours),
		title = 'Function'
		)
	);

seqKey <- legend.grob(
	legends = seq.legend,
	size = 1.5,
	label.cex = 0.8,
	title.cex = 0.8,
	title.just = 'left',
	layout = c(2,1),
	use.legacy.settings = TRUE
	);

clinicalKey <- legend.grob(
	legends = clinical.legend.trim,
	size = 1.5,
	label.cex = 0.8,
	title.cex = 0.8,
	title.just = 'left',
	layout = c(length(clinical.legend.trim),1),
	use.legacy.settings = TRUE
	);

## create individual plots:
# ATC tumours only
# create basechange barplot
basechangePlot <- create.barplot(
	proportion ~ Sample,
	basechangeSummary[basechangeSummary$Sample %in% c(primary, non.atc),],
	groups = basechangeSummary[basechangeSummary$Sample %in% c(primary, non.atc),]$variable,
	stack = TRUE,
	col = basechangeColours,
	border.col = 'transparent'
	);

# create functional consequence barplot
functionalPlot <- create.barplot(
	proportion ~ Sample,
	functionalSummary[functionalSummary$Sample %in% c(primary, non.atc),],
	groups = functionalSummary[functionalSummary$Sample %in% c(primary, non.atc),]$variable,
	stack = TRUE,
	col = functionalColours,
	border.col = 'transparent'
	);

# create mutation count plot
mutCountPlot <- create.scatterplot(
	Count ~ Sample,
	sampleCount[sampleCount$Sample %in% c(primary, non.atc),],
	cex = 0.5
	);

# create recurrence count barplot
genecountPlot <- create.barplot(
	Gene ~ ATC.Prop,
	recurrenceCount,
	xlimits = c(0,1),
	xat = seq(0,1,0.5),
	plot.horizontal = TRUE
	);

# create heatmap of sample x gene matrix (no clustering)
recurrenceHeatmap <- create.heatmap(
	t(recurrenceData[rev(as.character(recurrenceCount$Gene)), colnames(recurrenceData) %in% c(primary, non.atc)]),
	cluster.dimensions = 'none',
	yaxis.lab = NA,
	xaxis.lab = NA,
	colour.scheme = functionalColours,
	total.colours = 5,
	grid.row = TRUE,
	grid.col = TRUE,
	force.grid.row = TRUE,
	force.grid.col = TRUE,
	row.lwd = 2,
	col.lwd = 2,
	row.colour = 'white',
	col.colour = 'white',
	print.colour.key = TRUE,
	fill.colour = 'grey95',
	at = seq(0,4,1)
	);

cellProp <- create.barplot(
	Gene ~ Cell.Prop,
	recurrenceCount,
	plot.horizontal = TRUE
	);

cellHeatmap <- create.heatmap(
	t(recurrenceData[rev(as.character(recurrenceCount$Gene)), colnames(recurrenceData) %in% cell.line]),
	cluster.dimensions = 'none',
	colour.scheme = functionalColours,
	total.colours = 5,
	grid.row = TRUE,
	grid.col = TRUE,
	force.grid.row = TRUE,
	force.grid.col = TRUE,
	row.lwd = 2,
	col.lwd = 2,
	row.colour = 'white',
	col.colour = 'white',
	print.colour.key = FALSE,
	fill.colour = 'grey95',
	at = seq(0,4,1)
	);

### CREATE MULTIPLOT ##############################################################################
# make some lists!
myPlots <- list(
	atc.covariate, cell.covariate,
	recurrenceHeatmap, genecountPlot, 
	cellHeatmap, cellProp, 
	#functionalPlot,
	basechangePlot, mutCountPlot
	);

myAxes <- list(
	ylimits = list(
		NULL, NULL,
		NULL, NULL,
		NULL, NULL,
		c(0,1), c(-1,3.5)
		),
	xlimits = list(
		NULL, NULL,
		NULL, c(0,1),
		NULL, c(0,1),
		NULL, NULL
		),
	yat = list(
		1:ncol(heatmap.data), 1:ncol(heatmap.data),
		1:nrow(recurrenceData), 1:nrow(recurrenceCount), 
		1:nrow(recurrenceData), 1:nrow(recurrenceCount), 
		seq(0, 1, 0.5),	seq(-1, 3, 1)	
		),
	xat = list(
		1:nrow(heatmap.data[c(primary, non.atc),]), 1:nrow(heatmap.data[cell.line,]),
		1:ncol(recurrenceData[,c(primary, non.atc)]), seq(0, 1, 0.5),
		1:ncol(recurrenceData[,cell.line]), seq(0, 1, 0.5),
		1:ncol(recurrenceData[,c(primary, non.atc)]), 1:ncol(recurrenceData[,c(primary, non.atc)])
		),
	yaxis.lab = list(
		rev(colnames(heatmap.data)), rep('', ncol(heatmap.data)),
		rev(rownames(recurrenceData)), rep('', nrow(recurrenceCount)), 
		rep('', nrow(recurrenceData)), rep('', nrow(recurrenceCount)),
		seq(0, 1, 0.5), c(0.1,1,10,100,1000)
		),
	xaxis.lab = list(
		rep('', nrow(heatmap.data[c(primary, non.atc),])), rep('', nrow(heatmap.data[cell.line,])),
		rep('', ncol(recurrenceData[,c(primary, non.atc)])), seq(0, 1, 0.5),
		rep('', ncol(recurrenceData[,cell.line])), seq(0, 1, 0.5),
		rep('', ncol(recurrenceData[,c(primary, non.atc)])), rep('', ncol(recurrenceData[,c(primary, non.atc)]))
		)
	);

# combine individual plots
create.multiplot(
	plot.objects = myPlots,
	filename = generate.filename('AmpliSeq_filtered', 'RecurrentSNV_Profile', 'tiff'),
	plot.layout = c(4,4),
	layout.skip = c(FALSE,TRUE,FALSE,TRUE, rep(FALSE, 4), FALSE, rep(TRUE, 3), FALSE, rep(TRUE, 3)),
	panel.widths = c(10, 1, 1, 1),
	panel.heights = c(0.8, 1, 4.5, 0.6),
	resolution = 500, #cfg$resolution,
	height = 8,
	width = 11,
	x.spacing = -0.5,
	y.spacing = -0.5,
	xaxis.cex = 0.75,
	yaxis.cex = 0.6,
	ylab.label = '\t\t\t\t\t\t\t\t\t\t\t\t\t\tProportion   SNVs/Mbp    ',
	ylab.cex = 1,
	xlab.cex = 0.75,
	y.relation = 'free',
	x.relation = 'free',
	ylimits = myAxes$ylimits,
	xlimits = myAxes$xlimits,
	yat = myAxes$yat,
	xat = myAxes$xat,
	yaxis.lab = myAxes$yaxis.lab,
	xaxis.lab = myAxes$xaxis.lab,
	xaxis.tck = 0.2,
	yaxis.tck = 0.2,
	print.new.legend = TRUE,
	right.padding = 2,
	top.padding = 2,
	bottom.padding = 8,
	use.legacy.settings = TRUE,
	legend = list(
		inside = list(fun = seqKey, x = 0.76, y = 0.95),
		inside = list(fun = clinicalKey, x = 0, y = -0.03),
		inside = list(fun = draw.key, args = list(key = list(text = list(lab = parse(text = "bold('Proportion')"), cex = 0.75))), x = 0.752, y = 0.06),
		inside = list(fun = draw.key, args = list(key = list(text = list(lab = parse(text = "bold('Proportion')"), cex = 0.75))), x = 0.922, y = 0.06)
		)
	);

### WRITE SESSION INFO TO FILE ####################################################################
save.session.profile(generate.filename('VariantHeatmap_AmpliSeqFiltered', 'SessionInfo', 'txt'));
