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
phenodata <- read.delim(cfg$platform.map);

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
mutationRates <- read.delim('2018-03-19_ATC_MutationRates.tsv');
mutationRates <- mutationRates[!grepl('ANPT0173', mutationRates$Sample),];

# get SeqSig results
seqsig <- read.delim(cfg$seqsig.results);

# load covariate/clinical info
load('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/other/ATC_clinical_covariates.RData');

### ORGANIZE DATA #################################################################################
# add clinical
phenodata <- merge(
	phenodata[,c(1,4,6)],
	atc.clinical$covariates,
	by.x = 'CommonName',
	by.y = 'row.names'
	);

# clean up phenodata
phenodata <- phenodata[phenodata$WXS.Name %in% names(recurrenceData),];

# indicate sample groups
cell.line <- intersect(colnames(recurrenceData), phenodata[phenodata$Type == 'CellLine',]$WXS.Name);
non.atc <- intersect(colnames(recurrenceData), phenodata[grep('W$|M$|I$', phenodata$WXS.Name),]$WXS.Name);
primary <- intersect(colnames(recurrenceData), phenodata[phenodata$Type == 'ATC',]$WXS.Name);

# using only paired data, find recurrently altered genes
recurrenceCount <- data.frame(
	Gene  = rownames(recurrenceData),
	Count = apply(recurrenceData[,primary], 1, function(i) { length(i[!is.na(i)]) }),
	ATC.Prop = apply(recurrenceData[,primary], 1, function(i) { length(i[!is.na(i)])/length(i) }),
	Cell.Prop = apply(recurrenceData[,cell.line], 1, function(i) { length(i[!is.na(i)])/length(i) })
	);
recurrenceCount <- recurrenceCount[order(-recurrenceCount$Count),];

# filter out genes primarily mutated in contaminated samples
rec.genes <- as.character(recurrenceCount[recurrenceCount$Count >= 11,]$Gene);

recurrenceCount <- recurrenceCount[recurrenceCount$Gene %in% rec.genes,];
recurrenceCount <- recurrenceCount[order(-recurrenceCount$ATC.Prop),];

alt.genes <- as.character(recurrenceCount[grepl('^MUC|TTN', recurrenceCount$Gene),]$Gene);
keep.genes <- as.character(recurrenceCount[!grepl('^MUC|TTN', recurrenceCount$Gene),]$Gene);

sig.genes <- intersect(keep.genes, rownames(seqsig[which(seqsig$fdr.adjusted.p < 0.01),]));

# remove non-recurrent genes
recurrenceData  <- recurrenceData[sig.genes,];

# set gene order
recurrenceCount <- recurrenceCount[sig.genes,];
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
	Sample = mutationRates[mutationRates$SampleID %in% phenodata$WXS.Name,]$SampleID,
	Count = log10(mutationRates[mutationRates$SampleID %in% phenodata$WXS.Name,]$MutsPerMb)
	);

sampleCount$Sample <- factor(
	sampleCount$Sample,
	levels = rownames(tmpData)
	);
sampleCount <- sampleCount[order(sampleCount$Sample),];

# ensure same sample ordering
phenodata$WXS.Name <- factor(
	phenodata$WXS.Name,
	levels = sampleCount$Sample
	);
phenodata <- phenodata[order(phenodata$WXS.Name),];

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
colour.data <- atc.clinical$covariate.data; 
colour.data <- colour.data[as.character(phenodata$CommonName),];
rownames(colour.data) <- phenodata$WXS.Name;

# don't show age/event covariates for cell lines
colour.data[cell.line,c('Age','Event')] <- 'white';

# trim covariates
na.counts <- apply(colour.data,2,function(i) { length(i[grepl('grey80',i)])/length(i) } );
keep.covariates <- names(na.counts[which(na.counts < 0.4)]);

colour.data <- colour.data[,keep.covariates];
heatmap.data <- colour.data;
heatmap.colours <- c();

i <- 0;
for (covariate in 1:ncol(heatmap.data)) {
	heatmap.data[,covariate] <- as.numeric(as.factor(colour.data[,covariate])) + i;
	heatmap.colours <- c(heatmap.colours, levels(as.factor(colour.data[,covariate])));
	i <- max(heatmap.data[,covariate]);
	}

# update Event name
colnames(heatmap.data)[which(colnames(heatmap.data) == 'Event')] <- 'Status';

# update legend
clinical.legend <- atc.clinical$clinical.legend;
clinical.legend[[1]]$colours <- clinical.legend[[1]]$colours[1:4];
clinical.legend[[1]]$labels <- c('Metastasis','ATC','DTC', 'Cell line');

clinical.legend[[2]]$colours <- clinical.legend[[2]]$colours[1:2];
clinical.legend[[2]]$labels <- clinical.legend[[2]]$labels[1:2];

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
	at = 0:length(heatmap.colours),
	axes.lwd = 1,
	yaxis.lab = NA,
	yaxis.cex = 1
	);

symbol.locations <- list(
	  borders = list(
		  list(   
			  xright = 2.05, # bottom
			  xleft = 0.98, # top
			  ybottom = 1, # left
			  ytop = 13, # right
			  size = 0,
			  col = 'black'
			  ),
		  list(   
			  xright = 2.05, # bottom
			  xleft = 0.98, # top
			  ybottom = 1, # left
			  ytop = 13, # right
			  size = 1,
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
functionalColours	 <- c('darkseagreen4', 'darkturquoise', 'gold1'); #, 'orchid4');
names(functionalColours) <- c('Nonsynonymous', 'Stopgain/Stoploss', 'Splicing'); #, 'UTR');

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
	border.col = 'transparent',
	xaxis.lab = rep('', ncol(recurrenceData[,c(primary, non.atc)])),
	yaxis.lab = seq(0, 1, 0.5),
	xat = 1:ncol(recurrenceData[,c(primary, non.atc)]),
	yat = seq(0, 1, 0.5),
	ylab.label = '\nProportion',
	ylab.cex = 1,
	xlab.label = '',
	xaxis.cex = 0,
	yaxis.cex = 1,
	xaxis.tck = c(0.01,0),
	yaxis.tck = c(0.2,0)
	);

# create functional consequence barplot
#functionalPlot <- create.barplot(
#	proportion ~ Sample,
#	functionalSummary[functionalSummary$Sample %in% c(primary, non.atc),],
#	groups = functionalSummary[functionalSummary$Sample %in% c(primary, non.atc),]$variable,
#	stack = TRUE,
#	col = functionalColours,
#	border.col = 'transparent'
#	);

# create mutation count plot
mutCountPlot <- create.scatterplot(
	Count ~ Sample,
	sampleCount[sampleCount$Sample %in% c(primary, non.atc),],
	cex = 0.5,
	xaxis.lab = rep('', ncol(recurrenceData[,c(primary, non.atc)])),
	yaxis.lab = qw('0.01 0.1 1 10 100'),
	xat = 1:ncol(recurrenceData[,c(primary, non.atc)]),
	ylimits = c(-2,2),
	yat = seq(-2, 2, 1),
	ylab.label = 'SNVs/Mbp\n',
	ylab.cex = 1,
	xlab.label = '',
	xaxis.cex = 0,
	yaxis.cex = 1,
	xaxis.tck = c(0.01,0),
	yaxis.tck = c(0.2,0)
	);

# create recurrence count barplot (no longer used)
genecountPlot <- create.barplot(
	Gene ~ ATC.Prop,
	recurrenceCount,
	xlimits = c(0,1),
	xat = seq(0,1,0.5),
	xaxis.lab = c('0','0.5','1'),
	plot.horizontal = TRUE,
	yaxis.lab = rep('',nrow(recurrenceCount)),
	xaxis.tck = c(0.1,0),
	yaxis.tck = 0,
	xaxis.cex = 0.75,
	ylab.label = NULL,
	xlab.label = 'Proportion',
	xlab.cex = 0.65
	);

# create significance barplot
seqsig$Gene <- rownames(seqsig);
seqsig$Gene <- factor(
	seqsig$Gene,
	levels = rev(recurrenceCount$Gene)
	);
seqsig[which(seqsig$nlog10.fdr.adjusted.p == 'Inf'),]$nlog10.fdr.adjusted.p <- 20;

signifPlot <- create.barplot(
	Gene ~ nlog10.fdr.adjusted.p,
	seqsig[sig.genes,],
	xlimits = c(0,10),
	xat = c(0,2,10), #seq(0,10,2),
	#xaxis.lab = c(expression(bold('1')^phantom('1')), expression(bold('10')^bold('-2')), expression(bold('10')^bold('-8'))), #c('0','0.5','1'),
	xaxis.lab = c(expression(bold('0')),expression(bold('2')),expression(''>=bold('10'))), 
	plot.horizontal = TRUE,
	yaxis.lab = rep('',length(keep.genes)),
	xaxis.tck = c(0.1,0),
	yaxis.tck = 0,
	xaxis.cex = 0.75,
	ylab.label = NULL,
	xlab.label = expression(bold('-log')[bold('10')]*bold(' FDR')),
	xlab.cex = 0.65,
	abline.v = 2,
	abline.col = 'red'
	);

# create heatmap of sample x gene matrix (no clustering)
recurrenceHeatmap <- create.heatmap(
	t(recurrenceData[rev(as.character(recurrenceCount$Gene)), colnames(recurrenceData) %in% c(primary, non.atc)]),
	cluster.dimensions = 'none',
	yaxis.lab = NA,
	xaxis.lab = NULL,
	yaxis.cex = 1,
	colour.scheme = functionalColours,
	total.colours = 4,
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
	at = seq(0,3,1),
	ylab.label = '',
	xlab.label = '',
	axes.lwd = 1,
	yaxis.tck = c(0.1,0)
	);

cellProp <- create.barplot(
	Gene ~ Cell.Prop,
	recurrenceCount,
	plot.horizontal = TRUE,
	xlimits = c(0,1),
	xat = seq(0,1,0.5),
	xaxis.lab = c('0','0.5','1'),
	yaxis.lab = rep('',nrow(recurrenceCount)),
	xaxis.tck = c(0.1,0),
	yaxis.tck = 0,
	xaxis.cex = 0.75,
	ylab.label = NULL,
	xlab.label = 'Proportion',
	xlab.cex = 0.65
	);

cellHeatmap <- create.heatmap(
	t(recurrenceData[rev(as.character(recurrenceCount$Gene)), colnames(recurrenceData) %in% cell.line]),
	cluster.dimensions = 'none',
	colour.scheme = functionalColours,
	total.colours = 4,
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
	at = seq(0,3,1),
	xaxis.lab = NULL,
	xaxis.tck = 0,
	xaxis.top.tck = 0,
	yaxis.tck = c(0.1,0),
	yaxis.lab = NULL,
	axes.lwd = 1
	);

### CREATE MULTIPLOT ##############################################################################
# combine the plots
create.multipanelplot(
	plot.objects = list(mutCountPlot, basechangePlot, recurrenceHeatmap, genecountPlot, cellHeatmap, cellProp, atc.covariate, cell.covariate),
	filename = generate.filename('AmpliSeq_filtered', 'RecurrentSNV_Profile', 'tiff'),
	#resolution = 500,
	height = 7,
	width = 11,
	plot.objects.heights = c(1.2, 1.1, 4, 1),
	plot.objects.widths = c(9.8, 1, 0.8, 1),
	layout.width = 4,
	layout.height = 4,
	layout.skip = c(FALSE, rep(TRUE,3), FALSE, rep(TRUE,3), rep(FALSE,4), rep(c(FALSE,TRUE),2)),
	x.spacing = 1.5,
	y.spacing = -3.5,
	ylab.cex = 1,
	xlab.cex = 0.75,
	use.legacy.settings = TRUE,
	ylab.axis.padding = -3.5,
	left.legend.padding = 0,
	right.legend.padding = 1,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	top.padding = -1,
	legend = list(
		inside = list(fun = seqKey, x = 0.86, y = 0.87),
		bottom = list(fun = clinicalKey)
		)
	);

create.multipanelplot(
	plot.objects = list(mutCountPlot, basechangePlot, recurrenceHeatmap, signifPlot, cellHeatmap, cellProp, atc.covariate, cell.covariate),
	filename = generate.filename('AmpliSeq_filtered', 'RecurrentSNV_Profile_SeqSig', 'tiff'),
	#resolution = 500,
	height = 7,
	width = 11,
	plot.objects.heights = c(1.2, 1.1, 4, 1),
	plot.objects.widths = c(9.8, 1, 0.8, 1),
	layout.width = 4,
	layout.height = 4,
	layout.skip = c(FALSE, rep(TRUE,3), FALSE, rep(TRUE,3), rep(FALSE,4), rep(c(FALSE,TRUE),2)),
	x.spacing = c(1.5,3,2),
	y.spacing = c(-3.5,-3.5,-4.5),
	ylab.cex = 1,
	xlab.cex = 0.75,
	use.legacy.settings = TRUE,
	ylab.axis.padding = -3.5,
	left.legend.padding = 0,
	right.legend.padding = 1,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	top.padding = -1,
	legend = list(
		inside = list(fun = seqKey, x = 0.84, y = 0.85),
		bottom = list(fun = clinicalKey)
		)
	);

### WRITE SESSION INFO TO FILE ####################################################################
save.session.profile(generate.filename('VariantHeatmap_AmpliSeqFiltered', 'SessionInfo', 'txt'));
