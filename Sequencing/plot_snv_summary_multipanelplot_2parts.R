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

# get TERT status
tert.status <- read.delim('~/ATC/other/2019-09-19_tert_status.txt', as.is = TRUE);
tert.status <- tert.status[which(tert.status$WXS == 'Present' & tert.status$Type == 'ATC'),c('ID','Patient_ID','TERT')];
tert.status <- tert.status[!is.na(tert.status$TERT),];

### READ IN DATA FROM RECSNV ######################################################################
# get sample annotation
phenodata <- read.delim(cfg$platform.map);

has.normal <- as.character(phenodata[which(phenodata$Differentiation == 'Reference' & !is.na(phenodata$WXS.Site)),]$Individual);

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

# indicate samples to remove
to.remove <- c('ANPT0173');

# read in recurrence matrix
recurrenceData <- read.delim(gsub('filtered_pon', 'ampliseq_filtered', cfg$full.recurrence));
recurrenceData[recurrenceData == 5] <- NA;
recurrenceData <- recurrenceData[,!colnames(recurrenceData) %in% to.remove];

# get mutation counts per patient
mutationRates <- read.delim('2018-03-19_ATC_MutationRates.tsv');
mutationRates <- mutationRates[!mutationRates$Sample %in% to.remove,];

# get SeqSig results
seqsig <- read.delim(cfg$seqsig.results);

# load covariate/clinical info
load('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/other/2018-09-11_ATC_clinical_covariates.RData');

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
primary <- setdiff(primary, 'ANPT0175');

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

colour.data$TERT <- NA;
for (i in 1:nrow(colour.data)) {
	smp <- rownames(colour.data)[i];
	if (smp %in% tert.status$ID) {
		colour.data[smp,]$TERT <- tert.status[which(tert.status$ID == smp),]$TERT;
		}
	else if (smp %in% tert.status$Patient_ID) {
		colour.data[smp,]$TERT <- tert.status[which(tert.status$Patient_ID == smp),]$TERT;
		}
	else if (nchar(smp) == 9) {
		if (paste0(smp, 'P') %in% tert.status$ID) {
			colour.data[smp,]$TERT <- tert.status[which(tert.status$ID == paste0(smp,'P')),]$TERT;
			}
		}
	}		

tert.colours <- rep('grey80', nrow(colour.data));
tert.colours[which(colour.data$TERT == 0)] <- 'white';
tert.colours[which(colour.data$TERT %in% c('C228T','C250T'))] <- 'darkorange';

colour.data$TERT <- tert.colours;
keep.covariates <- c('TERT','Type','Sex','Age','Event');

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
	legend = list(
		colours = c('darkorange','white','grey80'),
		labels = c('promoter','WT','N/A'),
		title = 'TERT Status'
		),
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
	yaxis.cex = 1,
	xaxis.tck = 0,
	xaxis.top.tck = 0,
	yaxis.tck = c(0.1,0),
	grid.row = TRUE,
	force.grid.row = TRUE,
	row.lines = 4.5,
	row.lwd = 3,
	row.col = 'black'
	);

cell.covariate <- create.heatmap(
	t(heatmap.data[rownames(heatmap.data) %in% cell.line,1:3]),
	same.as.matrix = TRUE,
	cluster.dimensions = 'none',
	print.colour.key = FALSE,
	colour.scheme = heatmap.colours,
	total.colours = length(heatmap.colours)+1,
	at = 0:length(heatmap.colours),
	axes.lwd = 1,
	yaxis.lab = NA,
	xaxis.lab = qw('SW1736 BHT101 ASH3 THJ-29T 8505C C643 CAL62 THJ-21T THJ-11T KMH2 THJ-16T UHTh7 UHTh74'),
	yaxis.cex = 1,
	xaxis.cex = 0.9,
	xaxis.tck = 0.1,
	xaxis.top.tck = 0,
	yaxis.tck = c(0.1,0),
	grid.row = TRUE,
	force.grid.row = TRUE,
	row.lines = 2.5,
	row.lwd = 3,
	row.col = 'black'
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
		)
	);

clinical.legend.trim[[5]] <- list();
clinical.legend.trim[[5]]$colours <- functionalColours;
clinical.legend.trim[[5]]$labels <- names(functionalColours);
clinical.legend.trim[[5]]$title <- 'Function';
names(clinical.legend.trim)[5] <- 'legend';

seqKey <- legend.grob(
	legends = seq.legend,
	size = 1.5,
	label.cex = 0.8,
	title.cex = 0.8,
	title.just = 'left',
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
	pch = c(19,18)[match(sapply(sampleCount[sampleCount$Sample %in% c(primary, non.atc),]$Sample, substr, 0, 8) %in% has.normal, c('TRUE','FALSE'))],
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
	xlab.cex = 1,
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
	xaxis.cex = 1,
	ylab.label = NULL,
	xlab.label = 'Proportion',
	xlab.cex = 1
	);

cellHeatmap <- create.heatmap(
	t(recurrenceData[rev(as.character(recurrenceCount$Gene)), colnames(recurrenceData) %in% cell.line]),
	cluster.dimensions = 'none',
	yaxis.lab = NA,
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
	xaxis.lab = NULL,
	xaxis.tck = 0,
	xaxis.top.tck = 0,
	yaxis.tck = c(0.1,0)
	);

### CREATE MULTIPLOT ##############################################################################
## combine the plots
# primary tumours
create.multipanelplot(
	plot.objects = list(mutCountPlot, basechangePlot, recurrenceHeatmap, signifPlot, atc.covariate),
	filename = generate.filename('AmpliSeq_filtered', 'RecurrentSNV_Profile_primary', 'tiff'),
	#resolution = 500,
	height = 7,
	width = 11,
	plot.objects.heights = c(1.2, 1.2, 4, 1.2),
	plot.objects.widths = c(9, 1),
	layout.width = 2,
	layout.height = 4,
	layout.skip = c(FALSE, TRUE, FALSE, TRUE, rep(FALSE,3), TRUE),
	x.spacing = 1.5,
	y.spacing = c(-3.5,-3.5,-4.5),
	ylab.cex = 1,
	use.legacy.settings = TRUE,
	ylab.axis.padding = -3.5,
	left.legend.padding = 0,
	right.legend.padding = 1,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	top.padding = -1,
	legend = list(
		inside = list(fun = seqKey, x = 0.92, y = 0.87),
		bottom = list(fun = clinicalKey)
		)
	);

create.multipanelplot(
	plot.objects = list(cellHeatmap, cellProp, cell.covariate),
	filename = generate.filename('AmpliSeq_filtered', 'RecurrentSNV_Profile_cellline', 'tiff'),
	#resolution = 500,
	height = 6,
	width = 5,
	plot.objects.heights = c(4, 1.2),
	plot.objects.widths = c(2, 1),
	layout.width = 2,
	layout.height = 2,
	layout.skip = c(FALSE, FALSE, FALSE, TRUE),
	x.spacing = 0,
	y.spacing = -3,
	use.legacy.settings = TRUE,
	left.legend.padding = 0,
	right.legend.padding = 3,
	bottom.legend.padding = 0,
	top.legend.padding = 0,
	top.padding = -1,
	);

### WRITE SESSION INFO TO FILE ####################################################################
save.session.profile(generate.filename('VariantHeatmap_AmpliSeqFiltered', 'SessionInfo', 'txt'));
