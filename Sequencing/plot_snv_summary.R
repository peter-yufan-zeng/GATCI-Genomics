### plot_snv_summary.R ############################################################################
# Purpose: Whole-exome sequencing was performed on ATC tumour and normal samples to identify genes 
#	containing recurrent somatic SNVs in tumour samples.
#	SNVs were called using MuTect (or somaticSniper) and resulting vcfs were run through recSNV
#	to generate basechange data, functional mutation data and a recurrent mutation list.
#	Recurrent SNVs were identified via paired samples, with unpaired samples used for
#	validation.
# Initial write: 2016-02-03
# Modified:	 2016-05-25

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
basechangeSummary <- read.delim(cfg$basechange.data);

# read in functional change summary
functionalSummary <- read.delim(cfg$functional.data);

# reorder variable levels to match heatmap data
functionalSummary$variable <- factor(
	functionalSummary$variable,
	levels = c('Nonsynonymous', 'Stopgain-stoploss', 'Splicing', 'UTR')
	);
functionalSummary <- functionalSummary[!is.na(functionalSummary$variable),]

# read in recurrence matrix
recurrenceData <- read.delim(cfg$full.recurrence);
recurrenceData[recurrenceData == 5] <- NA;
recurrenceData <- recurrenceData[,!grepl('ANPT0173', colnames(recurrenceData))];

# get mutation counts per patient
mutationRates <- read.delim('2017-09-19_ATC_mutation_rates.tsv');
mutationRates <- mutationRates[!grepl('ANPT0173', mutationRates$Sample),];
 
### ORGANIZE DATA #################################################################################
# clean up phenodata
phenodata <- phenodata[!is.na(phenodata$WEX.Site),];
has.normal <- phenodata[phenodata$Type == 'Reference',]$Individual;
phenodata <- phenodata[!phenodata$Type == 'Reference',];

phenodata$Variants <- 'Tumour';
phenodata[phenodata$Individual %in% has.normal,]$Variants <- 'Somatic';
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

# fill in phenodata
inferred.sex <- read.delim('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/raw/OncoScan/inferred_sample_sex.txt');

for (i in 1:nrow(phenodata)) {
	if (!is.na(phenodata[i,]$Sex)) { next; }
	else {
		id <- as.character(phenodata[i,]$SampleName);
		if (!id %in% inferred.sex$SampleName) { next; }
		else {
			phenodata[i,]$Sex <- substr(inferred.sex[inferred.sex$SampleName == id,]$InferredSex,0,1);
			}
		}
	}

# indicate sample groups
tumour.only <- phenodata[phenodata$Variants == 'Tumour' & phenodata$Differentiation == 'ATC' & phenodata$Type != 'CellLine',]$ID;
cell.line <- phenodata[phenodata$Type == 'CellLine' & phenodata$Differentiation == 'ATC',]$ID;
non.atc <- phenodata[grep('W$|M$|I$', phenodata$ID),]$ID;
paired <- phenodata[phenodata$Variants == 'Somatic' & phenodata$Differentiation == 'ATC',]$ID;

# using only paired data, find recurrently altered genes
recurrenceCount <- data.frame(
	Gene  = rownames(recurrenceData),
	Count = apply(recurrenceData[,paired], 1, function(i) { length(i[!is.na(i)]) }),
	ATC.Prop = apply(recurrenceData[,c(paired, tumour.only)], 1, function(i) { length(i[!is.na(i)])/length(i) }),
	Cell.Prop = apply(recurrenceData[,cell.line], 1, function(i) { length(i[!is.na(i)])/length(i) })
	);
recurrenceCount <- recurrenceCount[order(-recurrenceCount$Count),];

# filter out genes primarily mutated in contaminated samples
rec.genes <- as.character(recurrenceCount[recurrenceCount$Count > 3 | recurrenceCount$ATC.Prop > 0.1,]$Gene);

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
# get colour schemes
sex.colours	  <- c('lightskyblue', 'lightpink');
age.colours	     <- default.colours(2,'spiral.afternoon');
sampleType.colours   <- default.colours(3);

covariate.data <- phenodata[!is.na(phenodata$ID),c('ID','Sex','Age')];
rownames(covariate.data) <- covariate.data$ID;
covariate.data <- covariate.data[,-1];

# fill in type colours
covariate.data$Type.Colour <- sampleType.colours[1];
covariate.data[tumour.only,]$Type.Colour <- sampleType.colours[2];
covariate.data[cell.line,]$Type.Colour   <- sampleType.colours[3];

# fill in sex colours
covariate.data$Sex.Colour <- 'darkgrey';
covariate.data[which(covariate.data$Sex == 'M'),]$Sex.Colour <- sex.colours[1];
covariate.data[which(covariate.data$Sex == 'F'),]$Sex.Colour <- sex.colours[2];

# fill in age colours
covariate.data$Age.Colour <- 'darkgrey';
covariate.data[which(covariate.data$Age == '0'),]$Age.Colour <- age.colours[1];
covariate.data[which(covariate.data$Age == '1'),]$Age.Colour <- age.colours[2];

# make the covariate heatmaps
heatmap.data <- covariate.data[,3:5];
names(heatmap.data) <- c('Type','Sex','Age');
heatmap.colours <- c();

i <- 0;
for (covariate in 1:ncol(heatmap.data)) {
	heatmap.data[,covariate] <- as.numeric(as.factor(covariate.data[,covariate+2])) + i;
	heatmap.colours <- c(heatmap.colours, levels(as.factor(covariate.data[,covariate+2])));
	i <- max(heatmap.data[,covariate]);
	}

# make a legend
clinical.legend <- list(
	legend = list(
		colours = sampleType.colours,
		labels = c('T/N pair', 'T only', 'Cell line'),
		title = 'Type'
		),
	legend = list(
		colours = c(sex.colours, 'darkgrey'),
		labels = c('Male','Female','N/A'),
		title = 'Sex'
		),
	legend = list(
		colours = c(age.colours, 'darkgrey'),
		labels = c(expression(''>='70'),expression(''<'70'),'N/A'),
		title = 'Age'
		)
	);

# make the heatmaps
# stop distinguising between paired/unpaired samples
atc.covariate <- create.heatmap(
	t(heatmap.data[rownames(heatmap.data) %in% c(paired,tumour.only),]),
	same.as.matrix = TRUE,
	cluster.dimensions = 'none',
	print.colour.key = FALSE,
	colour.scheme = heatmap.colours,
	total.colours = length(heatmap.colours)+1,
	at = 0:length(heatmap.colours)
	);

cell.covariate <- create.heatmap(
	t(heatmap.data[rownames(heatmap.data) %in% cell.line,]),
	same.as.matrix = TRUE,
	cluster.dimensions = 'none',
	print.colour.key = FALSE,
	colour.scheme = heatmap.colours,
	total.colours = length(heatmap.colours)+1,
	at = 0:length(heatmap.colours)
	);

# create colour schemes
# basechange
basechangeColours	<- force.colour.scheme(1:nlevels(basechangeSummary$variable), 'chromosome');
names(basechangeColours) <- levels(basechangeSummary$variable);

# mutational consequence
functionalColours	<- c('darkseagreen4', 'darkturquoise', 'gold1', 'orchid4');
names(functionalColours) <- c('Nonsynonymous', 'Stopgain/Stoploss', 'Splicing', 'UTR');

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
	layout = c(2,1)
	);

clinicalKey <- legend.grob(
	legends = clinical.legend,
	size = 1.5,
	label.cex = 0.8,
	title.cex = 0.8,
	title.just = 'left',
	layout = c(5,1) #c(3,2)
	);

## create individual plots:
# ATC tumours only
# create basechange barplot
basechangePlot <- create.barplot(
	proportion ~ Sample,
	basechangeSummary[basechangeSummary$Sample %in% c(paired,tumour.only),],
	groups = basechangeSummary[basechangeSummary$Sample %in% c(paired,tumour.only),]$variable,
	stack = TRUE,
	col = basechangeColours,
	border.col = 'transparent'
	);

# create functional consequence barplot
functionalPlot <- create.barplot(
	proportion ~ Sample,
	functionalSummary[functionalSummary$Sample %in% c(paired, tumour.only),],
	groups = functionalSummary[functionalSummary$Sample %in% c(paired, tumour.only),]$variable,
	stack = TRUE,
	col = functionalColours,
	border.col = 'transparent'
	);

# create mutation count plot
mutCountPlot <- create.scatterplot(
	Count ~ Sample,
	sampleCount[sampleCount$Sample %in% c(paired, tumour.only),],
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
	t(recurrenceData[rev(as.character(recurrenceCount$Gene)), colnames(recurrenceData) %in% c(paired, tumour.only)]),
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
	functionalPlot, basechangePlot, mutCountPlot
	);

myAxes <- list(
	ylimits = list(
		NULL, NULL,
		NULL, NULL,
		NULL, NULL,
		c(0,1), c(0,1), c(-1,3.5)
		),
	xlimits = list(
		NULL, NULL,
		NULL, c(0,1),
		NULL, c(0,1),
		NULL, NULL, NULL
		),
	yat = list(
		1:ncol(heatmap.data), 1:ncol(heatmap.data),
		1:nrow(recurrenceData), 1:nrow(recurrenceCount), 
		1:nrow(recurrenceData), 1:nrow(recurrenceCount), 
		seq(0, 1, 0.5),	seq(0, 1, 0.5),	seq(-1, 3, 1)	
		),
	xat = list(
		1:nrow(heatmap.data[c(paired, tumour.only),]), 1:nrow(heatmap.data[cell.line,]),
		1:ncol(recurrenceData[,c(paired, tumour.only)]), seq(0, 1, 0.5),
		1:ncol(recurrenceData[,cell.line]), seq(0, 1, 0.5),
		1:ncol(recurrenceData[,c(paired, tumour.only)]), 1:ncol(recurrenceData[,c(paired, tumour.only)]), 1:ncol(recurrenceData[,c(paired, tumour.only)])
		),
	yaxis.lab = list(
		rev(colnames(heatmap.data)), rep('', ncol(heatmap.data)),
		rev(rownames(recurrenceData)), rep('', nrow(recurrenceCount)), 
		rep('', nrow(recurrenceData)), rep('', nrow(recurrenceCount)),
		seq(0, 1, 0.5), seq(0, 1, 0.5), c(0.1,1,10,100,1000)
		),
	xaxis.lab = list(
		rep('', nrow(heatmap.data[c(paired, tumour.only),])), rep('', nrow(heatmap.data[cell.line,])),
		rep('', ncol(recurrenceData[,c(paired, tumour.only)])), seq(0, 1, 0.5),
		rep('', ncol(recurrenceData[,cell.line])), seq(0, 1, 0.5),
		rep('', ncol(recurrenceData[,c(paired, tumour.only)])), rep('', ncol(recurrenceData[,c(paired, tumour.only)])), rep('', ncol(recurrenceData[,c(paired, tumour.only)]))
		)
	);

# combine individual plots
create.multiplot(
	plot.objects = myPlots,
	filename = generate.filename(paste0(cfg$project.stem, '_', cfg$snv.caller), 'RecurrentSNV_Profile', 'tiff'),
	plot.layout = c(4,5),
	layout.skip = c(rep(c(FALSE, TRUE),2), rep(FALSE, 4), FALSE, rep(TRUE, 3), FALSE, rep(TRUE, 3), FALSE, rep(TRUE, 3)),
	panel.widths = c(10, 1, 1, 1), #c(7, 0.9, 4, 0.9, 1, 0.9),
	panel.heights = c(0.6, 1, 1, 6, 0.6),
	resolution = cfg$resolution,
	height = 9,
	width = 11,
	x.spacing = -0.5,
	y.spacing = -0.5,
	xaxis.cex = 0.75,
	yaxis.cex = 0.6,
	ylab.label = '\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t Proportion \t\tSNVs/Mbp',
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
	bottom.padding = 0,
	legend = list(
		inside = list(fun = seqKey, x = 0.76, y = 0.9),
		inside = list(fun = clinicalKey, x = 0.76, y = 1),
		inside = list(fun = draw.key, args = list(key = list(text = list(lab = parse(text = "bold('Proportion')"), cex = 0.75))), x = 0.75, y = 0.715),
		inside = list(fun = draw.key, args = list(key = list(text = list(lab = parse(text = "bold('Proportion')"), cex = 0.75))), x = 0.915, y = 0.715)
		)
	);

### WRITE SESSION INFO TO FILE ####################################################################
save.session.profile(generate.filename(paste0(cfg$project.stem, '_', cfg$snv.caller), 'SessionInfo', 'txt'));
