### plot_snv_summary_with_ampliseq.R ###############################################################
# Purpose: Whole-exome sequencing was performed on ATC tumour and normal samples to identify genes 
#	containing recurrent somatic SNVs in tumour samples.
#	SNVs were called using MuTect (or somaticSniper) and resulting vcfs were run through recSNV
#	to generate basechange data, functional mutation data and a recurrent mutation list.
#	AmpliSeq validation was used to filter FPs.
# Initial write: 2016-02-03
# Modified:	 2016-05-25
# Update: 	 2019-10-03

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

wgs.smps <- paste0('ATCWGS.',33:41);

### READ IN DATA FROM RECSNV ######################################################################
# get sample annotation
phenodata <- read.delim(cfg$platform.map);

has.normal <- c(
	as.character(phenodata[which(phenodata$Differentiation == 'Reference' & !is.na(phenodata$WXS.Site)),]$Individual),
	wgs.smps
	);

# read in basechange summary
basechangeSummary <- read.delim(gsub('filtered_pon', 'ampliseq_filtered', cfg$basechange.data));
tmp <- read.delim('~/ATC/WholeGenome/CallSNV/recSNV/most_recent_run/plots/basechange_proportion_data.txt');

wgs.smps.tmp <- c();
for (smp in wgs.smps) {
	wgs.smps.tmp <- c(wgs.smps.tmp, as.character(unique(tmp$Sample))[grepl(smp, unique(tmp$Sample))]);
	}
wgs.smps <- wgs.smps.tmp;
rm(wgs.smps.tmp);

tmp <- tmp[which(tmp$Sample %in% wgs.smps),];
basechangeSummary <- rbind(basechangeSummary, tmp);

# indicate samples to remove
to.remove <- c('ANPT0173');

# read in recurrence matrix
recurrenceData <- read.delim(gsub('filtered_pon', 'ampliseq_filtered', cfg$full.recurrence));
recurrenceData[recurrenceData == 5] <- NA;
recurrenceData <- recurrenceData[,!colnames(recurrenceData) %in% to.remove];

tmp <- read.delim('~/ATC/WholeGenome/CallSNV/recSNV/most_recent_run/filtered_variants_by_patient.tsv');
tmp <- tmp[!grepl('intergenic|intronic|upstream;downstream|ncRNA', tmp$Location),];
tmp <- tmp[which(!tmp$Function %in% c('synonymous SNV','unknown')),];

keep.fields <- names(tmp)[1:9];
for (smp in wgs.smps) {
	keep.fields <- c(keep.fields, names(tmp)[grep(sub('-','.',smp), names(tmp))]);
	}
tmp <- tmp[,keep.fields];

tmp$Count <- apply(tmp[,grepl('ATCWGS',names(tmp))],1,sum);
tmp <- tmp[which(tmp$Count > 0),];

for (smp in wgs.smps) {
	smp <- sub('-','.',smp);
	recurrenceData[,smp] <- NA;
	}

for (gene in rownames(recurrenceData)) {
	if (all(!grepl(gene, tmp$Gene))) { next; }
	x <- tmp[grepl(gene, tmp$Gene),];
	counts <- apply(x[,grepl('ATCWGS',names(x))],2,sum);

	for (smp in names(counts[which(counts > 0)])) {
		y <- x[which(x[,smp] == 1),];
		if (all(is.na(y$Function)) & any(y$Location %in% c('upstream','downstream'))) { recurrenceData[gene,smp] <- 6; }
		y <- y[!is.na(y$Function),];
		if (nrow(y) == 0) { next; }
		else if (any(y$Function == 'nonsynonymous SNV')) { recurrenceData[gene,smp] <- 1; }
		else if (any(grepl('stop', y$Function))) { recurrenceData[gene,smp] <- 2; }
		else if (any(y$Function == 'splicing')) { recurrenceData[gene,smp] <- 3; }
		else if (any(grepl('UTR', y$Function))) { recurrenceData[gene,smp] <- 4; }
		}
	}

for (smp in wgs.smps) {
	smp <- sub('-','.',smp);
	status <- recurrenceData['TERT',smp];
	if (is.na(status)) { tert.status[nrow(tert.status)+1,] <- c(smp, smp, 0); }
	else if (status == 6) { tert.status[nrow(tert.status)+1,] <- c(smp, smp, 'C228T'); }
	}

# manually fill this in (only necessary because I want to specify the mutations)
tert.status[which(tert.status$ID %in% c('ATCWGS.37P','ATCWGS.33P','ATCWGS.33A','ATCWGS.37A')),]$TERT <- 'C250T';
tert.status[nrow(tert.status)+1,] <- c('ANPT0014P','ANPT0014', '0'); # ATCWGS-17
tert.status[nrow(tert.status)+1,] <- c('ANPT0048P','ANPT0048', 'C228T'); # ATCWGS-18
tert.status[nrow(tert.status)+1,] <- c('ANPT0049P','ANPT0049', 'C228T'); # ATCWGS-19
tert.status[nrow(tert.status)+1,] <- c('ANPT0071P','ANPT0071', 'C228T'); # ATCWGS-26
tert.status[nrow(tert.status)+1,] <- c('ANPT0111P','ANPT0111', '0'); # ATCWGS-6
tert.status[nrow(tert.status)+1,] <- c('ANPT0118P','ANPT0118', '0'); # ATCWGS-13
tert.status[nrow(tert.status)+1,] <- c('ANPT0125P','ANPT0125', '0'); # ATCWGS-31

# get mutation counts per patient
mutationRates <- read.delim('2018-03-19_ATC_MutationRates.tsv');
mutationRates <- mutationRates[!mutationRates$Sample %in% to.remove,];

tmp <- read.delim('~/ATC/WholeGenome/plots/2019-09-11_ATCWGS_MutationRates.tsv');
tmp <- tmp[which(tmp$X %in% gsub('-','.',wgs.smps)),c('X','Filtered.snv.rate','MutationRate')];
colnames(tmp) <- colnames(mutationRates);
mutationRates <- rbind(mutationRates, tmp);

# get SeqSig results
seqsig <- read.delim('~/ATC/Exomes/recSNV/SeqSig/combined/ATC/2019-12-11_SeqSig_Driver_Gene_p_values.txt');

# get ampliseq results to add 7 additional smps
ampliseq <- read.delim('~/ATC/IonTorrent/readcounts/2018-03-06_VAF_sample_by_position.tsv');
ampli.cov <- read.delim('~/ATC/IonTorrent/readcounts/2018-03-06_COVERAGE_sample_by_position.tsv');

tmp <- read.delim('~/ATC/Exomes/recSNV/SeqSig/combined/2019-10-02_WXS_WGS_combined_filtVM.tsv');
tmp <- tmp[,1:9];

# add gene annotations
ampliseq <- merge(tmp, ampliseq);
ampli.cov <- merge(tmp, ampli.cov);

# apply a coverage filter
for (i in 1:nrow(ampliseq)) {
	chr <- as.character(ampliseq[i,]$chr);
	pos <- ampliseq[i,]$pos;
	smps <- ampli.cov[which(ampli.cov$chr == chr & ampli.cov$pos == pos),10:ncol(ampli.cov)];
	if (all(is.na(smps))) { ampliseq[i,names(smps)] <- NA; }
	else {
		smps <- smps[-which(smps >= 100)];
		ampliseq[i,names(ampliseq) %in% names(smps)] <- NA;
		}
	}

rm(ampli.cov, tmp);

# load covariate/clinical info
load('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/other/2018-09-11_ATC_clinical_covariates.RData');

### ORGANIZE DATA #################################################################################
# add clinical
phenodata <- merge(
	rbind(
		phenodata[,c(1,4,6)],
		data.frame(
			Individual = as.character(sapply(wgs.smps, substr, 0, 9)),
			CommonName = wgs.smps,
			WXS.Name = gsub('-','.',wgs.smps)
			)
		),
	atc.clinical$covariates,
	by.x = 'CommonName',
	by.y = 'row.names',
	all.x = TRUE
	);
phenodata[which(phenodata$CommonName %in% wgs.smps),]$Type <- sapply(
	phenodata[which(phenodata$CommonName %in% wgs.smps),]$CommonName,
	function(i) { 
		c('ATC','PTC','PTC')[match(substr(i, 10, 10), c('A','P','F'))];
		}
	);
phenodata[which(phenodata$CommonName %in% wgs.smps),]$Sex <- 'M';
phenodata[grepl('ATCWGS-33|ATCWGS-40', phenodata$Individual),]$Sex <- 'F';

# clean up phenodata
phenodata <- phenodata[phenodata$WXS.Name %in% names(recurrenceData),];

# indicate sample groups
cell.line <- intersect(colnames(recurrenceData), phenodata[phenodata$Type == 'CellLine',]$WXS.Name);
non.atc <- intersect(colnames(recurrenceData), phenodata[!phenodata$Type %in% c('CellLine','ATC'),]$WXS.Name);
primary <- intersect(colnames(recurrenceData), phenodata[phenodata$Type == 'ATC',]$WXS.Name);
primary <- setdiff(primary, 'ANPT0175');

# find recurrently altered genes
recurrenceCount <- data.frame(
	Gene  = rownames(recurrenceData),
	Count = apply(recurrenceData[,primary], 1, function(i) { length(i[!is.na(i)]) }),
	ATC.Prop = apply(recurrenceData[,primary], 1, function(i) { length(i[!is.na(i)])/length(i) }),
	DTC.Prop = apply(recurrenceData[,non.atc], 1, function(i) { length(i[!is.na(i)])/length(i) }),
	Diff.p = apply(
		recurrenceData,
		1,
		function(i) {
			atcs <- i[which(names(recurrenceData) %in% primary)];
			dtcs <- i[which(names(recurrenceData) %in% non.atc)];
			positive <- c(length(atcs[!is.na(atcs)]), length(dtcs[!is.na(dtcs)]));
			p.val <- prop.test(positive, c(length(atcs), length(dtcs)))$p.value;
			return(p.val);
			}
		),
	Cell.Prop = apply(recurrenceData[,cell.line], 1, function(i) { length(i[!is.na(i)])/length(i) })
	);
recurrenceCount <- recurrenceCount[order(-recurrenceCount$Count),];

primary <- c(primary, 'ANPT0175');

# identify recurrently altered genes
rec.genes <- as.character(recurrenceCount[recurrenceCount$Count >= 10,]$Gene);

recurrenceCount <- recurrenceCount[recurrenceCount$Gene %in% rec.genes,];
recurrenceCount <- recurrenceCount[order(-recurrenceCount$ATC.Prop),];

alt.genes <- as.character(recurrenceCount[grepl('^MUC|TTN', recurrenceCount$Gene),]$Gene);
keep.genes <- as.character(recurrenceCount[!grepl('^MUC|TTN', recurrenceCount$Gene),]$Gene);

sig.genes <- intersect(keep.genes, rownames(seqsig[which(seqsig$fdr.adjusted.p < 0.01),]));

# remove non-recurrent genes
recurrenceData  <- recurrenceData[sig.genes,];

# add in ampliseq smps
ampliseq.smps <- qw('ANPT0056P ANPT0057P ANPT0058P ANPT0059P ANPT0060P ANPT0061P ANPT0062P');
for (smp in ampliseq.smps) {
	recurrenceData[,smp] <- NA;
	}

for (gene in sig.genes) {
	if (nrow(ampliseq[which(ampliseq$Gene == gene),]) == 0) { next; }
	tmp <- apply(ampliseq[which(ampliseq$Gene == gene),ampliseq.smps],2,max,na.rm = TRUE);
	these.smps <- names(tmp)[which(tmp > 0.1)];
	tmp <- ampliseq[which(ampliseq$Gene == gene),c('Function',these.smps)];
	for (smp in these.smps) {
		funct <- as.character(tmp[order(-tmp[,smp]),]$Function[1]);
		recurrenceData[gene,smp] <- if (funct == 'nonsynonymous SNV') { 1 } else if (funct == 'splicing') { 3 } else { 2 };
		}
	}

ampliseq <- recurrenceData[,ampliseq.smps];
recurrenceData <- recurrenceData[,!names(recurrenceData) %in% ampliseq.smps];

# set gene order
recurrenceCount <- recurrenceCount[sig.genes,];
recurrenceCount$Gene <- factor(
	recurrenceCount$Gene,
	levels = rev(recurrenceCount$Gene)
	);

# sort samples by occurrence
tmpData <- t(ampliseq);
tmpData[!is.na(tmpData)] <- 1;
tmpData <- tmpData[do.call(order,transform(tmpData)),];
ampliseq.smps <- rownames(tmpData);

tmpData <- t(recurrenceData);
tmpData[!is.na(tmpData)] <- 1;
tmpData <- tmpData[do.call(order,transform(tmpData)),];

# get mutation counts
sampleCount <- rbind(
	data.frame(
		Sample = mutationRates[mutationRates$SampleID %in% phenodata$WXS.Name,]$SampleID,
		Count = log10(mutationRates[mutationRates$SampleID %in% phenodata$WXS.Name,]$MutsPerMb)
		),
	data.frame(
		Sample = ampliseq.smps,
		Count = NA
		)
	);

sampleCount$Sample <- factor(
	sampleCount$Sample,
	levels = c(rownames(tmpData),ampliseq.smps)
	);
sampleCount <- sampleCount[order(sampleCount$Sample),];

# ensure same sample ordering
basechangeSummary$Sample <- factor(
	gsub('-','.',basechangeSummary$Sample),
	levels = sampleCount$Sample
	);

for (smp in ampliseq.smps) {
	for (i in levels(basechangeSummary[,1])) {
		basechangeSummary[nrow(basechangeSummary)+1,] <- c(i, smp, 0, 0);
		}
	}
basechangeSummary$proportion <- as.numeric(basechangeSummary$proportion);
recurrenceData <- cbind(recurrenceData, ampliseq);
recurrenceData <- recurrenceData[,as.character(sampleCount$Sample)];

### VISUALIZATION #################################################################################
## make covariates
colour.data <- atc.clinical$covariate.data; 
colour.data <- colour.data[c(as.character(phenodata$CommonName),ampliseq.smps),];
rownames(colour.data) <- c(as.character(phenodata$WXS.Name), ampliseq.smps);

# don't show age/event covariates for cell lines
colour.data[cell.line,c('Age','Event')] <- 'white';

# trim covariates
#na.counts <- apply(colour.data,2,function(i) { length(i[grepl('grey80',i)])/length(i) } );
#keep.covariates <- names(na.counts[which(na.counts < 0.4)]);

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
tert.colours[which(colour.data$TERT == 'C228T')] <- 'darkorange';
tert.colours[which(colour.data$TERT == 'C250T')] <- 'orangered1';

colour.data$TERT <- tert.colours;

colour.data$Platform <- 'midnightblue'; # WXS
colour.data[grepl('ATCWGS',rownames(colour.data)),]$Platform <- 'mediumvioletred'; # WGS
colour.data[ampliseq.smps,]$Platform <- 'mediumspringgreen'; # Ampliseq

keep.covariates <- c('TERT','Platform','Type','Sex','Age','Event');
colour.data <- colour.data[,keep.covariates];

# fill in a few things
colour.data[is.na(colour.data$Event),c('Event','Age')] <- 'grey80';
colour.data[is.na(colour.data$Type),]$Type <- c('#65B4A2','#B1D39A','#B1D39A')[match(sapply(rownames(colour.data[is.na(colour.data$Type),]),substr,10,10),c('A','P','F'))];
colour.data[is.na(colour.data$Sex),]$Sex <- 'lightskyblue';
colour.data[grepl('ATCWGS.33|ATCWGS.40', rownames(colour.data)),]$Sex <- 'lightpink';

colour.data <- colour.data[colnames(recurrenceData),];

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

#clinical.legend[[2]]$colours <- clinical.legend[[2]]$colours[1:2];
#clinical.legend[[2]]$labels <- clinical.legend[[2]]$labels[1:2];

clinical.legend.trim <- list(
	legend = list(
		colours = c('darkorange','orangered1','white','grey80'),
		labels = c('C228T','C250T','WT','N/A'),
		title = 'TERT Status'
		),
	legend = list(
		colours = c('midnightblue','mediumvioletred','mediumspringgreen'),
		labels = c('WXS','WGS','IonTorrent'),
		title = 'Platform'
		),
	legend = clinical.legend[[1]],
	legend = clinical.legend[[2]],
	legend = clinical.legend[[3]],
	legend = clinical.legend[[8]]
	);

# make the heatmaps
atc.covariate <- create.heatmap(
	t(heatmap.data[rownames(heatmap.data) %in% c(primary, non.atc, ampliseq.smps),]),
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
	row.lines = 5.5,
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
	row.lines = 3.5,
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

clinical.legend.trim[[7]] <- list();
clinical.legend.trim[[7]]$colours <- functionalColours;
clinical.legend.trim[[7]]$labels <- names(functionalColours);
clinical.legend.trim[[7]]$title <- 'Function';
names(clinical.legend.trim)[7] <- 'legend';

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
	basechangeSummary[basechangeSummary$Sample %in% c(primary, non.atc, ampliseq.smps),],
	groups = basechangeSummary[basechangeSummary$Sample %in% c(primary, non.atc, ampliseq.smps),]$variable,
	stack = TRUE,
	col = basechangeColours,
	border.col = 'transparent',
	xaxis.lab = rep('', ncol(recurrenceData[,c(primary, non.atc, ampliseq.smps)])),
	yaxis.lab = seq(0, 1, 0.5),
	xat = 1:ncol(recurrenceData[,c(primary, non.atc, ampliseq.smps)]),
	yat = seq(0, 1, 0.5),
	ylab.label = '\nProportion',
	ylab.cex = 1,
	xlab.label = '',
	xaxis.cex = 0,
	yaxis.cex = 1,
	xaxis.tck = c(0.01,0),
	yaxis.tck = c(0.2,0)
	);

# create mutation count plot
mutCountPlot <- create.scatterplot(
	Count ~ Sample,
	sampleCount[sampleCount$Sample %in% c(primary, non.atc, ampliseq.smps),],
	cex = 0.5,
	pch = c(19,18)[match(
		sapply(
			sampleCount[sampleCount$Sample %in% c(primary, non.atc, ampliseq.smps),]$Sample,
			function(i) {
				if (grepl('ANPT',i)) { substr(i, 0, 8) } else { substr(i,0,9) }
				})
			%in% has.normal,
		c('TRUE','FALSE')
		)],
	xaxis.lab = rep('', ncol(recurrenceData[,c(primary, non.atc, ampliseq.smps)])),
	yaxis.lab = qw('0.01 0.1 1 10 100'),
	xat = 1:ncol(recurrenceData[,c(primary, non.atc, ampliseq.smps)]),
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
#seqsig[which(seqsig$nlog10.fdr.adjusted.p == 'Inf'),]$nlog10.fdr.adjusted.p <- 20;

signifPlot <- create.barplot(
	Gene ~ nlog10.fdr.adjusted.p,
	seqsig[sig.genes,],
	xlimits = c(0,10),
	xat = c(0,2,10), #seq(0,10,2),
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
	t(recurrenceData[rev(as.character(recurrenceCount$Gene)), colnames(recurrenceData) %in% c(primary, non.atc, ampliseq.smps)]),
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
	plot.objects.heights = c(1.2, 1.2, 4, 1.25),
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

write.table(
	recurrenceData,
	generate.filename('ATC_combined', 'RecurrentSNV_Profile', 'tsv'),
	row.names = TRUE,
	col.names = NA,
	sep = '\t'
	);

### WRITE SESSION INFO TO FILE ####################################################################
save.session.profile(generate.filename('VariantHeatmap_AmpliSeqFiltered', 'SessionInfo', 'txt'));
