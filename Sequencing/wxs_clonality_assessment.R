### clonality_assessment.R #########################################################################
# Perform 'clustering' of VAF data from IonTorrent, on progression samples, to assess clonality.

### LIBRARIES AND SETUP ###########################################################################
library(BoutrosLab.plotting.general);

# set working directory
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/Exomes/plots/');

# read in config file
cfg <- read.config.file();

### READ DATA  ####################################################################################
# get sample information
phenodata <- read.delim(cfg$phenodata, as.is = TRUE);
phenodata <- phenodata[!is.na(phenodata$WEX.Site),];
phenodata <- phenodata[!phenodata$Individual %in% c('ANPT0074','ANPT0104','ANPT0105','ANPT0126','ANPT0150','ANPT0173'),];
phenodata <- phenodata[!phenodata$SampleName == 'ANPT0165PI',];

# first, get variant calls
variants <- read.delim('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/Exomes/recSNV/ampliseq_filtered/filtered_variants_by_patient.tsv');

# move to data directory
setwd(cfg$readcount.directory);

these.files <- list.files(pattern = 'VAF_sample_by_position.tsv');
this.file <- rev(sort(these.files))[1];
vafData <- read.delim(this.file);

# and get coverage file
these.files <- list.files(pattern = 'COVERAGE_sample_by_position.tsv');
this.file <- rev(sort(these.files))[1];
cov.data <- read.delim(this.file);

### ORGANIZE DATA #################################################################################
# fix/clarify some sample names
phenodata$ID <- as.character(phenodata$SampleName);
phenodata[phenodata$SampleName == 'ANPT0021M',]$ID <- 'ANPT0021PM';
phenodata[phenodata$SampleName == 'ANPT0145PP',]$ID <- 'ANPT0145P';

# indicate sample groups
has.normal <- unique(substr(phenodata[phenodata$Type == 'Reference',]$SampleName, 0, 8));

paired <- intersect(
	unique(phenodata[phenodata$Individual %in% has.normal,]$ID),
	names(vafData)
	);
atc.samples <- paired[-grep('R$|W$|I$|M$', paired)];
ptc.samples <- c(paired[grep('W$|I$|M$', paired)], 'ANPT0174');

multi.region <- sort(intersect(
	substr(atc.samples,0,8),
	substr(ptc.samples,0,8)
	));

### VISUALIZATION ##################################################################################
# move to output directory
setwd('./clonality');

# make a function to determine colours
colour.function <- function(x,y) {
	scheme <- default.colours(3, palette = 'spiral.morning');
	diff <- x - y;
	colours <- rep(scheme[1],length(diff));
	colours[which(diff >= 0.2)]  <- scheme[2];
	colours[which(diff <= -0.2)] <- scheme[3];
        colours[which(y == 0)]       <- scheme[2];
	colours[which(x == 0)]       <- scheme[3];
	colours[which(x > 0.5 & y > 0.5)] <- scheme[1];
	colours[which(x < 0.2 & y < 0.2)] <- scheme[1];
	return(colours);
	}

plot.function <- function(formula, data, xlab.label = NULL, ylab.label = NULL, colours, text.data) {

	plot <- create.scatterplot(
		formula,
		data,
		alpha = 0.5,
		cex = 1,
		col = colours,
		xlab.label = xlab.label,
		ylab.label = ylab.label,
		xlab.cex = 1.5,
		ylab.cex = 1.5,
		xaxis.cex = 1,
		yaxis.cex = 1,
		xlimits = c(-0.05,0.8),
		xat = seq(0,1,0.2),
		ylimits = c(-0.05,0.8), 
		yat = seq(0,1,0.2),
		xaxis.tck = c(0.5,0),
		yaxis.tck = c(0.5,0),
		add.xyline = TRUE,
		xyline.col = 'grey50',
		use.legacy.settings = TRUE,
		strip.cex = 1.2,
		add.text = TRUE,
		text.cex = 0.8,
		#text.guess.labels = TRUE,
		text.labels = text.data$gene,
		text.x = text.data$x,
		text.y = text.data$y
		);

	return(plot);
	}

pointData <- list();

# extract plot data for each sample
for (i in multi.region) {

	if (i == 'ANPT0174') {
		# get variants
		myVars <- variants[,c('chr','pos','alt','Gene', 'ANPT0174','ANPT0175')];
		myVars$Count <- apply(myVars[,grep('ANPT',names(myVars))],1,sum);
		myVars <- myVars[which(myVars$Count > 0),-ncol(myVars)];
		tmp.ids <- c('ANPT0174P','ANPT0174R','ANPT0175P','ANPT0175R');
		atc <- 'ANPT0174P';
		ptc <- 'ANPT0175P';
		}

	else {
		# get variants
		myVars <- variants[,c('chr','pos','alt','Gene', names(variants)[grep(i,names(variants))])];
		myVars$Count <- apply(myVars[,grep('ANPT',names(myVars))],1,sum);
		myVars <- myVars[which(myVars$Count > 0),-ncol(myVars)];
		tmp.ids <- colnames(vafData)[grep(i, colnames(vafData))];
		atc <- tmp.ids[grepl('P$',tmp.ids)];
		ptc <- tmp.ids[!grepl('P$|R$',tmp.ids)];
		}

	norm <- tmp.ids[grepl('R$',tmp.ids)];

	if (length(ptc) > 1) { next; }

	# check coverage
	myCOV <- merge(
		myVars[,1:4],
		cov.data[,c('chr','pos',tmp.ids)]
		);
	myCOV <- myCOV[which(myCOV[,norm] >= 80),];
	myCOV <- myCOV[which(myCOV[,atc] >= 80 | myCOV[,ptc] >= 80),];

	# extract VAFs
	myVAFs <- merge(
		myCOV[,1:4],
		vafData[,c('chr','pos','alt',tmp.ids)]
		);

	if (nrow(myVAFs[!is.na(myVAFs[,atc]) & !is.na(myVAFs[,ptc]),]) == 0) { next; }

	pointData[[i]] <- myVAFs;

	}

genes.to.label <- list(
	'ANPT0015' = qw('MBOAT7 HIST1H1C HIST1H4J'),
	'ANPT0021' = qw('FLG GPRIN2 ASNSD1'),
	'ANPT0028' = qw('TP53 AIFM3 NFKBIB BRAF SCRAP'),
	'ANPT0137' = qw('IGSF22 IRF2BPL APC'),
	'ANPT0139' = qw('SLC35E2B MICAL3 FASTKD3'),
	'ANPT0147' = qw('SPSB2 ZNF626 KRTCAP3 ANO10 TLN1'),
	#'ANPT0148' = qw('SMARCA4'),
	#'ANPT0155' = qw('AIM1'),
	'ANPT0160' = qw('PIK3CB BICC1'),
	#'ANPT0177' = qw('SLC35E2B'),
	'ANPT0174' = qw('MUC6'),
	'ANPT0179' = qw('SLC35E2B')
	);

labelData <- list();

# extract label data for each sample
for (i in names(pointData)) {

	# determine points to label
	tmp <- pointData[[i]];

	if (i == 'ANPT0174') {
		tmp.ids <- c('ANPT0174P','ANPT0174R','ANPT0175P','ANPT0175R');
		atc <- 'ANPT0174P';
		ptc <- 'ANPT0175P';
		}

	else {
		tmp.ids <- colnames(tmp)[grep(i, colnames(tmp))];
		atc <- tmp.ids[grepl('P$',tmp.ids)];
		ptc <- tmp.ids[!grepl('P$|R$',tmp.ids)];
		}

	tmp <- tmp[which(tmp$Gene %in% genes.to.label[[i]]),c('Gene',atc,ptc)];

	if (nrow(tmp) > 0) {
		tmp2 <- tmp[tmp[,atc] > 0.05 | tmp[,ptc] > 0.05,];
		colnames(tmp2)[colnames(tmp2) == atc] <- 'y';
		colnames(tmp2)[colnames(tmp2) == ptc] <- 'x';
		colnames(tmp2)[1] <- 'gene';

		tmp2$x <- tmp2$x + 0.06;
		tmp2$y <- tmp2$y + 0.05;

		labelData[[i]] <- tmp2;
		}
	else { labelData[[i]] <- NULL }

	}

myPlots <- list();

# now, plot the data for each sample
for (i in names(pointData)) {
	
	# determine points to label
	tmp <- pointData[[i]];

	if (i == 'ANPT0174') {
		tmp.ids <- c('ANPT0174P','ANPT0174R','ANPT0175P','ANPT0175R');
		atc <- 'ANPT0174P';
		ptc <- 'ANPT0175P';
		}

	else {
		tmp.ids <- colnames(tmp)[grep(i, colnames(tmp))];
		atc <- tmp.ids[grepl('P$',tmp.ids)];
		ptc <- tmp.ids[!grepl('P$|R$',tmp.ids)];
		}

	ptc.label <- as.character(phenodata[phenodata$ID == ptc,]$Differentiation);

	# make the plots
	myPlots[[i]] <- plot.function(
		get(atc) ~ get(ptc) | i,
		data = pointData[[i]],
		colours = colour.function(y = tmp[,atc], x = tmp[,ptc]),
		xlab.label = paste0('VAF (', ptc.label, ')'),
		ylab.label = 'VAF (ATC)',
		text.data = labelData[[i]]
		);
	}

myPlots[['ANPT0015']] <- plot.function(
	ANPT0015PP ~ ANPT0015PW | 'ANPT0015',
	data = pointData[['ANPT0015']],
	colours = colour.function(y = pointData[['ANPT0015']]$ANPT0015PP, x = pointData[['ANPT0015']]$ANPT0015PW),
	xlab.label = 'VAF (PTC)',
	ylab.label = 'VAF (ATC)',
	text.data = data.frame(gene = labelData[['ANPT0015']]$gene, x = c(0.33, 0.45, 0.28), y = c(0.2, 0.11, 0.008))
	);

myPlots[['ANPT0021']] <- plot.function(
	ANPT0021P ~ ANPT0021PM | 'ANPT0021',
	data = pointData[['ANPT0021']],
	colours = colour.function(y = pointData[['ANPT0021']]$ANPT0021P, x = pointData[['ANPT0021']]$ANPT0021PM),
	xlab.label = 'VAF (Met)',
	ylab.label = 'VAF (ATC)',
	text.data = data.frame(gene = labelData[['ANPT0021']]$gene, x = c(0.22, 0.35, 0.1), y = c(0.09, 0.23, 0.31))
	);

myPlots[['ANPT0028']] <- plot.function(
	ANPT0028PP ~ ANPT0028PW | 'ANPT0028',
	data = pointData[['ANPT0028']],
	colours = colour.function(y = pointData[['ANPT0028']]$ANPT0028PP, x = pointData[['ANPT0028']]$ANPT0028PW),
	xlab.label = 'VAF (PTC)',
	ylab.label = 'VAF (ATC)',
	text.data = data.frame(gene = labelData[['ANPT0028']]$gene, x = c(0.07, 0.18, 0.445, 0.4), y = c(0.3, 0.05, 0.04, 0.305))
	);

myPlots[['ANPT0137']] <- plot.function(
	ANPT0137PP ~ ANPT0137PW | 'ANPT0137',
	data = pointData[['ANPT0137']],
	colours = colour.function(y = pointData[['ANPT0137']]$ANPT0137PP, x = pointData[['ANPT0137']]$ANPT0137PW),
	xlab.label = 'VAF (PTC)',
	ylab.label = 'VAF (ATC)',
	text.data = data.frame(gene = labelData[['ANPT0137']]$gene, x = c(0.42, 0.14, 0.19), y = c(0.18, 0.26, 0.11))
	);

myPlots[['ANPT0139']] <- plot.function(
	ANPT0139PP ~ ANPT0139PW | 'ANPT0139',
	data = pointData[['ANPT0139']],
	colours = colour.function(y = pointData[['ANPT0139']]$ANPT0139PP, x = pointData[['ANPT0139']]$ANPT0139PW),
	xlab.label = 'VAF (PTC)',
	ylab.label = 'VAF (ATC)',
	text.data = data.frame(gene = labelData[['ANPT0139']]$gene, x = c(0.6, 0.21, 0.34), y = c(0.36, 0.33, 0.17))
	);

myPlots[['ANPT0147']] <- plot.function(
	ANPT0147PP ~ ANPT0147PW | 'ANPT0147',
	data = pointData[['ANPT0147']],
	colours = colour.function(y = pointData[['ANPT0147']]$ANPT0147PP, x = pointData[['ANPT0147']]$ANPT0147PW),
	xlab.label = 'VAF (HTC)',
	ylab.label = 'VAF (ATC)',
	text.data = data.frame(gene = labelData[['ANPT0147']]$gene, x = c(0.25, 1, 0.58, 0.52, 0.37), y = c(0.65, 1, 0.68, 0.35, 0.5))
	);

myPlots[['ANPT0160']] <- plot.function(
	ANPT0160P ~ ANPT0160PI | 'ANPT0160',
	data = pointData[['ANPT0160']],
	colours = colour.function(y = pointData[['ANPT0160']]$ANPT0160P, x = pointData[['ANPT0160']]$ANPT0160PI),
	xlab.label = 'VAF (PDTC)',
	ylab.label = 'VAF (ATC)',
	text.data = data.frame(gene = labelData[['ANPT0160']]$gene, x = c(0.2, 0.52), y = c(0.62, 0.6))
	);

myPlots[['ANPT0179']] <- plot.function(
	ANPT0179PP ~ ANPT0179PW | 'ANPT0179',
	data = pointData[['ANPT0179']],
	colours = colour.function(y = pointData[['ANPT0179']]$ANPT0179PP, x = pointData[['ANPT0179']]$ANPT0179PW),
	xlab.label = 'VAF (PTC)',
	ylab.label = 'VAF (ATC)',
	text.data = data.frame(gene = labelData[['ANPT0179']]$gene, x = c(0.55), y = c(0.38))
	);

for (i in c(2,3,4,5,7,8,9,10)) {
	myPlots[[i]]$ylab$label <- '';
	}

create.multipanelplot(
	plot.objects = myPlots,
	filename = generate.filename('MultiRegion','comparisons','tiff'),
	#resolution = 500,
	height = 6,
	width = 11,
	layout.width = 5,
	layout.height = 2,
	plot.objects.heights = rep(1,2),
	plot.objects.widths = c(1.1,1,1,1,1),
	x.spacing = -0.4,
	y.spacing = -0.2,
	ylab.axis.padding = 1,
	xlab.axis.padding = 1,
	left.legend.padding = 0,
	right.legend.padding = 0,
	bottom.legend.padding = 1,
	top.legend.padding = 0,
	legend = list(
		bottom = list(fun = draw.key, args = list(
			key = list(
				points = list(col = default.colours(3, palette = 'spiral.morning'), cex = 1.5, alpha = 0.5, pch = 19),
				text = list(lab = c('clonal', 'subclonal (DTC)', 'subclonal (ATC)'), font = 2, cex = 1.5),
				columns = 3
				)
			), x = 0.75, y = 0.2
		))
	);

### SessionInfo ###################################################################################
save.session.profile(generate.filename('Clonality', 'Assessment','txt'));
