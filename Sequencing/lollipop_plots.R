### lollipop_plots.R ##############################################################################
# Generate lollipop plots demonstrating mutation data for key recurrently altered genes.

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

# read in recurrence matrix
#variantData <- read.delim('../recSNV/ampliseq_filtered/filtered_variants_by_patient.tsv');
variantData <- read.delim('~/ATC/Exomes/recSNV/SeqSig/combined/2019-10-02_WXS_WGS_combined_filtVM.tsv');
variantData <- variantData[,!grepl('ANPT0173', colnames(variantData))];
variantData <- variantData[,!grepl('ANPT0175', colnames(variantData))];
variantData <- variantData[-which(variantData$Location == 'intronic'),];

# load covariate/clinical info
load('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/other/2018-09-11_ATC_clinical_covariates.RData');

# get domain info
key.domains <- read.delim('2018-02-07_ProteinData_lolliplots.tsv', comment.char = '#');

### ORGANIZE DATA #################################################################################
# add clinical
phenodata <- merge(
	phenodata[,c(1,4,6)],
	atc.clinical$covariates,
	by.x = 'CommonName',
	by.y = 'row.names'
	);

# clean up phenodata
phenodata <- phenodata[phenodata$WXS.Name %in% names(variantData),];

# indicate sample groups
cell.line <- intersect(colnames(variantData), phenodata[phenodata$Type == 'CellLine',]$WXS.Name);
non.atc <- c(
	intersect(colnames(variantData), phenodata[grep('W$|M$|I$', phenodata$WXS.Name),]$WXS.Name),
	colnames(variantData)[grepl('ATCWGS', colnames(variantData)) & !grepl('A$', colnames(variantData))]
	);
primary <- c(
	intersect(colnames(variantData), phenodata[phenodata$Type == 'ATC',]$WXS.Name),
	colnames(variantData)[grepl('ATCWGS', colnames(variantData)) & grepl('A$', colnames(variantData))]
	);

# keep only recurrently altered genes
keep.genes <- as.character(key.domains$Gene);
variantData <- variantData[which(variantData$Gene %in% keep.genes),];

# get protein change
variantData$Change <- unlist(sapply(
	variantData$Result,
	function(i) {
		parts <- unlist(strsplit(unlist(strsplit(as.character(i),','))[1],':'));
		aa.change <- sub('p.','',parts[grep('^p\\.', parts)]);
		pos <- substr(aa.change,2,nchar(aa.change)-1);
		if (length(pos) > 0) { return(pos); } else { return(NA); }
		}
	));

variantData$AA <- unlist(sapply(
	variantData$Result,
	function(i) {
		parts <- unlist(strsplit(unlist(strsplit(as.character(i),','))[1],':'));
		aa.change <- sub('p.','',parts[grep('^p\\.', parts)]);
		#pos <- substr(aa.change,nchar(aa.change),nchar(aa.change));
		if (length(aa.change) > 0) { return(aa.change); } else { return(NA); }
		}
	));

# manually set the one splice variant
variantData[variantData$Function == 'splicing',]$Change <- 197;
variantData[variantData$Function == 'splicing',]$AA <- '';

### VISUALIZATION #################################################################################
# create colour schemes
functionalColours	 <- c('darkseagreen4', 'darkturquoise', 'darkgrey', 'gold1'); #, 'orchid4');
names(functionalColours) <- c('Nonsynonymous', 'Stopgain', 'Synonymous', 'Splicing'); #, 'UTR');

# create legend/colour key
seq.legend <- list(
	legend = list(
		colours = functionalColours,
		labels = names(functionalColours),
		title = 'Function'
		)
	);

seqKey <- legend.grob(
	legends = seq.legend,
	size = 2,
	label.cex = 1,
	title.cex = 1.2,
	title.just = 'left',
	use.legacy.settings = TRUE
	);

setwd('./hotspots');

# create a plot for each gene
for (gene in keep.genes) {

	# format data
	tmpData <- variantData[grep(gene,variantData$Gene),!grepl('filter',names(variantData))];
	tmpData <- tmpData[,c('Function','pos','Change','AA',names(tmpData)[grep('ANPT|ATCWGS',names(tmpData))])];

	gene.data <- reshape(
		tmpData,
		direction = 'long',
		varying = list(5:ncol(tmpData)),
		timevar = 'Sample',
		times = names(tmpData)[5:ncol(tmpData)],
		v.names = 'Type'
		);		
	rownames(gene.data) <- 1:nrow(gene.data);
	gene.data <- gene.data[which(gene.data$Type != 0),-ncol(gene.data)];

	atc.data <- gene.data[which(gene.data$Sample %in% primary),];

	data.to.plot <- aggregate(
		atc.data$Type,
		by = list(
			Function = atc.data$Function,
			POS = atc.data$Change,
			AA = atc.data$AA
			),
		FUN = function(i) { length(i[!is.na(i)]) } #/length(primary) }
		);

	data.to.plot$Function <- factor(
		data.to.plot$Function,
		levels = c('synonymous SNV','nonsynonymous SNV','stopgain','splicing'),
		labels = c('darkgrey','darkseagreen4','darkturquoise','gold1')
		);
	data.to.plot <- data.to.plot[order(data.to.plot$Function, data.to.plot$POS),];

	labelData <- data.to.plot[which(data.to.plot$x >= 3),];
	labelData <- labelData[!grepl('darkgrey',labelData$Function),];
		
	xat <- if (key.domains[which(key.domains$Gene == gene),]$Protein.Length > 10000) {
		c(seq(0,key.domains[which(key.domains$Gene == gene),]$Protein.Length,10000), key.domains[which(key.domains$Gene == gene),]$Protein.Length)
		} else if (key.domains[which(key.domains$Gene == gene),]$Protein.Length > 5000) {
		c(seq(0,key.domains[which(key.domains$Gene == gene),]$Protein.Length,2000), key.domains[which(key.domains$Gene == gene),]$Protein.Length) 
		} else if (key.domains[which(key.domains$Gene == gene),]$Protein.Length > 1000) {
		c(seq(0,key.domains[which(key.domains$Gene == gene),]$Protein.Length,500), key.domains[which(key.domains$Gene == gene),]$Protein.Length)
		} else if (key.domains[which(key.domains$Gene == gene),]$Protein.Length > 500) {
		c(seq(0,key.domains[which(key.domains$Gene == gene),]$Protein.Length,200), key.domains[which(key.domains$Gene == gene),]$Protein.Length)
		} else if (key.domains[which(key.domains$Gene == gene),]$Protein.Length > 200) {
		c(seq(0,key.domains[which(key.domains$Gene == gene),]$Protein.Length,100), key.domains[which(key.domains$Gene == gene),]$Protein.Length)
		} else if (key.domains[which(key.domains$Gene == gene),]$Protein.Length > 100) {
		c(seq(0,key.domains[which(key.domains$Gene == gene),]$Protein.Length,50), key.domains[which(key.domains$Gene == gene),]$Protein.Length)
		} else {
		c(seq(0,key.domains[which(key.domains$Gene == gene),]$Protein.Length,10), key.domains[which(key.domains$Gene == gene),]$Protein.Length)
		}

	# create individual plots
	create.lollipopplot(
		x ~ POS,
		data.to.plot,
		filename = generate.filename('ATCVariantRecurrence', gene, 'tiff'),
		resolution = 1600,
		width = 8,
		height = 4,
		col = as.character(data.to.plot$Function), #qw('darkseagreen4 darkturquoise darkgrey'),
		lwd = 0.5,
		cex = 1.5,
		pch = 21,
		xlimits = c(0,key.domains[which(key.domains$Gene == gene),]$Protein.Length),
		ylimits = c(0,max(data.to.plot$x)+1),
		xat = xat,
		yat = if (max(data.to.plot$x) > 15) { seq(0,ceiling(max(data.to.plot$x)/10)*10,5)
			} else if (max(data.to.plot$x) > 10) { seq(0,15,5)
			} else if (max(data.to.plot$x) > 5) { seq(0,10,2)
			} else { seq(0,max(data.to.plot$x) + 1, 1) },
		ylab.label = 'Number of ATC tumours',
		ylab.cex = 1.5,
		xlab.label = paste0(gene, ' Position'),
		xlab.cex = 1.5,
		xaxis.cex = 1,
		yaxis.cex = 1,
		xaxis.tck = c(0.5,0),
		yaxis.tck = c(0.5,0),
		regions.start = as.numeric(unlist(strsplit(as.character(key.domains[which(key.domains$Gene == gene),]$Domain.Start),';'))),
		regions.stop = as.numeric(unlist(strsplit(as.character(key.domains[which(key.domains$Gene == gene),]$Domain.End),';'))),
		regions.labels = unlist(strsplit(as.character(key.domains[which(key.domains$Gene == gene),]$Domain.Name),';')),
		regions.color = default.colours(6),
		regions.cex = 0.7,
		legend = if (gene == 'TP53') { list(inside = list(fun = seqKey, x = 0.05, y = 0.98)) } else { NULL },
		ylab.axis.padding = 2,
		add.text = (nrow(labelData) > 0),
		text.labels = labelData$AA,
		text.y = labelData$x,
		text.x = as.numeric(labelData$POS) - 200,
		use.legacy.settings = TRUE
		);
	}

### WRITE SESSION INFO TO FILE ####################################################################
save.session.profile(generate.filename('LollipopPlots', 'SessionInfo', 'txt'));
