### plot_progression_vafs.R ########################################################################
# identify progression candidates from multi-region tumours and plot VAFs

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general);

# general parameters
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/Exomes/plots');

# general parameters
date <- Sys.Date();
cfg <- read.config.file();

### READ DATA #####################################################################################
# get sample annotation
phenodata <- read.delim('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/other/phenodata.txt');
phenodata <- phenodata[!is.na(phenodata$WEX.Site),];
phenodata <- phenodata[!phenodata$Individual %in% c('ANPT0074','ANPT0104','ANPT0105','ANPT0126','ANPT0150','ANPT0173'),];
phenodata <- phenodata[!phenodata$SampleName == 'ANPT0165PI',];

# move to data directory
setwd(cfg$readcount.directory);

these.files <- list.files(pattern = 'VAF_sample_by_position.tsv');
this.file <- rev(sort(these.files))[1];
positions <- read.delim(this.file);

# and get coverage file
these.files <- list.files(pattern = 'COVERAGE_sample_by_position.tsv');
this.file <- rev(sort(these.files))[1];
cov.data <- read.delim(this.file);

# file to determine if variant is likely germline
these.files <- list.files(pattern = 'Readcounts_InNormals.tsv');
this.file <- rev(sort(these.files))[1];
filter.data <- read.delim(this.file);

# get variant annotations
annotation <- read.delim('../../recSNV/unfiltered_pon/unfiltered_variants_by_patient.tsv');
annotation <- annotation[,1:9];

# remove non-functional variants
annotation <- annotation[-which(annotation$Location == 'intergenic'),];
annotation <- annotation[-which(annotation$Function == 'synonymous SNV'),];

### ORGANIZE DATA #################################################################################
# fix/clarify some sample names
phenodata$ID <- as.character(phenodata$SampleName);
phenodata[phenodata$SampleName == 'ANPT0021M',]$ID <- 'ANPT0021PM';
phenodata[phenodata$SampleName == 'ANPT0145PP',]$ID <- 'ANPT0145P';

# indicate sample groups
has.normal <- unique(substr(phenodata[phenodata$Type == 'Reference',]$SampleName, 0, 8));

cell.line <- intersect(
	unique(phenodata[phenodata$Type == 'CellLine',]$SampleName),
	names(positions)
	);

paired <- intersect(
	unique(phenodata[phenodata$Individual %in% has.normal,]$ID),
	names(positions)
	);
paired <- paired[-grep('R$|W$|I$|M$', paired)];

tumour.only <- setdiff(
	names(positions),
	c(paired, cell.line)
	);
tumour.only <- tumour.only[-grep('R$|W$|I$|M$', tumour.only)];

atc.samples <- c(paired,tumour.only);
ptc.samples <- phenodata[grep('W$', phenodata$ID),]$ID;

# find coverage metrics for normal only
filter.data$Count <- apply(
	filter.data[,grep('ANPT', names(filter.data))],
	1,
	function(i) { length(i[which(i >= 17)]) }
	);

pos.to.keep <- filter.data[filter.data$Count <= 4,];
positions <- merge(positions, pos.to.keep[,c('chr','pos','ref','alt')]);

### STATISTICS ####################################################################################
# are any variants present at a higher/lower frequency in ATC/PTC
multi.region <- sort(intersect(
	substr(atc.samples,0,8),
	substr(ptc.samples,0,8)
	));

## counts only ##
tmpData <- positions[,1:4];

for (subsample in multi.region) {

	atc <- atc.samples[grep(subsample, atc.samples)];
	ptc <- ptc.samples[grep(subsample, ptc.samples)];

	tmp <- positions[,c('chr','pos','ref','alt',atc,ptc)];
	tmp$diff <- tmp[,atc] - tmp[,ptc];

	tmpData[,subsample] <- 0;
	tmpData[which(tmp$diff >= 0.05),subsample] <- 1;
	tmpData[which(tmp$diff <= -0.05),subsample] <- -1;
	}

# identify recurrent variants
tmpData$ATC.Count <- apply(tmpData[,grep('ANPT',names(tmpData))],1,function(i) { length(i[which(i == 1)]) } );
tmpData$PTC.Count <- apply(tmpData[,grep('ANPT',names(tmpData))],1,function(i) { length(i[which(i == -1)]) } );

tmpData <- tmpData[tmpData$ATC.Count > 1 | tmpData$PTC.Count > 1,];

# now filter by coverage
for (subsample in multi.region) {

	atc <- atc.samples[grep(subsample, atc.samples)];
	ptc <- ptc.samples[grep(subsample, ptc.samples)];

	tmp <- tmpData[,c('chr','pos','ref','alt',subsample)];
	cov <- cov.data[,c('chr','pos','ref',atc,ptc)];
	cov <- cov[which(cov[,atc] >= 17 & cov[,ptc] >= 17),];
	
	tmp[,subsample] <- 0;
	for (row in 1:nrow(cov)) {

		chr <- as.character(cov[row,]$chr);
		pos <- cov[row,]$pos;
	
		if (nrow(tmp[tmp$chr == chr & tmp$pos == pos,]) == 0) { next; }
		else { tmp[tmp$chr == chr & tmp$pos == pos,subsample] <- 1; }
		}
	tmpData[,subsample] <- tmpData[,subsample] * tmp[,subsample];
	}

# again, identify recurrent variants
tmpData$ATC.Count <- apply(tmpData[,grep('ANPT',names(tmpData))],1,function(i) { length(i[which(i == 1)]) } );
tmpData$PTC.Count <- apply(tmpData[,grep('ANPT',names(tmpData))],1,function(i) { length(i[which(i == -1)]) } );

tmpData <- tmpData[tmpData$ATC.Count > 1 | tmpData$PTC.Count > 1,];

# format output
to.write <- merge(
	annotation[,1:6],
	tmpData[,-grep('ANPT',names(tmpData))]
	);

# get data for plotting
data.to.plot <- merge(
	merge(
		annotation[,c(4,5,6,8,9,1,3)],
		positions
		),
	tmpData[tmpData$ATC.Count >= 3 | tmpData$PTC.Count >= 2,-grep('ANPT',names(tmpData))]
	);

# set low depth positions to NA
for (row in 1:nrow(data.to.plot)) {

	chr <- as.character(data.to.plot[row,]$chr);
	pos <- data.to.plot[row,]$pos;

	cov <- cov.data[cov.data$chr == chr & cov.data$pos == pos,grep('ANPT',names(cov.data))];
	cov[cov < 17] <- NA;
	cov[!is.na(cov)] <- 1;

	for (col in names(data.to.plot)[grep('ANPT',names(data.to.plot))]) {
		data.to.plot[row,col] <- data.to.plot[row,col] * cov[,col];
		}
	}

# order by difference in VAF between tumour/normal
atc <- phenodata[phenodata$Individual %in% multi.region & phenodata$Differentiation != 'Reference' & !grepl('W$',phenodata$ID),]$ID;
ptc <- phenodata[phenodata$Individual %in% multi.region & phenodata$Differentiation != 'Reference' & grepl('W$',phenodata$ID),]$ID;

# find FC of ATC/normal
data.to.plot$Diff <- apply(
	data.to.plot,
	1,
	function(i) {
		get.foldchange(
			i,
			group1 = which(names(data.to.plot) %in% paste0(has.normal, 'R')),
			group2 = which(names(data.to.plot) %in% atc),
			logged = FALSE
			)
		}
	);

# split into chunks
data.to.plot <- data.to.plot[order(-data.to.plot$ATC.Count, -data.to.plot$Diff),];
atc.to.plot <- data.to.plot[which(data.to.plot$ATC.Count >= 3),];

data.to.plot <- data.to.plot[order(-data.to.plot$PTC.Count),];
ptc.to.plot <- data.to.plot[which(data.to.plot$PTC.Count >= 2),];

# make the plots!
atc.objects <- list();
ptc.objects <- list();
panel.heights <- c();
for (subsample in rev(sort(multi.region))) {

	atc <- atc.samples[grep(subsample, atc.samples)];
	ptc <- ptc.samples[grep(subsample, ptc.samples)];

	ylab.labels <- c(
		subsample,
		as.character(phenodata[phenodata$ID == ptc,]$Differentiation)
		);

	if (subsample %in% has.normal) {
		ref <- paste0(subsample,'R');
		panel.heights <- c(panel.heights, 3);
		ylab.labels <- c(ylab.labels, 'normal');
		} else {
		ref <- NULL;
		panel.heights <- c(panel.heights, 2);
		}

	# ATC variants
	atc.objects[[subsample]] <- create.dotmap(
		t(atc.to.plot[,c(atc,ptc,ref)]),
		bg.data = t(atc.to.plot[,c(atc,ptc,ref)]),
		colour.scheme = c('white','red','red'),
		at = seq(0,1,0.001),
		colour.centering.value = 0.5,
		yaxis.cex = 0.8,
		yaxis.lab = ylab.labels,
		colourkey = TRUE,
		colourkey.cex = 1,
		bg.alpha = 1,
		spot.colour.function = function(x) { rep('transparent', length(x)) },
		spot.size.function = function(x) { rep(0.001, length(x)) },
		NA.spot.size = 1.5
		);

	# PTC variants
	ptc.objects[[subsample]] <- create.dotmap(
		t(ptc.to.plot[,c(atc,ptc,ref)]),
		bg.data = t(ptc.to.plot[,c(atc,ptc,ref)]),
		colour.scheme = c('white','red','red'),
		at = seq(0,1,0.001),
		colour.centering.value = 0.5,
		yaxis.lab = ylab.labels,
		yaxis.cex = 0.8,
		colourkey = TRUE,
		colourkey.cex = 1,
		bg.alpha = 1,
		spot.colour.function = function(x) { rep('transparent', length(x)) },
		spot.size.function = function(x) { rep(0.001, length(x)) },
		NA.spot.size = 1.5
		);
	}

# make a covariate for the variant types
recsnv.scheme <- list(
	names = c('nonsynonymous SNV', 'stopgain/stoploss', 'splicing', 'UTR', 'other'),
	colours = c('darkseagreen4', 'darkturquoise', 'gold1', 'orchid4', 'white')
	);

covariate.atc <- rep('white', nrow(atc.to.plot));
covariate.atc[which(atc.to.plot$Function == 'nonsynonymous SNV')]		<- recsnv.scheme$colours[1];
covariate.atc[which(atc.to.plot$Function %in% c('stopgain','stoploss'))]	<- recsnv.scheme$colours[2];
covariate.atc[which(atc.to.plot$Function == 'splicing')]			<- recsnv.scheme$colours[3];
covariate.atc[which(atc.to.plot$Function == 'UTR')]				<- recsnv.scheme$colours[4];
#covariate.atc[which(atc.to.plot$Function == 'unknown')]				<- recsnv.scheme$colours[5];
#covariate.atc[which(atc.to.plot$Function == 'synonymous SNV')]			<- recsnv.scheme$colours[6];
covariate.atc[which(atc.to.plot$Function == 'other')]				<- recsnv.scheme$colours[5];

covariate.ptc <- rep('white', nrow(ptc.to.plot));
covariate.ptc[which(ptc.to.plot$Function == 'nonsynonymous SNV')]		<- recsnv.scheme$colours[1];
covariate.ptc[which(ptc.to.plot$Function %in% c('stopgain','stoploss'))]	<- recsnv.scheme$colours[2];
covariate.ptc[which(ptc.to.plot$Function == 'splicing')]			<- recsnv.scheme$colours[3];
covariate.ptc[which(ptc.to.plot$Function == 'UTR')]				<- recsnv.scheme$colours[4];
#covariate.ptc[which(ptc.to.plot$Function == 'unknown')]				<- recsnv.scheme$colours[5];
#covariate.ptc[which(ptc.to.plot$Function == 'synonymous SNV')]			<- recsnv.scheme$colours[6];
covariate.ptc[which(ptc.to.plot$Function == 'other')]				<- recsnv.scheme$colours[5];

atcCovariates <- list(
	rect = list(
		col = 'black',
		fill = covariate.atc,
		lwd = 1
		)
	);

ptcCovariates <- list(
	rect = list(
		col = 'black',
		fill = covariate.ptc,
		lwd = 1
		)
	);

# make the legend for the covariate
myLegends <- list(
	legend = list(
		labels = recsnv.scheme$names,
		colours = recsnv.scheme$colours,
		title = 'SNV Type'
		),
	legend = list(
		colours = c('white','red'),
		labels = c(0, expression(''>='0.5')),
		title = 'VAF',
		continuous = TRUE,
		height = 15,
		width = 1.5,
		pos.x = -0.2,
		tck = 0.5,
		tck.number = 1
		)
	);

# combine the plots
create.multiplot(
	atc.objects,
	filename = generate.filename('MultiRegion','variantAlleleFreq_ATC','tiff'),
	resolution = 800,
	height = 10,
	width = 10,
	plot.layout = c(1,length(atc.objects)),
	panel.heights = rev(panel.heights),
	y.relation = 'free',
	y.spacing = 0.5,
	yaxis.cex = 0.8,
	yaxis.tck = 0.5,
	right.padding = 25,
	left.padding = 5,
	top.padding = 3,
	xaxis.lab = atc.to.plot$Gene,
	xaxis.rot = 90,
	xaxis.cex = 0.8,
	axes.lwd = 2,
	print.new.legend = TRUE,
	legend = list(
		inside = list(fun = covariates.grob(atcCovariates, ord = 1:length(covariate.atc), side = 'top'), x = 0.5, y = 1.04),
		inside = list(fun = legend.grob(myLegends, size = 2, title.just = 'left', between.row = 4), x = 1.03, y = 1.01)
		)
	);

create.multiplot(
	ptc.objects,
	filename = generate.filename('MultiRegion','variantAlleleFreq_PTC','tiff'),
	resolution = 800,
	height = 10,
	width = 10,
	plot.layout = c(1,length(ptc.objects)),
	panel.heights = rev(panel.heights),
	y.relation = 'free',
	y.spacing = 0.5,
	yaxis.cex = 0.8,
	yaxis.tck = 0.5,
	right.padding = 25,
	left.padding = 5,
	top.padding = 3,
	xaxis.lab = ptc.to.plot$Gene,
	xaxis.rot = 90,
	xaxis.cex = 0.8,
	axes.lwd = 2,
	print.new.legend = TRUE,
	legend = list(
		inside = list(fun = covariates.grob(ptcCovariates, ord = 1:length(covariate.ptc), side = 'top'), x = 0.5, y = 1.04),
		inside = list(fun = legend.grob(myLegends, size = 2, title.just = 'left', between.row = 4), x = 1.03, y = 1.01)
		)
	);

write.table(
	to.write,
	file = generate.filename('MultiRegionVariants','topProgressionCandidates','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	quote = FALSE,
	sep = '\t'
	);

# not run!
{

### are proportions alt reads/total reads different between groups?
tmpData <- positions[,1:4];

for (subsample in multi.region) {

	tmpData[,subsample] <- NA;

	for (i in 1:nrow(tmpData)) {

		chrom <- as.character(positions[i,]$chr);
		pos   <- positions[i,]$pos;
		ref.allele <- as.character(positions[i,]$ref);
		alt.allele <- as.character(positions[i,]$alt);
	
		tmp <- readCounts[readCounts$chr == chrom & readCounts$pos == pos,];
		tmp <- tmp[grep(subsample, tmp$Sample),];

		if (nrow(tmp) < 2) { next; }
		if (all(tmp$depth == 0)) { next; }
		if (all(tmp[,alt.allele] == 0)) { next; }
		if (all(tmp[1,c('depth',alt.allele)] == 0)) { next; }
		if (all(tmp[2,c('depth',alt.allele)] == 0)) { next; }

		tmpData[i,subsample] <- prop.test(
			tmp[,alt.allele],
			tmp$depth
			)$p.value;
		}
	}

# identify recurrent variants
tmpData$Count <- apply(tmpData[,grep('ANPT',names(tmpData))],1,function(i) { length(i[which(i < 0.1)]) });

# filter variants
tmpData2 <- tmpData[tmpData$Count > 0,1:4];
tmpData2$effect.size <- NA;
tmpData2$p.value <- NA;

for (i in 1:nrow(tmpData2)) {

	chrom <- as.character(tmpData2[i,]$chr);
	pos   <- tmpData2[i,]$pos;
	ref.allele <- as.character(tmpData2[i,]$ref);
	alt.allele <- as.character(tmpData2[i,]$alt);

	vafs <- positions[positions$chr == chrom & positions$pos == pos & positions$alt == alt.allele,];

	tmp <- data.frame(
		Sample = c(atc,ptc),
		VAF = c(as.numeric(vafs[,atc]),as.numeric(vafs[,ptc]))
		);

	if (all(tmp[!is.na(tmp$VAF),]$VAF == 0)) { next; }
	
	tmpData2[i,c('p.value','effect.size')] <- get.ttest.p.and.foldchange(
		tmp$VAF,
		group1 = grepl('W$',tmp$Sample),
		group2 = !grepl('W$',tmp$Sample),
		paired = TRUE,
		logged = TRUE
		);
	}

tmpData2$padj.value <- p.adjust(tmpData2$p.value, 'fdr');

# apply filter
filtered.variants <- tmpData2[which(tmpData2$p.value < 0.1),];

### VISUALIZATIONS ################################################################################
# get data for plotting
data.to.plot <- merge(
	positions,
	filtered.variants
	);

data.to.plot$chr <- factor(data.to.plot$chr, levels = paste0('chr',c(1:22, 'X','Y')));
data.to.plot <- data.to.plot[order(data.to.plot$chr, data.to.plot$pos),];
data.to.plot$Label <- paste0(data.to.plot$chr,'_',data.to.plot$pos,'_',data.to.plot$ref,'>',data.to.plot$alt);

# make the plots!
plot.objects <- list();
for (subsample in rev(sort(multi.region))) {

	atc <- atc.samples[grep(subsample, atc.samples)];
	ptc <- ptc.samples[grep(subsample, ptc.samples)];

	# create the heatmap
	plot.objects[[subsample]] <- create.heatmap(
		data.to.plot[,c(ptc,atc)],
		cluster.dimensions = 'none',
		colour.scheme = c('white','red'),
		at = c(seq(0,0.5,0.05),1),
		yaxis.lab = NA,
		yaxis.cex = 0.8,
		colourkey.cex = 1
		);
	}

# make a covariate for significance
spot.colour.function <- function(x) {
	colours <- rep('white', length(x));
	colours[sign(x) == 1] <- 'darkorange';
	colours[sign(x) == -1] <- 'dodgerblue';
	colours[sign(x) == 0] <- 'transparent';
	return(colours);
	}

plot.objects[[length(plot.objects)]] <- create.dotmap(
	t(data.to.plot$effect.size),
	bg.data = t(-log10(data.to.plot$p.value)),
	colour.scheme = c('white','black'),
	spot.size.function = 1,
	spot.colour.function = spot.colour.function,
	bg.alpha = 1,
	at = seq(0,3,0.001),
	pch = 21,
	yaxis.lab = 't-test'
	);

create.multiplot(
	plot.objects,
	filename = generate.filename('MultiRegion','variantAlleleFreq_byStats','tiff'),
	resolution = 500,
	height = 10,
	width = 8,
	plot.layout = c(1,length(plot.objects)),
	y.relation = 'free',
	yaxis.cex = 0.8,
	yaxis.tck = 0.5,
	left.padding = 5,
	bottom.padding = 3,
	xaxis.lab = data.to.plot$Label,
	xaxis.rot = 90,
	xaxis.cex = 0.8
	);

### ALTERNATE METHOD ###############################################################################
tmpData <- positions[,1:4];

for (subsample in multi.region) {

	atc <- atc.samples[grep(subsample, atc.samples)];
	ptc <- ptc.samples[grep(subsample, ptc.samples)];

	tmp <- positions[,c(atc,ptc)];
	tmp$diff <- tmp[,atc] - tmp[,ptc];

	tmpData[,subsample] <- 0;
	tmpData[which(tmp$diff >= 0.01),subsample] <- 1;
	}

# identify recurrent variants
tmpData$Count <- apply(tmpData[,grep('ANPT',names(tmpData))],1,sum);

atc <- phenodata[phenodata$Individual %in% multi.region & phenodata$Differentiation != 'Reference' & !grepl('W$',phenodata$ID),]$ID;
ptc <- phenodata[phenodata$Individual %in% multi.region & phenodata$Differentiation != 'Reference' & grepl('W$',phenodata$ID),]$ID;

# filter variants
tmpData2 <- tmpData[tmpData$Count > 1,-grep('ANPT',names(tmpData))];
tmpData2$effect.size <- NA;
tmpData2$p.value <- NA;

for (i in 1:nrow(tmpData2)) {

	chrom <- as.character(tmpData2[i,]$chr);
	pos   <- tmpData2[i,]$pos;
	ref.allele <- as.character(tmpData2[i,]$ref);
	alt.allele <- as.character(tmpData2[i,]$alt);

	vafs <- positions[positions$chr == chrom & positions$pos == pos & positions$alt == alt.allele,];

	tmp <- data.frame(
		Sample = c(atc,ptc),
		VAF = c(as.numeric(vafs[,atc]),as.numeric(vafs[,ptc]))
		);

	if (all(tmp[!is.na(tmp$VAF),]$VAF == 0)) { next; }
	
	tmpData2[i,c('p.value','effect.size')] <- get.ttest.p.and.foldchange(
		tmp$VAF,
		group1 = grepl('W$',tmp$Sample),
		group2 = !grepl('W$',tmp$Sample),
		paired = TRUE,
		logged = TRUE
		);
	}

filtered.variants <- tmpData2[tmpData2$p.value < 0.1 & tmpData2$effect.size > 0.01,];

# get data for plotting
data.to.plot <- merge(
	positions,
	filtered.variants
	);

data.to.plot$chr <- factor(data.to.plot$chr, levels = paste0('chr',c(1:22, 'X','Y')));
data.to.plot <- data.to.plot[order(data.to.plot$chr, data.to.plot$pos),];
data.to.plot$Label <- paste0(data.to.plot$chr,'_',data.to.plot$pos,'_',data.to.plot$ref,'>',data.to.plot$alt);

# make the plots!
plot.objects <- list();
for (subsample in rev(sort(multi.region))) {

	atc <- atc.samples[grep(subsample, atc.samples)];
	ptc <- ptc.samples[grep(subsample, ptc.samples)];

	# create the heatmap
	plot.objects[[subsample]] <- create.heatmap(
		data.to.plot[,c(ptc,atc)],
		cluster.dimensions = 'none',
		colour.scheme = c('white','red'),
		at = c(seq(0,0.5,0.05),1),
		yaxis.lab = NA,
		yaxis.cex = 0.8,
		colourkey.cex = 1
		);
	}

# make a covariate for significance
spot.colour.function <- function(x) {
	colours <- rep('white', length(x));
	colours[sign(x) == 1] <- 'darkorange';
	colours[sign(x) == -1] <- 'dodgerblue';
	colours[sign(x) == 0] <- 'transparent';
	return(colours);
	}

plot.objects[[length(plot.objects)+1]] <- create.dotmap(
	t(data.to.plot$effect.size),
	bg.data = t(-log10(data.to.plot$p.value)),
	colour.scheme = c('white','black'),
	spot.size.function = 1,
	spot.colour.function = spot.colour.function,
	bg.alpha = 1,
	at = seq(0,2,0.001),
	pch = 21,
	yaxis.lab = 't-test'
	);

create.multiplot(
	plot.objects,
	filename = generate.filename('MultiRegion','variantAlleleFreq_byRecurrence','tiff'),
	resolution = 500,
	height = 10,
	width = 8,
	plot.layout = c(1,length(plot.objects)),
	y.relation = 'free',
	yaxis.cex = 0.8,
	yaxis.tck = 0.5,
	left.padding = 5,
	bottom.padding = 3,
	xaxis.lab = data.to.plot$Label,
	xaxis.rot = 90,
	xaxis.cex = 0.8
	);

### counts only ###
tmpData <- positions[,1:4];

for (subsample in multi.region) {

	atc <- atc.samples[grep(subsample, atc.samples)];
	ptc <- ptc.samples[grep(subsample, ptc.samples)];

	tmp <- positions[,c(atc,ptc)];
	tmp$diff <- tmp[,atc] - tmp[,ptc];

	tmpData[,subsample] <- 0;
	tmpData[which(tmp$diff >= 0.05),subsample] <- 1;
	}

# identify recurrent variants
tmpData$Count <- apply(tmpData[,grep('ANPT',names(tmpData))],1,sum);

annotation <- read.delim('../../recSNV/unfiltered_pon/unfiltered_variants_by_patient.tsv');

to.write <- merge(
	annotation[,1:6],
	tmpData[tmpData$Count > 1, -grep('ANPT',names(tmpData))]
	);

to.write <- to.write[order(-to.write$Count),];

# get data for plotting
data.to.plot <- merge(
	merge(
		annotation[,c(4,5,6,8,9)],
		positions
		),
	tmpData[tmpData$Count > 6,-grep('ANPT',names(tmpData))]
	);

data.to.plot <- data.to.plot[order(-data.to.plot$Count),];
data.to.plot$Label <- paste0(data.to.plot$chr,'_',data.to.plot$pos,'_',data.to.plot$ref,'>',data.to.plot$alt);

# make the plots!
plot.objects <- list();
for (subsample in rev(sort(multi.region))) {

	atc <- atc.samples[grep(subsample, atc.samples)];
	ptc <- ptc.samples[grep(subsample, ptc.samples)];

	# create the heatmap
	plot.objects[[subsample]] <- create.heatmap(
		data.to.plot[,c(ptc,atc)],
		cluster.dimensions = 'none',
		colour.scheme = c('white','red'),
		at = c(seq(0,0.5,0.05),1),
		yaxis.lab = NA,
		yaxis.cex = 0.8,
		colourkey.cex = 1
		);
	}

create.multiplot(
	plot.objects,
	filename = generate.filename('MultiRegion','variantAlleleFreq_byCount','tiff'),
	resolution = 500,
	height = 10,
	width = 8,
	plot.layout = c(1,length(plot.objects)),
	y.relation = 'free',
	yaxis.cex = 0.8,
	yaxis.tck = 0.5,
	left.padding = 5,
	bottom.padding = 3,
	xaxis.lab = data.to.plot$Gene,
	xaxis.rot = 45,
	xaxis.cex = 0.8
	);

write.table(
	to.write,
	file = generate.filename('MultiRegionVariants','topProgressionCandidates','tsv'),
	row.names = FALSE,
	col.names = TRUE,
	quote = FALSE,
	sep = '\t'
	);
}

### WRITE SESSION INFO TO FILE ####################################################################
save.session.profile(generate.filename('PlotVAFs', 'SessionInfo', 'txt'));
