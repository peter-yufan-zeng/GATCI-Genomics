### plot_gaddygram.R ##############################################################################
# Using the PCAWG dataset, create a gaddygram plot demonstrating SNV rates across various tumour
# types; here we will include our ATC cohort for comparison.

### LIBRARIES AND SETUP ###########################################################################
library(BoutrosLab.plotting.general);

# set working directory
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/Exomes/plots/');

# read in config file
cfg <- read.config.file();

### READ DATA  ####################################################################################
# get sample information
phenodata <- read.delim(cfg$phenodata);

# get SNV counts (ATC)
mutationCounts <- read.delim(sub('filtered_pon','ampliseq_filtered',cfg$mutation.counts));
unfiltered.mutationCounts <- read.delim(sub('filtered_pon','ampliseq_unfiltered',cfg$mutation.counts));

has.normal <- unique(phenodata[phenodata$Type == 'Reference' & !is.na(phenodata$WEX.Site),]$Individual);

phenodata <- read.delim(cfg$platform.map);
phenodata <- phenodata[which(phenodata$WXS.Name %in% names(mutationCounts)),];
cell.lines <- as.character(phenodata[phenodata$Type == 'CellLine',]$WXS.Name);
pdtc.samples <- names(mutationCounts)[grep('W|I',names(mutationCounts))];
atc.samples <- setdiff(names(mutationCounts),c(pdtc.samples,cell.lines,'ANPT0021M'));

# get PCAWG SNV rates
pcawgData <- read.delim(
	'/.mounts/labs/boutroslab/private/CancerBiology/MutationDensity/Integration/Correlation/2017-10-16_SupplementaryTable1_MutationDensityMetrics.tsv'
	);
pcawgData <- pcawgData[,c('SampleID','TumourType','Sex','Age','MutsPerMb')];

# get callable bases (used for calculating snvs / mbp)
callable.bases <- read.delim('/u/sprokopec/ATC/Exomes/callable_bases/2017-09-19_CallableBases_unfiltered.tsv',
	header = FALSE
	);

# get estimated rates for GATCI PTC
ptc.rates <- read.delim('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/PTC/data/2018-09-19_EstimatedMutationRates_PTC.tsv');

wgs.rates <- read.delim('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/WholeGenome/plots/2019-09-11_ATCWGS_MutationRates.tsv', row.names = 1);

wgs.rates <- wgs.rates[c(
	paste0('ATCWGS.', paste0(33:41, 'A')),
	paste0('ATCWGS.', paste0(c(33:38,40), 'P')),
	'ATCWGS.39F','ATCWGS.41F'
	),];
wgs.rates <- wgs.rates[order(wgs.rates$Filtered.snv.rate),];
wgs.rates$TumourType <- 'co-occurring DTC';
wgs.rates[grepl('A$',rownames(wgs.rates)),]$TumourType <- 'ATC';

### ORGANIZE DATA #################################################################################
# remove pcawg samples without sex/age info
pcawgData <- droplevels(
	pcawgData[!is.na(pcawgData$Sex) & !is.na(pcawgData$Age),]
	);
pcawgData <- pcawgData[!is.na(pcawgData$MutsPerMb),];
pcawgData$Paired <- NA;

# format ATC data
sample.counts <- apply(mutationCounts[,colnames(mutationCounts) %in% atc.samples],2,sum);
sample.rates <- sample.counts;

for (i in 1:nrow(callable.bases)) {

	sample.name <- as.character(callable.bases[i,1]);
	count <- as.numeric(sample.counts[sample.name]);

	if (is.na(count)) { next; }
	sample.rates[sample.name] <- count/(callable.bases[i,2]/10**6);
	}

# add ATC to pcawgData
tmp <- data.frame(
	SampleID = names(sample.rates),
	TumourType = rep('ATC',length(sample.rates)),
	Sex = NA,
	Age = NA,
	MutsPerMb = as.numeric(sample.rates),
	Paired = 0
	);
tmp[which(tmp$SampleID %in% as.character(phenodata[phenodata$Individual %in% has.normal,]$WXS.Name)),]$Paired <- 1;

allData <- data.frame(rbind(pcawgData,tmp));

# format our W/PDTC data
sample.counts <- apply(mutationCounts[,colnames(mutationCounts) %in% pdtc.samples],2,sum);
sample.rates <- sample.counts;

for (i in 1:nrow(callable.bases)) {

	sample.name <- as.character(callable.bases[i,1]);
	count <- as.numeric(sample.counts[sample.name]);

	if (is.na(count)) { next; }
	sample.rates[sample.name] <- count/(callable.bases[i,2]/10**6);
	}

# add ATC to pcawgData
tmp <- data.frame(
	SampleID = names(sample.rates),
	TumourType = rep('co-occurring DTC',length(sample.rates)),
	Sex = NA,
	Age = NA,
	MutsPerMb = as.numeric(sample.rates),
	Paired = 0
	);
tmp[which(tmp$SampleID %in% as.character(phenodata[phenodata$Individual %in% has.normal,]$WXS.Name)),]$Paired <- 1;

allData <- data.frame(rbind(allData,tmp));

# add PTC (GATCI) data
tmp <- data.frame(
	SampleID = ptc.rates$Sample,
	TumourType = rep('PTC (GATCI)', nrow(ptc.rates)),
	Sex = NA,
	Age = NA,
	MutsPerMb = as.numeric(ptc.rates$Filtered),
	Paired = NA
	);

allData <- data.frame(rbind(allData,tmp));
	

# find median mutation load per tumour type (for ordering)
type.medians <- aggregate(
	allData$MutsPerMb,
	by = list(
		Type = allData$TumourType
		),
	median,
	na.rm = TRUE
	);
type.medians <- type.medians[order(-type.medians$x),];

# order tumour types
allData$TumourType <- factor(
	allData$TumourType,
	levels = type.medians$Type
	);

# organize by tumour type
data.to.plot <- allData[order(allData$TumourType, allData$MutsPerMb),];

# remove samples with no mutations
data.to.plot <- data.to.plot[data.to.plot$MutsPerMb > 0,];

data.to.plot$Order <- NA;
counter <- 0;

for (i in 1:nlevels(data.to.plot$TumourType)) {

	ttype <- levels(data.to.plot$TumourType)[i];

	n.samples <- nrow(data.to.plot[data.to.plot$TumourType == ttype,]);

	locs <- seq(counter, counter+1, length.out = n.samples + 2);

	data.to.plot[data.to.plot$TumourType == ttype,]$Order <- locs[2:(length(locs)-1)];

	counter <- counter + 1;

	}

# indicate yaxis labels
yaxis.labels <- c(
	#expression(paste(bold('10'^'-6'))),
	expression(paste(bold('10'^'-4'))),
	expression(paste(bold('10'^'-2'))),
	expression(paste(bold('10'^'0'))),
	expression(paste(bold('10'^'2'))),
	expression(paste(bold('10'^'4')))
	#expression(paste(bold('10'^'6'))),
	#expression(paste(bold('10'^'8'))),
	#expression(paste(bold('10'^'10')))
	);

# set up median lines
curves.list <- list();
for (i in 1:nrow(type.medians)) {
	curves.list[[i]] <- eval(
		parse(
			text = paste0(
				"function(x) { return(log10(type.medians[",
				i,
				",]$x)) }"
				)
			)
		);
	}

# set colours
data.to.plot$COL <- 'black';
data.to.plot[data.to.plot$TumourType == 'ATC',]$COL <- 'darkorchid4';
data.to.plot[data.to.plot$TumourType == 'co-occurring DTC',]$COL <- 'mediumorchid2';
data.to.plot[data.to.plot$TumourType == 'PTC (GATCI)',]$COL <- as.character(pcawg.colours(scheme = 'tumour.subtype', x = 'thy.adenoca'));
data.to.plot[data.to.plot$TumourType == 'Thy-AdenoCA',]$COL <- as.character(pcawg.colours(scheme = 'tumour.subtype', x = 'thy.adenoca'));

data.to.plot$PCH <- 19;
data.to.plot[which(data.to.plot$Paired == 0),]$PCH <- 4;

### CREATE THE PLOT ###############################################################################
# gaddy gram
create.scatterplot(
	log10(MutsPerMb) ~ Order,
	data.to.plot,
	filename = generate.filename('GaddyGram', 'MutsPerMbp', 'tiff'),
	#resolution = 800,
	height = 4,
	width = 10,
	cex = 0.5,
	alpha = 0.65,
	pch = data.to.plot$PCH,
	col = data.to.plot$COL,
	top.padding = 1,
	xlab.label = NULL,
	xlimits = c(0,34),
	xlab.cex = 1.5,
	xat = seq(0.5,35,1),
	xaxis.lab = sub('Thy-AdenoCA','PTC (PCAWG)',levels(data.to.plot$TumourType)),
	xaxis.rot = 90,
	xaxis.cex = 0.7,
	xaxis.tck = c(0.5,0),
	yaxis.tck = 0.5,
	ylab.label = 'SNVs/Mbp',
	ylab.cex = 1.25,
	ylimits = c(-4,4),
	yat = seq(-4,4,2),
	yaxis.cex = 1,
	yaxis.lab = yaxis.labels,
	add.curves = TRUE,
	curves.exprs = curves.list,
	curves.from = seq(0,38,1)+0.1,
	curves.to = seq(1,39,1)-0.1,
	curves.col = 'red',
	curves.lwd = 1,
	add.rectangle = TRUE,
	xleft.rectangle = seq(0,38,2),
	xright.rectangle = seq(1,39,2),
	ybottom.rectangle = -500,
	ytop.rectangle = 500,
	col.rectangle = 'grey90',
	alpha.rectangle = 0.5,
	use.legacy.settings = TRUE,
	add.points = TRUE,
	points.x = c(seq(26.1,26.9, length.out = 9),seq(18.1,18.9,length.out = 9)), #[match(wgs.rates$TumourType, c('ATC','co-occurring DTC'))],
	points.y = log10(wgs.rates$Filtered.snv.rate),
	points.col = 'lightseagreen',
	points.cex = 0.6,
	points.pch = 18
	);

# add cell lines and save
sample.counts <- apply(mutationCounts[,colnames(mutationCounts) %in% cell.lines],2,sum);
sample.rates <- sample.counts;

for (i in 1:nrow(callable.bases)) {

	sample.name <- as.character(callable.bases[i,1]);
	count <- as.numeric(sample.counts[sample.name]);

	if (is.na(count)) { next; }
	sample.rates[sample.name] <- count/(callable.bases[i,2]/10**6);
	}

# add ATC to pcawgData
tmp <- data.frame(
	SampleID = names(sample.rates),
	TumourType = rep('CellLine',length(sample.rates)),
	Sex = NA,
	Age = NA,
	MutsPerMb = as.numeric(sample.rates)
	);

allData <- data.frame(rbind(allData,tmp));
data.to.plot <- allData[order(allData$TumourType, allData$MutsPerMb),];
	
snv.rates <- data.to.plot[data.to.plot$TumourType %in% c('ATC','co-occurring DTC','CellLine'),c('SampleID','MutsPerMb')];

## for each sample, add the pre-recSNV filter mutation rates
snv.rates$UnfilteredRate <- NA;

for (group in c('atc.samples', 'pdtc.samples','cell.lines')) {

	smps <- as.character(get(group));
	sample.counts <- apply(unfiltered.mutationCounts[,colnames(unfiltered.mutationCounts) %in% smps],2,sum);
	sample.rates <- sample.counts;

	for (smp in smps) {

		count <- as.numeric(sample.counts[smp]);
		n.bases <- callable.bases[which(callable.bases$V1 == smp),]$V2;

		if (is.na(count) | is.na(n.bases)) { next; }
		snv.rates[which(snv.rates$SampleID == smp),]$UnfilteredRate <- count/(n.bases/10**6);

		}
	}

# write mutation rates to file
write.table(
	snv.rates,
	file = generate.filename('ATC', 'MutationRates', 'tsv'),
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

### SessionInfo ###################################################################################
save.session.profile(generate.filename('Plot', 'GaddyGram','txt'));
