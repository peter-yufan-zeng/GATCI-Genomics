### plot_progression_variants.R ####################################################################
# Purpose: Whole-exome sequencing was performed on ATC tumour and normal samples to identify genes 
#	containing recurrent somatic SNVs in tumour samples.
#	SNVs were called using MuTect (or somaticSniper) and resulting vcfs were run through recSNV
#	to generate basechange data, functional mutation data and a recurrent mutation list.
#	Recurrent SNVs were identified via multi-region samples (those that are exclusive to ATC).

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general);

# general parameters
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/Exomes/plots/');

# general parameters
date <- Sys.Date();
cfg <- read.config.file();

### READ IN DATA FROM RECSNV ######################################################################
setwd(cfg$readcount.directory);

# get progression variants
files <- list.files(pattern = 'MultiRegionVariants_topProgressionCandidates.tsv');
files <- rev(sort(files));
prog.vars <- read.delim(files[1]);
prog.genes <- c(
	as.character(unique(prog.vars[which(prog.vars$ATC.Count > 2),]$Gene)),
	as.character(unique(prog.vars[which(prog.vars$PTC.Count == 2),]$Gene))
	);

# read in recurrence matrix
recurrenceData <- read.delim(cfg$full.recurrence);
recurrenceData[recurrenceData == 5] <- NA;

# get sample annotation
phenodata <- read.delim('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/other/phenodata.txt');
phenodata <- phenodata[!is.na(phenodata$WEX.Site),];
phenodata <- phenodata[!phenodata$Individual %in% c('ANPT0074','ANPT0104','ANPT0105','ANPT0126','ANPT0150','ANPT0173'),];
phenodata <- phenodata[!phenodata$SampleName == 'ANPT0165PI',];

# clean up phenodata
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
phenodata[phenodata$SampleName == 'ANPT0177PP',]$ID <- 'ANPT0177P';
phenodata[phenodata$SampleName == 'ANPT0177PW',]$ID <- 'ANPT0177W';
phenodata[phenodata$SampleName == 'ANPT0179PP',]$ID <- 'ANPT0179P';
phenodata[phenodata$SampleName == 'ANPT0179PW',]$ID <- 'ANPT0179W';
phenodata[phenodata$SampleName == 'ANPT0182PP',]$ID <- 'ANPT0182PP';
phenodata[phenodata$SampleName == 'ANPT0182PW',]$ID <- 'ANPT0182PW';

### ORGANIZE DATA #################################################################################
# indicate sample groups
has.normal <- unique(substr(phenodata[phenodata$Type == 'Reference',]$SampleName, 0, 8));

cell.line <- intersect(
	unique(substr(phenodata[phenodata$Type == 'CellLine',]$SampleName, 0, 8)),
	names(recurrenceData)
	);

paired <- intersect(
	unique(phenodata[phenodata$Individual %in% has.normal,]$ID),
	names(recurrenceData)
	);
#paired <- paired[-grep('W$|I$|M$', paired)];

tumour.only <- setdiff(
	names(recurrenceData)[grep('ANPT',names(recurrenceData))],
	c(paired, cell.line)
	);
#tumour.only <- tumour.only[-grep('W$|I$|M$', tumour.only)];

non.atc <- unique(phenodata[grep('W$|I$|M$', phenodata$ID),]$ID);

# indicate samples with multi-region sampling
multi.region <- phenodata[phenodata$Individual %in% phenodata[phenodata$ID %in% non.atc,]$Individual & phenodata$Type != 'Reference',]$ID;
# remove samples with no poor-diff component
#multi.region <- multi.region[-grep('ANPT0127|ANPT0143', multi.region)];

# order samples
tmp.order <- aggregate(
	phenodata$Age,
	by = list(
		Sample = phenodata$Individual
		),
	FUN = length
	);
colnames(tmp.order)[2] <- 'N';
tmp.order <- tmp.order[order(-tmp.order$N),];

tmp.order$Group <- NA;
tmp.order[tmp.order$Sample %in% substr(multi.region,0,8),]$Group <- 1;

# remove reference samples
phenodata <- phenodata[!phenodata$Type == 'Reference',];

# continue describing sample order
tmp.order[tmp.order$Sample %in% paired & is.na(tmp.order$Group),]$Group <- 2;
tmp.order[tmp.order$Sample %in% tumour.only & is.na(tmp.order$Group),]$Group <- 2;
tmp.order[tmp.order$Sample %in% cell.line & is.na(tmp.order$Group),]$Group <- 3;

tmp.order <- tmp.order[!is.na(tmp.order$Group),];

phenodata.trim <- merge(
	phenodata[,c('ID','Individual','Differentiation','Source','Sex','Age','WEX.Site')],
	tmp.order[,c(1,3)],
	by.x = 'Individual',
	by.y = 'Sample'
	);

#rm(tmp.order);

# and further sort samples
phenodata.trim$Differentiation <- factor(
	phenodata.trim$Differentiation,
	levels = c('PTC','FTC','HTC','PDTC','ATC','Met')
	);

phenodata.trim$Tmp <- phenodata.trim$Group;
phenodata.trim <- phenodata.trim[order(
	phenodata.trim$Group,
	phenodata.trim$Differentiation,
	phenodata.trim$Individual
	),];
phenodata.trim$ID <- factor(
	phenodata.trim$ID,
	levels = phenodata.trim$ID
	);
phenodata.trim <- phenodata.trim[,-ncol(phenodata.trim)];

# filter data to desired genes
recurrenceData  <- recurrenceData[rownames(recurrenceData) %in% as.character(prog.genes),];

# set gene order
recurrenceCount <- apply(recurrenceData[,c(paired, tumour.only)], 1, function(i) { length(i[!is.na(i)]) });
recurrenceCount <- sort(recurrenceCount);

# filter/sort recurrence data
data.to.plot <- recurrenceData[names(recurrenceCount),as.character(phenodata.trim$ID)];

### VISUALIZATION #################################################################################
## make covariates
# get colour schemes
sex.colours		<- c('lightskyblue', 'lightpink');
age.colours		<- default.colours(2,'spiral.afternoon');
sampleGroup.colours	<- default.colours(3);
sampleType.colours	<- rev(default.colours(5,'spiral.sunrise'))[2:5];
seqSite.colours		<- default.colours(6,'spiral.morning');
sampleSource.colours	<- default.colours(5,'spiral.dusk')[c(2,3,4)];

covariate.data <- phenodata.trim[,c('ID','Differentiation','Sex','Age','Source','WEX.Site')];
rownames(covariate.data) <- covariate.data$ID;
covariate.data <- covariate.data[,-1];

# fill in type colours
covariate.data$Type.Colour	<- sampleType.colours[3];
covariate.data[which(covariate.data$Differentiation == 'PDTC'),]$Type.Colour			<- sampleType.colours[1];
covariate.data[which(covariate.data$Differentiation %in% c('PTC','HTC','FTC')),]$Type.Colour	<- sampleType.colours[2];
covariate.data[which(covariate.data$Differentiation == 'Met'),]$Type.Colour			<- sampleType.colours[4];

# fill in group colours
covariate.data$Group.Colour <- sampleGroup.colours[1];
covariate.data[tumour.only,]$Group.Colour <- sampleGroup.colours[2];
covariate.data[cell.line,]$Group.Colour   <- sampleGroup.colours[3];

# fill in sex colours
covariate.data$Sex.Colour <- 'darkgrey';
covariate.data[which(covariate.data$Sex == 'M'),]$Sex.Colour <- sex.colours[1];
covariate.data[which(covariate.data$Sex == 'F'),]$Sex.Colour <- sex.colours[2];

# fill in age colours
covariate.data$Age.Colour <- 'darkgrey';
covariate.data[which(covariate.data$Age == '0'),]$Age.Colour <- age.colours[1];
covariate.data[which(covariate.data$Age == '1'),]$Age.Colour <- age.colours[2];

# fill in sample source colours
covariate.data$Source.Colour <- 'darkgrey';
covariate.data[which(covariate.data$Source == 'Fresh'),]$Source.Colour  <- sampleSource.colours[1];
covariate.data[which(covariate.data$Source == 'Frozen'),]$Source.Colour <- sampleSource.colours[2];
covariate.data[which(covariate.data$Source == 'FFPE'),]$Source.Colour   <- sampleSource.colours[3];

# fill in sequencing site colours
covariate.data$Site.Colour <- seqSite.colours[1];
covariate.data[covariate.data$WEX.Site == 'Broad',]$Site.Colour    <- seqSite.colours[2];
covariate.data[covariate.data$WEX.Site == 'Columbia',]$Site.Colour <- seqSite.colours[3];
covariate.data[covariate.data$WEX.Site == 'Hopkins',]$Site.Colour  <- seqSite.colours[4];
covariate.data[covariate.data$WEX.Site == 'Sickkids',]$Site.Colour <- seqSite.colours[5];
covariate.data[covariate.data$WEX.Site == 'Yale',]$Site.Colour     <- seqSite.colours[6];

# make the covariate heatmaps
heatmap.data <- covariate.data[,6:10];
names(heatmap.data) <- c('MultiRegion','Type','Sex','Age','Source'); #,'Site');
heatmap.colours <- c();

i <- 0;
for (covariate in 1:ncol(heatmap.data)) {
	heatmap.data[,covariate] <- as.numeric(as.factor(covariate.data[,covariate+5])) + i;
	heatmap.colours <- c(heatmap.colours, levels(as.factor(covariate.data[,covariate+5])));
	i <- max(heatmap.data[,covariate]);
	}

# make a legend
clinical.legend <- list(
	legend = list(
		colours = rev(sampleType.colours),
		labels = c('Metastasis', 'ATC', 'PTC/HTC/FTC', 'PDTC'),
		title = 'Region'
		),
	legend = list(
		colours = sampleGroup.colours,
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
		),
	legend = list(
		colours = c(sampleSource.colours, 'darkgrey'),
		labels = c('Fresh','Frozen','FFPE','N/A'),
		title = 'Sample Source'
		)
#	legend = list(
#		colours = c(seqSite.colours),
#		labels = c('Baylor', 'Broad', 'Columbia', "John's Hopkins", 'Sickkids', 'Yale'),
#		title = 'Site'
#		)
	);

clinicalKey <- legend.grob(
	legends = clinical.legend,
	size = 1.5,
	label.cex = 0.8,
	title.cex = 0.8,
	title.just = 'left',
	layout = c(6,1) #c(3,2)
	);

# make the heatmaps
covariateHeatmaps <- list();
for (group in levels(as.factor(phenodata.trim$Group))) {

	these.samples <- phenodata.trim[phenodata.trim$Group == group,]$ID;

	covariateHeatmaps[[group]] <- create.heatmap(
		t(heatmap.data[rownames(heatmap.data) %in% these.samples,]),
		same.as.matrix = TRUE,
		cluster.dimensions = 'none',
		print.colour.key = FALSE,
		colour.scheme = heatmap.colours,
		total.colours = length(heatmap.colours)+1,
		at = 0:length(heatmap.colours)
		);
	}

# create colour schemes
# mutational consequence
functionalColours        <- c('darkseagreen4', 'darkturquoise', 'gold1', 'orchid4');
names(functionalColours) <- c('Nonsynonymous', 'Stopgain/Stoploss', 'Splicing', 'UTR');

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
	size = 1.5,
	label.cex = 0.8,
	title.cex = 0.8,
	title.just = 'left',
	layout = c(1,1)
	);

# create heatmap of sample x gene matrix (no clustering)
recurrenceHeatmaps <- list();
for (group in levels(as.factor(phenodata.trim$Group))) {

	these.samples <- phenodata.trim[phenodata.trim$Group == group,]$ID;

	recurrenceHeatmaps[[group]] <- create.heatmap(
		t(data.to.plot[,as.character(these.samples)]),
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
	}

### CREATE MULTIPLOT ##############################################################################
# make some lists!
myPlots <- list(
	covariateHeatmaps[[1]], covariateHeatmaps[[2]], covariateHeatmaps[[3]], #covariateHeatmaps[[4]],
	recurrenceHeatmaps[[1]], recurrenceHeatmaps[[2]], recurrenceHeatmaps[[3]] #, recurrenceHeatmaps[[4]]
	);

myAxes <- list(
	ylimits = list(
		NULL, NULL, NULL, 
		NULL, NULL, NULL
		),
	xlimits = list(
		NULL, NULL, NULL,
		NULL, NULL, NULL
		),
	yat = list(
		1:ncol(heatmap.data), 1:ncol(heatmap.data), 1:ncol(heatmap.data),
		1:nrow(data.to.plot), 1:nrow(data.to.plot), 1:nrow(data.to.plot)
		),
	xat = list(
		1:nrow(phenodata.trim[phenodata.trim$Group == 1,]),
		1:nrow(phenodata.trim[phenodata.trim$Group == 2,]),
		1:nrow(phenodata.trim[phenodata.trim$Group == 3,]),
		1:nrow(phenodata.trim[phenodata.trim$Group == 1,]),
		1:nrow(phenodata.trim[phenodata.trim$Group == 2,]),
		1:nrow(phenodata.trim[phenodata.trim$Group == 3,])
		),
	yaxis.lab = list(
		rev(colnames(heatmap.data)), rep('', ncol(heatmap.data)), rep('', ncol(heatmap.data)),
		rownames(data.to.plot), rep('', nrow(data.to.plot)), rep('', nrow(data.to.plot))
		),
	xaxis.lab = list(
		rep('', nrow(phenodata.trim[phenodata.trim$Group == 1,])),
		rep('', nrow(phenodata.trim[phenodata.trim$Group == 2,])),
		rep('', nrow(phenodata.trim[phenodata.trim$Group == 3,])),
		rep('', nrow(phenodata.trim[phenodata.trim$Group == 1,])),
		rep('', nrow(phenodata.trim[phenodata.trim$Group == 2,])),
		rep('', nrow(phenodata.trim[phenodata.trim$Group == 3,]))
		)
	);

# combine individual plots
create.multiplot(
	plot.objects = myPlots,
	filename = generate.filename(paste0(cfg$project.stem, '_', cfg$snv.caller), 'ProgressionVariants', 'tiff'),
	plot.layout = c(3,2),
	panel.widths = c(3, 8, 1),
	panel.heights = c(4, 0.7),
	resolution = cfg$resolution,
	height = 8,
	width = 11,
	bottom.padding = 10,
	x.spacing = -0.5,
	y.spacing = -0.5,
	xaxis.cex = 0.75,
	yaxis.cex = 0.6,
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
	legend = list(
		inside = list(fun = seqKey, x = 0.8, y = -0.03),
		inside = list(fun = clinicalKey, x = 0, y = -0.03)
		)
	);

### WRITE SESSION INFO TO FILE ####################################################################
save.session.profile(generate.filename('MultiRegion_recurrencePlot', 'SessionInfo', 'txt'));
