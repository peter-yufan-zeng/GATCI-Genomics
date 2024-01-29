### plot_region_overlap.R #########################################################################
# Whole exome sequencing was performed at a number of facilities using various capture kits.
# Region overlap across kits was evaluated using bedtools multiinter.

### USAGE #########################################################################################
# Rscript plot_coverage.R WORKINGDIR

### LIBRARIES AND SETUP ###########################################################################
library(BoutrosLab.plotting.general);

# set working directory
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/Exomes/exome_capture/');

### READ DATA  ####################################################################################
# get sample information
regions <- read.delim('region_overlap.bed');
# save header
header <- names(regions);

# but rename for ease
colnames(regions) <- c('chrom','start','end','num','list', 'A','B','C','D','E','F');

# get region size
regions$size <- apply(
	regions[,c('start', 'end')],
	1,
	function(i)
		i[2] - i[1]
	);

# save kit names
kit.names <- data.frame(cbind(header, colnames(regions)[-ncol(regions)]));

### PAIRWISE COMPARISON ###########################################################################
# get pairwise overlap
# make a data frame for results
myData <- data.frame(matrix(nrow = 6, ncol = 6));
rownames(myData) <- c('A','B','C','D','E','F');
colnames(myData) <- c('A','B','C','D','E','F');

# for each interval file
for (i in colnames(myData)) {

	x <- regions[regions[,i] > 0,];

	for (j in rownames(myData)) {

		if (i == j) { next; }
		if (!is.na(myData[i,j])) { next; }

		myData[j,i] <- round(sum(x[x[,j] > 0,]$size)/1000000, 1);
		
		}
	}

# get text
myText <- myData;

# for each interval file
for (i in rownames(myText)) {
	for (j in colnames(myText)) {

		if (i == j) { myText[i,j] <- as.character(i); }
		else if (is.na(myData[i,j])) { myText[i,j] <- as.character(''); }
		else { myText[i,j] <- paste0(myData[i,j], 'Mbp'); }
		
		}
	}

# create heatmap
create.heatmap(
	myData,
	filename = generate.filename('pairwise_region', 'comparison', 'tiff'),
	resolution = 800,
	cluster.dimensions = 'none',
	print.colour.key = FALSE,
	cell.text = as.matrix(myText),
	text.col = 'black',
	text.cex = 1,
	text.fontface = 2,
	col.pos = which(myText != '1', arr.ind = TRUE)[,1],
	row.pos = which(myText != '1', arr.ind = TRUE)[,2],
	fill.colour = 'white',
	colour.scheme = c('lightblue1','dodgerblue3'),
	axes.lwd = 0	
	);

### SessionInfo ###################################################################################
save.session.profile(generate.filename('SessionInfo', 'ShowOverlap','txt'));
