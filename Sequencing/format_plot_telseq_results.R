### format_plot_telseq_results.R ##################################################################
# Collate telSeq results, identify outliers (if multiple lanes), calculate mean/median and plot
# results.

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general);
library(outliers);

date <- Sys.Date();

setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/WholeGenome/Telomere/');

### READ DATA #####################################################################################
# find files
telseq.output <- list.files(pattern = 'telseq_tumour.txt|telseq_normal.txt', recursive = TRUE);

# initiate list to hold data
my.data <- list();

# read in each file
for (i in 1:length(telseq.output)) {
	file <- telseq.output[i];
	my.data[[i]] <- read.delim(file, as.is = TRUE);
	}

# combine results
combined <- do.call(rbind, my.data);

# remove repeat entries (23N)
combined <- unique(combined);

### OUTLIERS ######################################################################################
# write function to check for outliers
grubbs.flag <- function(x) {
	outliers <- NULL;
	test <- x;
	if (length(test) < 4) { return(NA); }
	else {
		grubbs.result <- grubbs.test(test);
		pv <- grubbs.result$p.value;
		while(pv < 0.01) {
			outliers <- c(outliers,as.numeric(strsplit(grubbs.result$alternative," ")[[1]][3]));
			test <- x[!x %in% outliers];
			grubbs.result <- grubbs.test(test);
			pv <- grubbs.result$p.value;
			}
		if (any(x %in% outliers)) { return(TRUE) } else { return(FALSE) }
		}
	}

outliers <- aggregate(LENGTH_ESTIMATE ~ Sample, combined, grubbs.flag);

for (i in 1:nrow(outliers)) {
	if (outliers[i,2]) {
		id <- outliers[i,1];
		combined <- combined[!grepl(id, combined$Readgroup),];
		}
	}

### FORMAT DATA ###################################################################################
# calculate mean/median
to.write <- aggregate(LENGTH_ESTIMATE ~ Sample, combined, mean);
colnames(to.write) <- c('id','mean');
to.write$median <- aggregate(LENGTH_ESTIMATE ~ Sample, combined, median)[,2];
to.write$sd <- aggregate(LENGTH_ESTIMATE ~ Sample, combined, sd)[,2];

# save results
write.table(
	to.write[order(as.character(to.write$id)),],
	file = 'cohort_telseq_estimates.txt',
	row.names = FALSE,
	col.names = TRUE,
	sep = '\t'
	);

# format for plotting
to.plot <- to.write;

to.plot$Type <- 'A';
to.plot[grepl('N',to.plot$id),]$Type <- 'N';
to.plot[grepl('P$|F$',to.plot$id),]$Type <- 'P';

to.plot$Sample <- sapply(
	to.plot$id,
	function(i) {
		i <- gsub('ATCWGS-','',as.character(i));
		tmp <- substr(i,0,nchar(i)-1);
		tmp <- sub('A','',tmp);
		return(tmp);
		}
	);

to.plot$Type <- factor(to.plot$Type, levels = c('A','P','N'));
to.plot <- to.plot[order(to.plot$Type, -to.plot$mean),];
to.plot$Sample <- factor(to.plot$Sample, levels = unique(to.plot$Sample));
to.plot <- to.plot[order(to.plot$Sample, to.plot$Type),];
to.plot$id <- factor(to.plot$id, levels = to.plot$id);

### PLOT DATA #####################################################################################
# make/move to plot directory
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/WholeGenome/plots');

label.pos <- c(get.line.breaks(to.plot$Sample),nrow(to.plot)+0.5) - 
	(c(get.line.breaks(to.plot$Sample),nrow(to.plot)+0.5) - c(0.5, get.line.breaks(to.plot$Sample)))/2;

# create a barplot
create.barplot(
	mean ~ id,
	to.plot,
	filename = generate.filename('ATCWGS','TelSeq_results','tiff'),
	width = 9,
	height = 4.5,
	use.legacy.settings = TRUE,
	xaxis.tck = c(1,0),
	yaxis.tck = c(1,0),
	col = c('#65B4A2','#B1D39A','#FFE1EE')[match(to.plot$Type, c('A','P','N'))],
	xat = label.pos,
	xaxis.lab = paste0(c('','\n'), as.character(unique(to.plot$Sample))),
	#xaxis.rot = 90,
	xaxis.cex = 0.9,
	xlab.label = NULL,
	ylab.label = 'TelSeq Estimate',
	ylab.axis.padding = 2,
	ylimits = c(0,12),
	yat = seq(0,12,2),
	y.error.up = to.plot$sd,
	y.error.down = to.plot$sd,
	error.whisker.width = 0.01,
	y.error.bar.col = 'black', #c('#65B4A2','#B1D39A','#FFE1EE')[match(to.plot$Type, c('N','P','A'))],
	add.rectangle = TRUE,
	xleft.rectangle = c(0, get.line.breaks(to.plot$Sample)),
	xright.rectangle = c(get.line.breaks(to.plot$Sample),nrow(to.plot)+0.5),
	ytop.rectangle = 15,
	ybottom.rectangle = -5,
	col.rectangle = c('grey80','white'),
	alpha.rectangle = 0.8,
	key = list(
		x = 0.72, y = 0.9,
		rect = list(
			col = c('#65B4A2','#B1D39A','#FFE1EE'),
			size = 2
			),
		text = list(
			lab = c('ATC','PTC/FTC','Normal'),
			cex = 1,
			font = 2
			),
		background = 'white'
		)
	);

### SessionInfo ###################################################################################
save.session.profile(generate.filename('PlotTelSeq', 'SessionInfo', 'txt'));
