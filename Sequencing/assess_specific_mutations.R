### assess_specific_mutations.R ###################################################################
# Examine the overlap between patients with mutations in certain genes 
# [ those that typically co-occur ]

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general);
library(VennDiagram);

# general parameters
setwd('/.mounts/labs/boutroslab/private/Collaborators/AnthonyNichols/AnaplasticThyroid/analysis_sdp/WholeGenome/plots/');

# general parameters
date <- Sys.Date();
cfg <- read.config.file();

### READ IN DATA FROM RECSNV ######################################################################
# get mutation calls
snv.data <- read.delim(cfg$recsnv.results);
snv.data <- snv.data[,!grepl('filter',names(snv.data))];
snv.data <- snv.data[-which(snv.data$Location == 'intergenic'),];
snv.data <- snv.data[-which(snv.data$Location == 'intronic'),];
snv.data <- snv.data[-which(snv.data$Location == 'upstream;downstream'),];

# trim down gene ids where necessary
snv.data$Symbol <- sapply(
        snv.data$Gene,
        function(i) {
                i <- as.character(i);
                i <- unlist(strsplit(i,'\\(|,|\\='))[1];
                return(i);
                }
        );

# extract ncRNAs
ncrnas <- unique(snv.data[grepl('ncRNA', snv.data$Location),]$Symbol);
ncrna.data <- snv.data[which(snv.data$Symbol %in% ncrnas),];
snv.data <- snv.data[-which(snv.data$Symbol %in% ncrnas),];

tmp <- aggregate(
        snv.data[,grepl('ANPT|ATCWGS',names(snv.data))],
        by = list(
                Symbol = snv.data$Symbol
                ),
        FUN = max
        );
rownames(tmp) <- tmp$Symbol;
tmp <- tmp[,-1];

for (gene in rownames(tmp)) {

        if (all(tmp[gene,] == 0)) { next; }

        for (smp in colnames(tmp)[which(tmp[gene,] == 1)]) {

                funct <- as.character(snv.data[which(snv.data$Symbol == gene & snv.data[,smp] == 1),]$Function);

                # if they're all NA, probably means upstream/downstream
                if (all(is.na(funct))) {
                        type <- as.character(snv.data[which(snv.data$Symbol == gene & snv.data[,smp] == 1),]$Location);
                        if (any(type %in% c('upstream','downstream'))) { tmp[gene,smp] <- 6; }
                        next;
                        }

                # remove any NAs now if necessary
                funct <- funct[!is.na(funct)];
                if (any(funct == 'nonsynonymous SNV')) { tmp[gene,smp] <- 1; }
                else if (any(grepl('stop', funct))) { tmp[gene,smp] <- 2; }
                else if (any(funct == 'splicing')) { tmp[gene,smp] <- 3; }
                else if (any(grepl('UTR', funct))) { tmp[gene,smp] <- 4; }
                else if (any(funct == 'synonymous SNV')) { tmp[gene,smp] <- 5; }
                else if (any(funct == 'unknown')) { tmp[gene,smp] <- 0; }
                else { stop('cant tell the variant type!!') }
                }
        }

gene.counts <- apply(tmp,1,function(i) { length(i[which(i > 0)]) - length(i[which(i == 5)]) });
snv.data <- tmp[which(gene.counts > 0),];

rm(gene.counts, tmp, funct, type, smp, gene);

### COMPARE SAMPLES ###############################################################################
# get rid of cell lines/non-atc
atc.data <- snv.data[,!grepl('P$|F$',colnames(snv.data))];
ptc.data <- snv.data[,grepl('P$|F$',colnames(snv.data))];

get.any.overlap <- function(i) {
	if (all(i == 0)) { return(0) } else { return(min(i[which(i > 0)])) }
	}

atc.data$ATCWGS.21 <- apply(atc.data[,grepl('21A',names(atc.data))], 1, get.any.overlap);
atc.data$ATCWGS.23 <- apply(atc.data[,grepl('23A',names(atc.data))], 1, get.any.overlap);
atc.data$ATCWGS.25 <- apply(atc.data[,grepl('25A',names(atc.data))], 1, get.any.overlap);
atc.data$ATCWGS.42 <- apply(atc.data[,grepl('42A',names(atc.data))], 1, get.any.overlap);
atc.data$ATCWGS.51 <- apply(atc.data[,grepl('51A',names(atc.data))], 1, get.any.overlap);

atc.data <- atc.data[,!grepl('21A',names(atc.data))];
atc.data <- atc.data[,!grepl('23A',names(atc.data))];
atc.data <- atc.data[,!grepl('25A',names(atc.data))];
atc.data <- atc.data[,!grepl('42A',names(atc.data))];
atc.data <- atc.data[,!grepl('51A',names(atc.data))];

# get sample names for each gene
sample.list <- list();

for (gene in c('EIF1AX','NRAS','HRAS','KRAS','BRAF','PIK3CA','DNAH17','TERT')) {
	counts <- apply(atc.data[rownames(atc.data) == gene,], 2, sum);
	these.samples <- names(counts)[which(counts > 0)];
	sample.list[[gene]] <- these.samples;
	}

# make plots and do calculations
venn.diagram(
	x = list(
		EIF1AX = sample.list[['EIF1AX']],
		RAS = unique(c(sample.list[['NRAS']], sample.list[['KRAS']], sample.list[['HRAS']]))
		),
	filename = generate.filename('eif1ax_ras', 'overlap', 'tiff'),
	resolution = 1000,
	height = 6,
	width = 6,
	units = 'in',
	fill = default.colours(2, is.venn = TRUE)[[1]],
	col = default.colours(2, is.venn = TRUE)[[2]],
	alpha = 0.5,
	cex = 4,
	fontface = 2,
	cat.cex = 3,
	cat.fontfamily = 'sans',
	fontfamily = 'sans',
	cat.fontface = 2,
	ext.text = FALSE,
	cat.dist = c(0.06,0.1),
	cat.pos = rev(c(210, 130)),
	margin = 0.15,
	hyper.test = TRUE,
	total.population = ncol(atc.data),
	sub.cex = 2,
	lower.tail = FALSE
	);

venn.diagram(
	x = list(
		BRAF = sample.list[['BRAF']],
		RAS = unique(c(sample.list[['NRAS']], sample.list[['KRAS']], sample.list[['HRAS']]))
		),
	filename = generate.filename('braf_ras', 'overlap', 'tiff'),
	resolution = 1000,
	height = 6,
	width = 6,
	units = 'in',
	fill = default.colours(2, is.venn = TRUE)[[1]],
	col = default.colours(2, is.venn = TRUE)[[2]],
	alpha = 0.5,
	cex = 4,
	fontface = 2,
	cat.cex = 3,
	cat.fontface = 2,
	cat.fontfamily = 'sans',
	fontfamily = 'sans',
	ext.text = FALSE,
	cat.dist = 0.05,
	cat.pos = 0,
	margin = 0.1,
	hyper.test = TRUE,
	total.population = ncol(atc.data),
	lower.tail = FALSE
	);

venn.diagram(
	x = list(
		EIF1AX = sample.list[['EIF1AX']],
		BRAF = sample.list[['BRAF']],
		NRAS = sample.list[['NRAS']]
		),
	filename = generate.filename('braf_ras_eif', 'overlap', 'tiff'),
	fill = default.colours(3, is.venn = TRUE)[[1]],
	col = default.colours(3, is.venn = TRUE)[[2]],
	alpha = 0.5,
	cex = 2,
	fontface = 2,
	cat.cex = 2,
	cat.fontface = 2,
	ext.text = FALSE,
	cat.dist = 0.05,
	margin = 0.1,
	hyper.test = TRUE,
	total.population = ncol(atc.data),
	lower.tail = FALSE
	);

venn.diagram(
	x = list(
		BRAF = sample.list[['BRAF']],
		PIK3CA = sample.list[['PIK3CA']]
		),
	filename = generate.filename('braf_pik3ca', 'overlap', 'tiff'),
	resolution = 1000,
	height = 6,
	width = 6,
	units = 'in',
	fill = default.colours(2, is.venn = TRUE)[[1]],
	col = default.colours(2, is.venn = TRUE)[[2]],
	alpha = 0.5,
	cex = 4,
	sub.cex = 2,
	fontface = 2,
	cat.cex = 3,
	cat.fontface = 2,
	cat.fontfamily = 'sans',
	fontfamily = 'sans',
	ext.text = FALSE,
	cat.dist = c(0.08,0.06),
	cat.pos = c(-40, 25),
	margin = 0.22,
	hyper.test = TRUE,
	total.population = ncol(atc.data),
	lower.tail = FALSE
	);

venn.diagram(
	x = list(
		DNAH17 = sample.list[['DNAH17']],
		PIK3CA = sample.list[['PIK3CA']]
		),
	filename = generate.filename('dnah17_pik3ca', 'overlap', 'tiff'),
	fill = default.colours(2, is.venn = TRUE)[[1]],
	col = default.colours(2, is.venn = TRUE)[[2]],
	alpha = 0.5,
	cex = 2,
	fontface = 2,
	cat.cex = 2,
	cat.fontface = 2,
	ext.text = FALSE,
	cat.dist = 0.05,
	cat.pos = c(-40, 40),
	margin = 0.1,
	hyper.test = TRUE,
	total.population = length(all.samples),
	lower.tail = FALSE
	);

### WRITE SESSION INFO TO FILE ####################################################################
save.session.profile(generate.filename('GeneOverlap_venndiagram', 'SessionInfo', 'txt'));
