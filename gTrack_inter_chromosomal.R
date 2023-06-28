## This script generates gTrack plots with interactions occuring between two chromosomes with genome track 
## and any bigwig track.

#You can have a seventh column indicating strength of loop...

#install.packages("devtools")

## allows dependencies that throw warnings to install
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) 

## install with  the necessary gtrack packages with devtools
#devtools::install_github("mskilab/gTrack")
#devtools::install_github("mskilab/gUtils")

# ======= gTrack=========

# Load gTrack
library(gTrack)
library(gUtils)

rm(list=ls())
options(scipen = 100)

# Define and set working directory
# EDIT HERE #
data_path<-"~/Dropbox/HiChIP_KM/gTrack" 
setwd(data_path)

# Read in HiChIP/interaction/loop data

# Chromosome 1 interactions file
# EDIT HERE #
intra_1 <- read.table("./data/MCF7_allValidPairs.chr20.5kb.counts.combined.chr20.bedpe")

# Chromosome 2 interactions file
# EDIT HERE #
intra_2 <- read.table("./data/MCF7_allValidPairs.chr17.5kb.counts.combined.chr17.bedpe")

# Chromosome 1 and 2 interactions file
# EDIT HERE #
inter <- read.table("./data/MCF7_allValidPairs.chr17.chr20.5kb.counts.combined.bedpe")

# Merge intra-chromosomal interactions
intra <- rbind(intra_1,intra_2)

intra <- intra[!(intra$V2==intra$V5),] # remove the ones that left and right anchors are the same
intra <- intra[!(intra$V3==intra$V5),] # remove the ones that left and right anchors are immediately adjacent

inter <- inter[!(inter$V2==inter$V5),] # remove the ones that left and right anchors are the same
inter <- inter[!(inter$V3==inter$V5),] # remove the ones that left and right anchors are immediately adjacent


# Define regions of interest from two different chromosomes
# EDIT HERE #
#Chr8
chr14 = "chr8"
chr14.start <- 128000000 # Start chromosome location for viewing
chr14.end <- 130000000 # End chromosome location for viewing
#chr14
chr7 = "chr17"
chr7.start <- 37000000 # Start chromosome location for viewing
chr7.end <- 39000000 # End chromosome location for viewing

# Define bigwig files
# EDIT HERE #
H3K27ac_file <- "./data/MCF7_E2_H3K27ac.hg19.bw"
#DNase_file <- "./data/CUTLL1_BRD4.hg19.bw"
#CN_file <- "./data/ZR75_CN.bw"

# Get intra-chromosomal loops on chr7 within defined range
intra1.7 <- intra[intra$V1==chr7&intra$V2>=chr7.start&intra$V5<=chr7.end,]
# Get intra-chromosomal loops on chr14 within defined range
intra1.14 <- intra[intra$V1==chr14&intra$V2>=chr14.start&intra$V5<=chr14.end ,]
# Combine loops from chr7 and chr14
intra1 <- rbind(intra1.7,intra1.14)

# Get inter-chromosomal loops between chr7 and chr14 within the defined range
inter1.7_14 <- inter[inter$V1==chr7&inter$V2>=chr7.start&inter$V3<=chr7.end&inter$V4==chr14&inter$V5>=chr14.start&inter$V6<=chr14.end,]
inter1.14_7 <- inter[inter$V1==chr14&inter$V2>=chr14.start&inter$V3<=chr14.end&inter$V4==chr7&inter$V5>=chr7.start&inter$V6<=chr7.end,]
# Combine loops from chr7 and chr14
inter1 <- rbind(inter1.7_14,inter1.14_7)

# Create combined list of unique inter- and intra- chromosomal loops within our region of interest
ins <- rbind(intra1,inter1)
ins <- unique(ins)

# Create chr postition id for each loop anchor and associate the PET value to ins1 variable
ins$fea1 <- paste(ins$V1,ins$V2,ins$V3,sep="_")
ins$fea2 <- paste(ins$V4,ins$V5,ins$V6,sep="_")
ins1 <- ins[,c(8,9,7)]
#ins1$V8 <- scale(ins1$V8)

fea.all<- unique(c(ins$fea1,ins$fea2))
fea.all <- fea.all[order(fea.all,decreasing = T)]

# Create a matrix of zeros 
mt <- matrix(0, nrow = length(fea.all), ncol = length(fea.all))
#mt <- matrix(-0.3, nrow = length(fea.all), ncol = length(fea.all))

# Assign chromatin anchor position ID to the rows and columns of matrix
rownames(mt) <- fea.all
colnames(mt) <- fea.all

# Assign PET values for loops within the matrix
for (i in 1:nrow(ins1)){
  NN1 <- which(rownames(mt)==ins1$fea1[i])
  NN2 <- which(colnames(mt)==ins1$fea2[i])
  mt[NN1,NN2] <- log10(ins1$V7[i])
}

# Create gene range variable for plotting
ss.chr <- unname(sapply(fea.all,function(x){strsplit(x,"_")[[1]][1]}))
ss.start <- as.numeric(unname(sapply(fea.all,function(x){strsplit(x,"_")[[1]][2]})))
ss.end <- as.numeric(unname(sapply(fea.all,function(x){strsplit(x,"_")[[1]][3]})))
gr <- GRanges(seqnames = Rle(ss.chr), 
              ranges = IRanges(start = ss.start,
                               end = ss.end))
chr7.pos <- as.numeric(unname(sapply(fea.all[grepl(chr7,fea.all)],function(x){strsplit(x,"_")[[1]][2:3]})))
chr14.pos <- as.numeric(unname(sapply(fea.all[grepl(chr14,fea.all)],function(x){strsplit(x,"_")[[1]][2:3]})))
chr7.range <- paste0(chr7,":",chr7.start,"-",max(chr7.pos))
chr14.range <- paste0(chr14,":",chr14.start,"-",chr14.end)
wins <- GRanges(c(chr7.range,chr14.range))

# Create data for plot using gTrack
# Here you can change the color scheme, size, and label of the heatmap
dat.gt <- gTrack(gr, mdata = mt, height = 20, colormaps = c("white","orange","red","black"), name = 'HiChIP')
par(mar=c(0,0,0,0))

# Create gene annotation data
# Here you can change the size of the bigwig tracks and even define the genes you want to label
ge = track.gencode(build = "hg19",gene.collapse = T , labels.suppress.gr = T, height = 3)

# Creates variables containing bigwig files
# Here you can change the color scheme, size, and label of the bigwig tracks
gt.enh = gTrack(H3K27ac_file, name = 'H3K27Ac', bar = TRUE, col = "blue", height = 5)
#gt.DNase = gTrack(DNase_file, name = 'BRD4', bar = TRUE, col = "red", height = 6)
#gt.cn = gTrack(CN_file, name = 'CN', bar = TRUE, col = "red", height = 5)

# Create plot
# EDIT HERE #
# ADD gt.enh, gt.DNase, or gt.cn depending on bigwig tracks being presented
pdf("MCF7_5kb_log10_interaction_H3K27ac.pdf",height = 30,width = 20)
plot(c(ge,gt.enh,dat.gt), wins)
dev.off()

png("MCF7_5kb_log10_interaction_H3K27ac.png",units="in", width=20, height=30, res=300)
plot(c(gt.enh,dat.gt), wins)
dev.off()
