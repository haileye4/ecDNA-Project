
# Install and load package
#install.packages("RCircos")
library(RCircos)

# Set working directory
setwd("/Users/haileyjohnson/Desktop/ecDNA_Project")

######################################################################
# Input a bedpe file that contains the looping/interaction data of interest
# It is important here to ensure all non-canoncial chromosomes are filtered out (i.e. chrM, etc)

#Translocations <- read.table("BREAST_PCAWG_chr8_chr11_translocations.bedpe", #This I need to make
                             #sep="\t",  quote="", header=F)

Translocations <- MYCandFriends
                             
# These are the headers for the bedpe file. There should only be 6 columns with the chromosome, start position, and stop
# position for both anchors of a chromatin interaction or WGS breakpoints
names(Translocations) <- c("Chromosome", "chromStart", "chromEnd", "Chromosome.1", "chromStart.1", "chromEnd.1")

Translocation <- as.data.frame(Translocations)

# Exclude chromosomes that you are not interested in plotting
chr.exclude <- c("chrY")

# Use this if you want to plot the entire genome
#chr.exclude = NULL

# This inputs the cytobands coloring the different regions of each chromosome
data(UCSC.HG19.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;

# Define the number of tracks you want inside the circle/chromosome ring (I think you can have as many tracks as you wish)
tracks.inside <- 2

# Define the number of tracks you want outside the circle/chromosome ring
tracks.outside <- 0

# Define the circos plot
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

# Begin building circos plot
out.file <- "MYCandFriendsCircos.pdf"
 pdf(file=out.file, height=10, width=10, compress=TRUE)
 RCircos.Set.Plot.Area();
 RCircos.Chromosome.Ideogram.Plot();

# This track is the chromsomal interaction
track.num <- 2;
RCircos.Link.Plot(Translocation, track.num, TRUE, is.sorted=FALSE)
 ribbon.data <- Translocations
 colnames(ribbon.data) <- c("chromA", "chromStartA", "chromEndA", "chromB", "chromStartB", "chromEndB")

# Now you can start introducing additional tracks from bed (tile based) or bedgraph files (different values for 
# separate regions, similar to bigwigs just a less efficient version of those files, use histogram) 
 
# Bed track example
#RCircos.Tile.Data <- read.table("ZR75_CCND1_enhancers.bed", sep="\t",  quote="", header=F)
#track.num <- 1;
#side <- "in";
#RCircos.Tile.Plot(RCircos.Tile.Data, track.num, side);

# Bedgraph track example
#RCircos.Histogram.Data <- read.table("DMS273_CN_chr8.bedgraph", 
                                     #sep="\t",  quote="", header=F)

# hist.colors <- rep("black", nrow(RCircos.Histogram.Data))
#rows <- which(RCircos.Histogram.Data$V4< -1)
#rows_positive <- which(RCircos.Histogram.Data$V4>1)
#hist.colors[rows] <- "blue";
#hist.colors[rows_positive] <- "red";
#RCircos.Histogram.Data["PlotColor"] <- hist.colors;
#data.col <- 4;
#track.num <- 2;
#side <- "in";
# RCircos.Line.Plot(RCircos.Histogram.Data,
                 #  +     data.col, track.num, side);

dev.off()
 
 
