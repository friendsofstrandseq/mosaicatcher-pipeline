#' author David Porubsky, Maryam Ghareghani

#required packages
library("RColorBrewer")
library("gtable")
library("grid")
library("reshape2")
library("ggplot2")
library("gridExtra")

#Helper functions
getlegend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

setPanelHeights <- function(g, heights){
  g$heights <- grid:::unit.list(g$heights)
  id_panels <- unique(g$layout[g$layout$name=="panel", "t"])
  g$heights[id_panels] <- heights
  g
}

#Main heatmap function
plotHeatmapSegment <- function(dataFrame, plot.log=FALSE, file=NULL, aggProbs=F, maximumCN=4, CNV=3) {
  #get rid of leading X character after reading data in
  colnames(dataFrame) <- gsub(colnames(dataFrame), pattern = 'X', replacement = '')
  
  #Filter probs for certain CNVs
  probs.names <- colnames(dataFrame[(maximumCN+9):ncol(dataFrame)])
  probs.names.l <- sapply(probs.names, function(x) strsplit(x, ''))
  probs.names.cnv <- sapply(probs.names.l, function(x) sum(as.numeric(x)))
  probs2filt <- names(probs.names.cnv[probs.names.cnv > CNV])
  dataFrame <- dataFrame[,!colnames(dataFrame) %in% probs2filt]
  
  probs <- as.matrix(dataFrame[,c(8:ncol(dataFrame))])
  
  #cluter probs hclust
  #get order of rows
  #sort dataFrame rows based on hclust order
  ord <- order.dendrogram(as.dendrogram(hclust(dist(probs, method = "euclidean"), method = "ward.D")))
  dataFrame <- dataFrame[ord,]
  probs <- probs[ord,]
  
  #sort input data.frame
  rowIds <- dataFrame$cells
  if (aggProbs) {
    rowIds <- c(rowIds[which.max(rowIds)], rowIds[-which.max(rowIds)])
    dataFrame$cells <- factor(dataFrame$cells, levels=rowIds)
  } else {
    dataFrame$cells <- factor(dataFrame$cells, levels=rowIds)
  }
  
  #tranform wide table format into a long format for plotting
  #tab.long <- melt(dataFrame, id.vars=c('cells', 'types', 'Wcount', 'Ccount','chr'), measure.vars=c("CN0","CN1","CN2","CN3","CN4","CN5","X00","X01","X10","X02","X11","X20","X03","X12","X21","X30","X04","X13","X22","X31","X40","X05","X14","X23","X32","X41","X50"))
  tab.long <- melt(dataFrame, id.vars=c('chr','start','end','cells', 'types', 'Wcount', 'Ccount'), measure.vars=colnames(probs))
  
  #set theme for the main heatmap
  heatmap_theme <- theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.y=element_blank(),
    plot.margin = unit(c(-0.5,-0.5,-0.5,-0.5),"mm"),
    legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
    legend.position="bottom"
  )
  
  #plot the main heatmap
  if (plot.log) {
    plt <- ggplot(tab.long) + geom_tile(aes(x=variable, y=cells, fill=as.numeric(value))) + heatmap_theme + scale_fill_gradientn(colours =c("white","blue","red"), name="", trans='log')
  } else {
    plt <- ggplot(tab.long) + geom_tile(aes(x=variable, y=cells, fill=as.numeric(value))) + heatmap_theme + scale_fill_gradient(low = "white", high = "red", name="")
  }
  
  #set the theme for description columns
  header_theme <- theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.length = unit(0,"null"),
    plot.margin = unit(c(-0.5,-0.5,-0.5,-0.5),"mm"),
    legend.margin=margin(t=0, r=0, b=0, l=0, unit="cm"),
    legend.position="bottom"
  )
  
  #set the colors for the upper description row
  colColors <- brewer.pal(n=maximumCN+1, name="Set1")
  names(colColors) <- paste0("CN",0:maximumCN) #c("CN0","CN1","CN2","CN3","CN4","CN5")
  CNV.states <- rep(names(colColors), table(probs.names.cnv))
  CNV.states <- CNV.states[gsub(CNV.states, pattern = 'CN', replacement = '') <= CNV]
  colAnnot.df <- data.frame(ID=factor(levels(tab.long$variable), levels=levels(tab.long$variable)), type = c(names(colColors), CNV.states))
  
  #plot the upper description row
  header <- ggplot(colAnnot.df) + geom_tile(aes(x=ID, y=1, fill=type)) + scale_fill_manual(values = colColors) + header_theme + guides(fill = guide_legend(nrow = 1)) + ggtitle(paste0(dataFrame[1,1], "_", dataFrame[1,2], "_", dataFrame[1,3]))
  header.dummy <- ggplot(data.frame(size=1)) + geom_tile(aes(x=size, y=1), fill='white') + header_theme + theme(legend.position="none")
  
  #plot the right side description column
  celltypes <- as.character(dataFrame$types)
  
  if (aggProbs) {
    celltypes <- c(celltypes[which(celltypes=='all')], celltypes[-which(celltypes=='all')])
  }  
  colType.df <- data.frame(ID=celltypes, level=c(1:length(celltypes)))
  cellType <- ggplot(colType.df) + geom_tile(aes(x=1, y=factor(level), fill=ID)) + scale_fill_manual(values = c('cc'="paleturquoise4",  'cw'="blue", 'wc'="olivedrab",'ww'="sandybrown", 'all'="red")) + header_theme
  
  #extract legends from the plots
  plt.leg <- getlegend(plt)
  header.leg <- getlegend(header)
  cellType.leg <- getlegend(cellType)
  
  #organize plots into a single row
  legends <- gtable_cbind(plt.leg, header.leg, cellType.leg)
  
  #remove legends from the plots
  plt <-  plt + theme(legend.position='none')
  header <- header + theme(legend.position='none')
  cellType <- cellType + theme(legend.position='none')
  
  #extract ggplot table
  g1 <- ggplotGrob(header)
  g2 <- ggplotGrob(header.dummy)
  g3 <- ggplotGrob(plt)
  g4 <- ggplotGrob(cellType)
  
  #connect plot components by row
  g.col1 <- rbind(g1, g3, size = "last")
  g.col2 <- rbind(g2, g4, size = "last")
  
  #set the relative heigth of plot elements
  g.col1 <- setPanelHeights(g.col1, unit.c(unit(1,"line"), unit(nrow(g3),"line")))
  g.col2 <- setPanelHeights(g.col2, unit.c(unit(1,"line"), unit(nrow(g4),"line")))
  
  #bind plot columns into a plot matrix
  g.matrix <- gtable_cbind(g.col1, g.col2)
  
  #set relative widths of the final plot matrix
  g.matrix$widths[c(4,11)] <- unit.c(unit(ncol(g3),"null"), unit(1, "line")) 
  
  #add legend rows at the bottom of the final plot
  final.plot <- arrangeGrob(g.matrix, legends, ncol = 1)
  #set the position of the legend
  final.plot$heights <- unit.c(unit(14,"line"), unit(1, "line")) 
  
  #grid.newpage()
  #grid.draw(final.plot)
  return(list(heatmap.plt=final.plot, hc.ord=rev(dataFrame$cells)))
}


hapProbsFile = snakemake@input[["haplotypeProbs"]]
GTprobsFile = snakemake@input[["genotypeProbs"]]

GTprobs <- read.table(GTprobsFile, stringsAsFactors = F, header = T)
hapProbs <- read.table(hapProbsFile, stringsAsFactors = F, header = T)

GTprobs.l <- split(GTprobs, paste(GTprobs$chr, GTprobs$start, GTprobs$end))
hapProbs.l <- split(hapProbs, paste(hapProbs$chr, hapProbs$start, hapProbs$end))

pdf(snakemake@output[[1]])
lapply(1:length(GTprobs.l), function(x) grid.arrange(plotHeatmapSegment(GTprobs.l[[x]])$heatmap.plt, plotHeatmapSegment(hapProbs.l[[x]])$heatmap.plt))
dev.off()
