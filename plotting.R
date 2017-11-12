#' Plot reads based on directionality
#'
#' @param input.reads A \code{\link[GenomicRanges]{GRanges}} object with Strand-specific read data
#' @param bin.size An \code{integer} (e.g. 2000) used to bin data for plotting
#' @param col_a Assigns mreads (i.e. Watson) col (e.g "grey", or rgb(1,2,3, max=255))
#' @param col_b Assigns preads (i.e. Crick) col (e.g "red", or rgb(1,2,3, max=255))
#' @author Ashley Sanders
#' @export
#' 

read.plotR <- function(input.reads, bin.size=2000, col_a=rgb(243,165,97, max=255), col_b=rgb(103,139,139, max=255))
{  
  
  # bin the data for plotting as histogram
  chr<- as.character(seqnames(input.reads)[1])
  max_pos <- max(start(input.reads))
  numbins <- floor(max_pos/bin.size) #calculate number of bins
  ir <- successiveIRanges(rep(bin.size,numbins)) #create continuous ranges based on number of bins and the binsize
  lastBinWidth <- max_pos - (bin.size*numbins)
  lastBinStart <- (max_pos - lastBinWidth) + 1
  lastRange <- IRanges(start = lastBinStart, width = lastBinWidth) #calculate last range
  ir <- c(ir,lastRange)
  #rows <- rep(0,length(ir))
  gr <- GRanges(seqnames=chr, ranges=ir) #initialize GRanges object to store bin read counts
  
  gr$midpoint <- start(gr) + (width(gr)/2)  #get the midposition of each bin
  
  #counts overlaps between bins and reads in input.reads
  gr$mreads <- GenomicRanges::countOverlaps(gr, input.reads[strand(input.reads)=='-']) 
  gr$preads <- GenomicRanges::countOverlaps(gr, input.reads[strand(input.reads)=='+'])
  gr$readNo <- gr$mreads + gr$preads
  
  ## prepare plots data  
  dfplot.reads <- as.data.frame(gr)
  dfplot.reads$preads <- - dfplot.reads$preads	# negative non-ref (preads) counts
  
  # plotting function 
  ### p1 <- barplot of plus/minus reads
  p1 <- ggplot(dfplot.reads) + geom_rect(aes(ymin=0, ymax=mreads, xmin=start, xmax=end), fill=col_a, size=0.5)+ # minus (Watson) reads on top
    scale_fill_identity() +
    scale_x_continuous(limits=c( min(start(input.reads)), max(end(input.reads)) )) + #, breaks=NULL, labels(NULL))    
    #scale_y_continuous(limits=c( min(start(input.reads)), max(end(input.reads)) )) + #, breaks=NULL, labels(NULL))    
    geom_rect(aes(ymin=0, ymax=preads, xmin=start, xmax=end), fill=col_b, size=0.5) ### + ylab("depth") # + labs(x=chr) # plus (Crick) reads on btm
  p1 <- p1  #+ theme_light()
  
  
  out<-list(p1)
  
  return(out)
}


#' Create seperate PDF files for different chromosomes and in each file, create plots for W, C, and total read coverage in all cells with the fitted NB distribution.
#' 
#' @inheritParams read.bams
#' @inheritParams estimateR
#' @inheritParams read.plotR
#' @param cellTypes TODO ...
#' @param dispersion The estimated dispersion parameters.
#' @author Maryam Ghareghani
#' @export
#' 

NBfitplots <- function(directory, counts, cellTypes, p, dispersion, bin.size)
{
  #setwd(directory)
  numCells = (ncol(counts[[1]])-3)/2
  for (i in 1:length(counts))
  {
    pdf(paste0(directory,"chr", i, "NBfit.pdf"))
    par(mfrow=c(2,3))
    for (j in 1:numCells)
    {
      Wcount = counts[[i]][,2*j+2]
      Ccount = counts[[i]][,2*j+3]
      
      r = dispersionPar(cellTypes[i,j], dispersion[i,j], binLength = bin.size)
      r = c(r, dispersion[i,j]) # adding the dispersion par for the total number of read counts
      
      allCounts = list(Wcount, Ccount, Wcount + Ccount)
      name = c("Watson", "Crick", "Total")
      
      for (q in 1:3)
      {
        ymax = max(table(allCounts[[q]]))/sum(table(allCounts[[q]]))
        hist(allCounts[[q]], breaks = (-1:max(allCounts[[q]]))+0.5, freq = FALSE, xlab = "", ylab = "", 
             main = paste("cell", j, name[q], ", t =", cellTypes[i,j]), ylim = c(0, ymax))
        x = 0:max(allCounts[[q]])
        v = dnbinom(x, size = r[q], prob = p)
        lines(x,v,col="red")
      }
    }
    dev.off()
  }
}




#' Plot probability table as a heatmap.
#'
#' @param dataFrame A \code{\link{data.frame}} object that containing genotype probabilities.
#' @param plot.log A logical indicating whether or not to plot in logarithmic scale.
#' @param file A file to export the plot to.
#' @param aggProbs A logical indicating whether or not to plot aggregate probability values.
#' @param CNV A copy number value until which the probability values are plotted.
#' @author David Porubsky
#' @export
#' 

plotHeatmapSegment <- function(dataFrame, plot.log=FALSE, file=NULL, aggProbs=FALSE, CNV=3) {

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

  #get rid of leading X character after reading data in
  colnames(dataFrame) <- gsub(colnames(dataFrame), pattern = 'X', replacement = '')
  
  #Filter probs for certain CNVs
  probs.names <- colnames(dataFrame[14:ncol(dataFrame)])
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
  colColors <- brewer.pal(n=6, name="Set1")
  names(colColors) <- c("CN0","CN1","CN2","CN3","CN4","CN5")
  CNV.states <- rep(c("CN0","CN1","CN2","CN3","CN4","CN5"), table(probs.names.cnv))
  CNV.states <- CNV.states[gsub(CNV.states, pattern = 'CN', replacement = '') <= CNV]
  colAnnot.df <- data.frame(ID=factor(levels(tab.long$variable), levels=levels(tab.long$variable)), type = c("CN0","CN1","CN2","CN3","CN4","CN5", CNV.states))
  
  #plot the upper description row 
  header <- ggplot(colAnnot.df) + geom_tile(aes(x=ID, y=1, fill=type)) + scale_fill_manual(values = colColors) + header_theme + guides(fill = guide_legend(nrow = 1))
  header.dummy <- ggplot(data.frame(size=1)) + geom_tile(aes(x=size, y=1), fill='white') + header_theme + theme(legend.position="none")
  
  #plot the right side description column
  celltypes <- as.character(dataFrame$types)
  
  if (aggProbs) {
    celltypes <- c(celltypes[which(celltypes=='all')], celltypes[-which(celltypes=='all')])
  }  
  colType.df <- data.frame(ID=celltypes, level=c(1:length(celltypes)))
  cellType <- ggplot(colType.df) + geom_tile(aes(x=1, y=factor(level), fill=ID)) + scale_fill_manual(values = c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown", 'all'="red")) + header_theme

  #extract legends from the plots
  plt.leg <- getlegend(plt)
  header.leg <- getlegend(header)
  cellType.leg <- getlegend(cellType)
  
  #organize plots into a single row
  legends <- cbind(plt.leg, header.leg, cellType.leg)
  
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
  g.matrix <- cbind(g.col1, g.col2)
  
  #set relative widths of the final plot matrix
  g.matrix$widths[c(4,11)] <- unit.c(unit(ncol(g3),"null"), unit(1, "line")) 
  
  #add legend rows at the bottom of the final plot
  final.plot <- arrangeGrob(g.matrix, legends, ncol = 1)
  #set the position of the legend
  final.plot$heights <- unit.c(unit(14,"line"), unit(1, "line")) 
  
  grid.newpage()
  grid.draw(final.plot)
}

