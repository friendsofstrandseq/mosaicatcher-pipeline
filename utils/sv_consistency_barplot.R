## Jun 27, 2018 ADS
## function to generate SV_Consistency barplots from Mosaicatcher outputs
SVplotting <- function(fileLoc, File){
  
  library(cowplot)
  library(ggplot2)
  library(data.table)
  library(GenomicRanges)
  
  if (!file.exists("SV_ConsistencyCheck")) {
    dir.create("./SV_ConsistencyCheck")
  }
  
  # ***********************************************
  ash12rainbow= c("#77AADD", "#4477AA", "#114477", "#CC99BB", "#AA4488", "#771155", "#DDDD77", "#AAAA44", "#777711", "#DDAA77", "#AA7744", "#774411")
  # blues = del; pinks = dup; greens = inv; browns = complex
  
  #################################################   
  ##### Helper FUNCTIONS 1:
  plt_f.a <- function(in_df){
    title <- paste0("CF range: ", min(sort(in_df$cellFrac)), "-", max(sort(in_df$cellFrac)))
    ggplot(in_df) + geom_bar(aes(x=regions, y=value, fill=Var1), stat="identity")+ xlab("Region of predicted SV") + ylab("Total cells with SV (N)") + theme_bw()+
      #theme(axis.text.x = element_text(size=10, angle=90), axis.text.y = element_text(size=4))+   
      scale_fill_manual(values = ash12rainbow, drop=FALSE, name="SV class") + coord_flip()+
      ggtitle(title)
  }
  #tmp <- df[which(df$cellCount < 100 & df$cellCount > 20),]
  #plt_f.a(tmp)
  
  ##### Helper FUNCTIONS 2:
  plt_f.b <- function(in_df){ 
    in_df.2<-as.data.frame(setorder(setDT(in_df), -prop)[, head(.SD, 1), keyby = regions]) # pulls out the 'best-fit' SV type
    title <- paste0("agreement: ", round(min(sort(in_df.2$prop)),1), "-", round(max(sort(in_df.2$prop)),1))
    ggplot(in_df.2) + geom_bar(aes(x=regions, y=as.numeric(prop), fill=Var1), stat="identity") + xlab("Region of predicted SV") + ylab("Best-fit SV called (%)") + theme_bw()+
      #theme(axis.text.x = element_text(angle=90))  + 
      scale_fill_manual(values = ash12rainbow, drop=FALSE, name="SV class")+ coord_flip()+
      ggtitle(title)
  }
  # plt.b<- plt_f.b(tmp)
  
  ##### Helper FUNCTIONS 3:
  plot2x2 <- function(plot.a, plot.b){
    prow <- plot_grid(
      plot.a + theme(legend.position="none"),
      plot.b + theme(legend.position="none", axis.text.y = element_blank(), axis.title.y = element_blank())+geom_hline(yintercept = 50, linetype="dotted"),
      nrow = 1,
      ncol=2,
      rel_widths = c(1, .5)
    )
    legend <- get_legend(plot.a)
    p <- plot_grid(prow, legend, ncol = 2, rel_widths = c(1, .1))
    return(p)
  }
  # plt.a<- plt_f.a(tmp)
  # plt.b<- plt_f.b(tmp)
  # plot2x2(plt.a, plt.b)
  
  #################################################     
  
  #read in data
  input<- GRanges(read.table(paste0(fileLoc,File), header=T, stringsAsFactors = F))
  input$sv_call_name <- factor(input$sv_call_name, levels=c("del_h1",  "del_h2",  "del_hom", "dup_h1",  "dup_h2",  "dup_hom", "inv_h1",  "inv_h2",  "inv_hom", "idup_h1", "idup_h2", "complex"))
  message(paste0(File, "loaded"))
  
  cellNo <- length(unique(input$cell))
  cellCount<- table(as.factor(input))
  regions<- split(input, as.factor(input)) # split data into the regions (SV locations)
  
  # calculate the number of each SV class at each region 
  svType <- sapply(1:length(regions), function(x) table(regions[[x]]$sv_call_name)) 
  
  # prepare df for plotting
  df<- melt(svType) # var2=region
  df$regions <- unlist(names(regions))[df$Var2] # add region information
  df$cellCount <- cellCount[df$Var2] # add No of cells with SV called at region
  df$cellFrac <- round(df$cellCount/cellNo, 2) # calculate cell fraction
  df$prop <- round((df$value/df$cellCount)*100, 4) # calculate proportion of each SV_type
  # head(df)
  
  message("dataframe generated")
  # ***********************************************
  ## PLOT 1
  ### Scatter plot of SV consistency: plot of cell fraction vs SV type
  df.2<-as.data.frame(setorder(setDT(df), -prop)[, head(.SD, 1), keyby = regions]) # orders based on prop and pulls out TOP (#1) hit
  
  scatterPlt<- ggplot(df.2) +
    #geom_point(aes(x=as.numeric(cellCount), y=as.numeric(prop), col=Var1), size=2)+ xlab("number of cells with SV (N)") + 
    geom_point(aes(x=as.numeric(cellFrac), y=as.numeric(prop), col=Var1), size=3, alpha=0.5)+ xlab("fraction of cells with SV (VAF)") + 
    ylab("level of agreement for SV type (%)") + theme_bw()+
    scale_color_manual(values = ash12rainbow, drop=FALSE) 
  # scatterPlt ## CURRENTLY NOT SAVED
  
  # ***********************************************
  ## PLOT 2
  ### complex plots of total SV types called in each region
  
  breaks <- quantile(df[df$cellCount >1,]$cellCount) # used to divide into "high-rare" SVs
  
  ## GENERATE PLOTS based on breaks
  ## PRODUCE 2x2 PLOTs AND SAVE
  sv_high <- df[which(df$cellCount > breaks[4]),]
  plt_high.a <- plt_f.a(sv_high)
  plt_high.b <- plt_f.b(sv_high)
  plt.high.save<- plot2x2(plt_high.a, plt_high.b)
  ggsave(plt.high.save, file=paste0("./SV_ConsistencyCheck/", File, "_highCF_barplot.pdf" ), width=15, height=(length(unique(sv_high$regions)))*0.15, onefile = T)
  
  sv_med <- df[which(df$cellCount <= breaks[4] & df$cellCount > breaks[3]),]
  plt_med.a <- plt_f.a(sv_med)
  plt_med.b <- plt_f.b(sv_med)
  plt.med.save<- plot2x2(plt_med.a, plt_med.b)
  ggsave(plt.med.save, file=paste0("./SV_ConsistencyCheck/", File, "_medCF_barplot.pdf" ), width=15, height=(length(unique(sv_med$regions)))*0.15, onefile = T)
  
  sv_low <- df[which(df$cellCount <=  breaks[3] & df$cellCount >  breaks[2]),]
  plt_low.a <- plt_f.a(sv_low)
  plt_low.b <- plt_f.b(sv_low)
  plt.low.save<- plot2x2(plt_low.a, plt_low.b)
  ggsave(plt.low.save, file=paste0("./SV_ConsistencyCheck/", File, "_lowCF_barplot.pdf" ), width=15, height=(length(unique(sv_low$regions)))*0.15, onefile = T)
  
  sv_rare <- df[which(df$cellCount <=  breaks[2] & df$cellCount > 1),]
  plt_rare.a <- plt_f.a(sv_rare)
  plt_rare.b <- plt_f.b(sv_rare)
  plt.rare.save<- plot2x2(plt_rare.a, plt_rare.b)
  ggsave(plt.rare.save, file=paste0("./SV_ConsistencyCheck/", File, "_rareCF_barplot.pdf" ), width=15, height=(length(unique(sv_rare$regions)))*0.15, onefile = T)
  message("pdfs saved successfully")
  
  ############ SAVED as PDFs  
} 
