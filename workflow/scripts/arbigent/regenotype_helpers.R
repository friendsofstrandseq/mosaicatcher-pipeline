# Make a table
make_table_sc_separated <- function(pg_f){
  library(reshape2)
  # library(reshape)
    
  n_sites = dim(unique(pg_f[,c('chrom','start','end')]))[1]
  
  haps_to_consider = c('1010','0101','1001','0110', '1000', '0010', '0000', '1110','1011', '1111')

  # Cut down to interesting haplotypes
  pgi2 = pg_f[pg_f$haplotype %in% haps_to_consider,]
  
  supertable = cast(pgi2, formula = cell + chrom + start +end + class + expected + W + C ~ haplo_name, value='logllh')

  top_pred = rep('None', dim(supertable)[1])
  top_score = rep(1, dim(supertable)[1])
  second_pred = rep('None', dim(supertable)[1])
  second_score = rep(1, dim(supertable)[1])
  for (row in 1:dim(supertable)[1]){
    #print(colnames(
      relevant_cols = supertable[row,9:dim(supertable)[2]]
      
      # Top
      best_names = colnames(relevant_cols[which(relevant_cols == max(relevant_cols))])
      top_pred[row] = paste(best_names, collapse='/')
      top_score[row] = as.numeric(relevant_cols[which(relevant_cols == max(relevant_cols))])[1]
      # Second
      remaining_cols = relevant_cols[which(relevant_cols != max(relevant_cols))]
      second_best_names = colnames(remaining_cols[which(remaining_cols == max(remaining_cols))])
      second_pred[row] = paste(second_best_names, collapse='/')
      second_score[row] = as.numeric(remaining_cols[which(remaining_cols == max(remaining_cols))])[1]
      
  }
  #supertable$top_score = top_score
  #supertable$second_score =second_score
  supertable$top_pred = top_pred
  supertable$second_pred = second_pred
  supertable$llr_1st_to_2nd =  round(log10((10**top_score) / (10**second_score)),5)
  supertable$llr_1st_to_ref =  round(top_score,5)
  
  supertable$expected = round(supertable$expected, 5)

  return(supertable)

  }

# Whoeps, 31th May 2020
# Admittedly the name 'standalone' is not wisely chosen. Anyway these are functions needed for the 
# three outputs of standalone. R

#' Function A: make a dumbbell plot
#'
#' @param call_llhs
#' @param 
#' @author Wolfram Hoeps
#' @export
make_dumbbell <- function(segs_llhs_f, groupname='Inv', run_shiny=F){
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  # take the segs_llhs (which are call_llhs), and add important information:
  # maxval, secval, len, maxname, interesting_minval and sum_hetinv_max
  segs_llhs = attach_max_sec_name_to_call_llhs(segs_llhs_f)
  
  # I think len is no longer needed but I could be wrong.
  segs_llhs$len = segs_llhs$end - segs_llhs$start
  segs_llhs$interesting_minval = sign(segs_llhs$interesting_minval) * log10(abs(segs_llhs$interesting_minval)+1)
  segs_llhs$maxval = sign(segs_llhs$maxval) *log10(abs(segs_llhs$maxval)+1)
  segs_llhs$sum_hetinv_max = sign(segs_llhs$sum_hetinv_max) *log10(abs(segs_llhs$sum_hetinv_max)+1)
  segs_llhs$sum_0101 = sign(segs_llhs$sum_0101) *log10(abs(segs_llhs$sum_0101)+1)
  segs_llhs$maxname_short = substr(segs_llhs$maxname,5,9)
  
  # SORT
  segs_llhs=segs_llhs[order(segs_llhs$len),]
  rownames(segs_llhs) <- NULL
  ### PLOT ###
  g <- ggplot(segs_llhs, aes(x=(len/1000))) +scale_fill_brewer(palette = "Spectral")
  g = g  +    
    # Add a few grid lines first of all
    geom_hline(yintercept=log1p(0), alpha=1) +
    geom_hline(yintercept=2, alpha=0.2) +
    geom_hline(yintercept=-2, alpha=0.2) +
    
    # Add the vertical lines that go from the lowest to the highest value for each segment
    geom_segment(aes(xend=(len/1000), y=interesting_minval, yend=maxval), 
                 size=0.5,  colour="grey70") + 
    
    # Add the dots for HOM, HET and HIGHEST. They will all lie on the vertial lines we just created
    #geom_point(aes(y=maxval, colour='score OTHER', group=segs_llhs$start, text=paste0(groupname,' #',row.names(segs_llhs),'\n',segs_llhs$maxname_short)), 
    #           alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[1]) +
    geom_point(aes(y=maxval, colour='score OTHER', group=segs_llhs$start, text=paste0('start:',segs_llhs$start,'\n ',segs_llhs$maxname_short)),
                      alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[1]) +
    geom_point(aes(y=sum_hetinv_max, colour='score HET',group=segs_llhs$start),
               alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[3]) +
    geom_point(aes(y=sum_0101, colour='score HOM',group=segs_llhs$start),
               alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[2]) +

    # Change colors. Entries are sorted alphabetically, so this is a bit a mess. 
    # alphabet: [het, hom, other], hue_pal: red, green, blue
    # desired order: het - green, other - red, hom - blue
    # therefore: c(2,3,1)
    scale_color_manual(values=hue_pal()(3)[c(2,3,1)]) +
    
    # Add labels
    labs(title="Inversion re-genotyping", 
         subtitle="") +
    xlab('SV length [kb]') +
    ylab('LLR over REF') +
    
    # Make a nice x axis with log and logticks
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x)
      #limits=c(1,5000)
      #labels = scales::trans_format("log10", scales::x)
    ) +
    # scale_y_log10(
    #   breaks = scales::trans_breaks("log10", function(y) 10^(y+1))
    #   #limits=c(1,5000)
    #   #labels = scales::trans_format("log10", scales::x)
    # ) +
    annotation_logticks(sides = 'b',colour = 'gray20' ) +
    scale_y_continuous(breaks = (seq(-4,4,by=1))) +
    # Remove title for all legends
    #theme(legend.title=element_blank()) +
    
    # Add text for other
    #geom_text(aes(y= maxval, label = paste0('start:',segs_llhs$start,'\n bestpred',segs_llhs$maxname_short)),
    #          size = 2,
    #           hjust = "left") +
#    geom_text(aes(y= maxval, label = paste0(groupname, ' #',row.names(segs_llhs),'\n',segs_llhs$maxname_short)),
#              size = 2,
#              hjust = "left") +


    
    # EEEEEHm, not sure. 
    guides(color = guide_legend(reverse = TRUE))
  
  # Okay, we created 'g' now. We CAN pass it to shiny if we want to
  if (run_shiny){
    # embed in plotly #
    ui <- fluidPage(
      plotlyOutput("distPlot")
    )
    server <- function(input, output) {
      output$distPlot <- renderPlotly({
        g
      })
    }
    # Run a shiny app
    shinyApp(ui = ui, server = server)
  }
  
  return(g)
}

#' Function A.2: make a dumbbell plot but only hom and het
#'
#' @param call_llhs
#' @param 
#' @author Wolfram Hoeps
#' @export
make_dumbbell_2 <- function(segs_llhs_f, groupname='Inv', run_shiny=F){
  library(ggplot2)
  library(plotly)
  library(shiny)
  library(RColorBrewer)
  library(scales)
  # take the segs_llhs (which are call_llhs), and add important information:
  # maxval, secval, len, maxname, interesting_minval and sum_hetinv_max
  segs_llhs = attach_max_sec_name_to_call_llhs(segs_llhs_f)
  
  # I think len is no longer needed but I could be wrong.
  segs_llhs$len = segs_llhs$end - segs_llhs$start
  segs_llhs$interesting_minval = sign(segs_llhs$interesting_minval) * log10(abs(segs_llhs$interesting_minval)+1)
  segs_llhs$maxval = sign(segs_llhs$maxval) *log10(abs(segs_llhs$maxval)+1)
  segs_llhs$sum_hetinv_max = sign(segs_llhs$sum_hetinv_max) *log10(abs(segs_llhs$sum_hetinv_max)+1)
  segs_llhs$sum_0101 = sign(segs_llhs$sum_0101) *log10(abs(segs_llhs$sum_0101)+1)
  segs_llhs$maxname_short = substr(segs_llhs$maxname,5,9)
  
  # SORT
  segs_llhs=segs_llhs[order(segs_llhs$len),]
  rownames(segs_llhs) <- NULL
  
  # Get the max value of hom and het. This might be not very elegant but WELL
  df = segs_llhs[,c('sum_0101', 'sum_hetinv_max')]
  segs_llhs[,'hardpriormax'] = apply(df,1,max)
  
  ### PLOT ###
  g <- ggplot(segs_llhs, aes(x=(len/1000))) +scale_fill_brewer(palette = "Spectral")
  g = g  +    
    # Add a few grid lines first of all
    geom_hline(yintercept=log1p(0), alpha=1) +
    geom_hline(yintercept=2, alpha=0.2) +
    geom_hline(yintercept=-2, alpha=0.2) +
    
    # Add the vertical lines that go from the lowest to the highest value for each segment

  
    geom_segment(aes(xend=(len/1000), y=interesting_minval, yend=hardpriormax), 
                 size=0.5,  colour="grey70") + 
    
    # Add the dots for HOM, HET and HIGHEST. They will all lie on the vertial lines we just created
    #geom_point(aes(y=maxval, colour='score OTHER', group=segs_llhs$start, text=paste0(groupname,' #',row.names(segs_llhs),'\n',segs_llhs$maxname_short)), 
    #           alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[1]) +
    geom_point(aes(y=sum_hetinv_max, colour='score HET',group=segs_llhs$start),
               alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[3]) +
    geom_point(aes(y=sum_0101, colour='score HOM',group=segs_llhs$start),
               alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[2]) +
    
    # Change colors. Entries are sorted alphabetically, so this is a bit a mess. 
    # alphabet: [het, hom, other], hue_pal: red, green, blue
    # desired order: het - green, other - red, hom - blue
    # therefore: c(2,3,1)
    scale_color_manual(values=hue_pal()(3)[c(2,3)]) +
    
    # Add labels
    labs(title="Inversion re-genotyping", 
         subtitle="") +
    xlab('SV length [kb]') +
    ylab('LLR over REF') +
    
    # Make a nice x axis with log and logticks
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x)
      #limits=c(1,5000)
      #labels = scales::trans_format("log10", scales::x)
    ) +
    # scale_y_log10(
    #   breaks = scales::trans_breaks("log10", function(y) 10^(y+1))
    #   #limits=c(1,5000)
    #   #labels = scales::trans_format("log10", scales::x)
    # ) +
    annotation_logticks(sides = 'b',colour = 'gray20' ) +
    scale_y_continuous(breaks = (seq(-4,4,by=1))) +
    # Remove title for all legends
    #theme(legend.title=element_blank()) +
    

    # EEEEEHm, not sure. 
    guides(color = guide_legend(reverse = TRUE))
  
  # Okay, we created 'g' now. We CAN pass it to shiny if we want to
  if (run_shiny){
    # embed in plotly #
    ui <- fluidPage(
      plotlyOutput("distPlot")
    )
    server <- function(input, output) {
      output$distPlot <- renderPlotly({
        g
      })
    }
    # Run a shiny app
    shinyApp(ui = ui, server = server)
  }
  
  return(g)
}

#' Function A.3: make a dumbbell plot but with probabilities
#'
#' @param call_llhs
#' @param 
#' @author Wolfram Hoeps
#' @export
make_dumbbell_probs <- function(segs_llhs_f, groupname='Inv', run_shiny=F){
  library(ggplot2)
  library(plotly)
  library(shiny)
  library(RColorBrewer)
  library(scales)
  # take the segs_llhs (which are call_llhs), and add important information:
  # maxval, secval, len, maxname, interesting_minval and sum_hetinv_max
  segs_llhs = attach_max_sec_name_to_call_llhs(segs_llhs_f)
  
  # I think len is no longer needed but I could be wrong.
  segs_llhs$len = segs_llhs$end - segs_llhs$start
  #segs_llhs$interesting_minval = sign(segs_llhs$interesting_minval) * log10(abs(segs_llhs$interesting_minval)+1)
  #segs_llhs$maxval = sign(segs_llhs$maxval) *log10(abs(segs_llhs$maxval)+1)
  #segs_llhs$sum_hetinv_max = sign(segs_llhs$sum_hetinv_max) *log10(abs(segs_llhs$sum_hetinv_max)+1)
  #segs_llhs$sum_0101 = sign(segs_llhs$sum_0101) *log10(abs(segs_llhs$sum_0101)+1)
  segs_llhs$maxname_short = substr(segs_llhs$maxname,5,9)
  
  # SORT
  segs_llhs=segs_llhs[order(segs_llhs$len),]
  rownames(segs_llhs) <- NULL
  ### PLOT ###
  g <- ggplot(segs_llhs, aes(x=(len/1000))) +scale_fill_brewer(palette = "Spectral")
  g = g  +    
    # Add a few grid lines first of all
    geom_hline(yintercept=0, alpha=1) +
    geom_hline(yintercept=0.5, alpha=0.2) +
    geom_hline(yintercept=1, alpha=0.2) +
    
    # Add the vertical lines that go from the lowest to the highest value for each segment
    geom_segment(aes(xend=(len/1000), y=interesting_minval, yend=maxval), 
                 size=0.5,  colour="grey70") + 
    
    # Add the dots for HOM, HET and HIGHEST. They will all lie on the vertial lines we just created
    #geom_point(aes(y=maxval, colour='score OTHER', group=segs_llhs$start, text=paste0(groupname,' #',row.names(segs_llhs),'\n',segs_llhs$maxname_short)), 
    #           alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[1]) +
    geom_point(aes(y=maxval, colour='score OTHER', group=segs_llhs$start, text=paste0('start:',segs_llhs$start,'\n ',segs_llhs$maxname_short)),
               alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[1]) +
    geom_point(aes(y=sum_hetinv_max, colour='score HET',group=segs_llhs$start),
               alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[3]) +
    geom_point(aes(y=sum_0101, colour='score HOM',group=segs_llhs$start),
               alpha=1, size=2.5) +#, col=brewer.pal(n=3, name='Set1')[2]) +
    geom_point(aes(y=sum_1010, colour='score REF',group=segs_llhs$start),
               alpha=1, size=2.5) +
    # Change colors. Entries are sorted alphabetically, so this is a bit a mess. 
    # alphabet: [het, hom, other], hue_pal: red, green, blue
    # desired order: het - green, other - red, hom - blue
    # therefore: c(2,3,1)
    scale_color_manual(values=hue_pal()(4)[c(2,3,1,4)]) +
    
    # Add labels
    labs(title="Inversion re-genotyping", 
         subtitle="") +
    xlab('SV length [kb]') +
    ylab('Probability') +
    
    # Make a nice x axis with log and logticks
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x)
      #limits=c(1,5000)
      #labels = scales::trans_format("log10", scales::x)
    ) +
    # scale_y_log10(
    #   breaks = scales::trans_breaks("log10", function(y) 10^(y+1))
    #   #limits=c(1,5000)
    #   #labels = scales::trans_format("log10", scales::x)
    # ) +
    annotation_logticks(sides = 'b',colour = 'gray20' ) +
    #scale_y_continuous(breaks = (seq(-4,4,by=1))) +
    # Remove title for all legends
    #theme(legend.title=element_blank()) +
    
    # Add text for other
    #geom_text(aes(y= maxval, label = paste0('start:',segs_llhs$start,'\n bestpred',segs_llhs$maxname_short)),
    #          size = 2,
    #           hjust = "left") +
    #    geom_text(aes(y= maxval, label = paste0(groupname, ' #',row.names(segs_llhs),'\n',segs_llhs$maxname_short)),
    #              size = 2,
    #              hjust = "left") +
  
  
  
  # EEEEEHm, not sure. 
  guides(color = guide_legend(reverse = TRUE))
  
  # Okay, we created 'g' now. We CAN pass it to shiny if we want to
  if (run_shiny){
    # embed in plotly #
    ui <- fluidPage(
      plotlyOutput("distPlot")
    )
    server <- function(input, output) {
      output$distPlot <- renderPlotly({
        g
      })
    }
    # Run a shiny app
    shinyApp(ui = ui, server = server)
  }
  
  return(g)
}



#' Function C: save a beeswarm plot for each segment
#'
#' @param call_llhs
#' @param 
#' @author Wolfram Hoeps
#' @export
save_beeswarms <- function(pg_f, call_llhs_f, outdir, testrun=F, compositemode=F){

  if (compositemode){
    outdir = paste0(outdir, 'bulk_')
  }
  suppressMessages(dir.create(outdir))
  suppressMessages(dir.create(paste0(outdir, 'beeplots')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/0101')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/0110')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/1001')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/1010')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/else')))
  suppressMessages(dir.create(paste0(outdir, 'beeplots/reject')))
  
  only_main_haps=F
  
  # If compositemode, we want larger dots. 
  dotsize_l = 0.5 + (1.5*compositemode)
  dotsize_s = 0.2 + (1.8*compositemode)
  
  # also outdir has to change then. 

  # Iterate over all segments
  if (testrun==T){
    print('Running beewarm in testmode. Only first 10 segments considered.')
    n_pics = 3
  } else {
    n_pics = dim(call_llhs_f)[1]
  }
  for (n in 1:n_pics){
    
    # Extract the information of for this segment
    pgi = suppressMessages(select_group(pg_f, n))
    # Return a list with the most likely classifications
    topnames_hap_presort = get_top_scoring_haplotypes_standalone(call_llhs_f, pgi, 70)
    # ...and extract the absolute winner FIRST
    most_likely_pred = topnames_hap_presort[1]
    # ... and now sort it
    topnames_hap = c(c('0101','1001','0110'), topnames_hap_presort[!topnames_hap_presort %in% c('0101', '1001','0110')])
    # ... and return the 30 top ones
    topnames_hap = topnames_hap[1:50]
    ## In order to make the 'winner cell' larger, we need to do some magic here. ##
    # pgi2 are all data
    pgi2 = pgi[pgi$haplotype %in% topnames_hap,]
    # pgi3 are SV llhs of only the winning SV per cell.
    pgi3 = pgi %>% group_by(cell) %>% top_n(1, logllh)
    pgi3 = pgi3[pgi3$haplotype %in% topnames_hap,]
    
    ## magic end ##
    
    #### PLOT ####
    g <- ggplot(data=pgi2, aes(x=haplotype, y=logllh)) +  
      # Plot first the fat dots
      geom_beeswarm(data=pgi3,
                    size=dotsize_l,
                    cex=0.2,
                    aes(color=class)) +
      # ... then add the slim ones
      geom_beeswarm(data=pgi2,
                    size=dotsize_s,
                    cex=0.1,
                    alpha=0.5,
                    aes(color=class)
      ) +
      # Add median bars
      stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, 
                   geom="crossbar", width=0.7, size=0.1) +
      
      # Sort entries by our topnames list
      scale_x_discrete(limits = topnames_hap) +
      theme(axis.text.x = element_text(angle = 90)) +

      # Add title
      labs(title=paste(pgi2$chrom, ':',pgi2$start, '-',pgi2$end, ' ', ' LEN=',pgi2$len/1000, '  ', ' INV_no=', n, sep=''), subtitle="") +
      
      # Dotted line for 0
      geom_hline(yintercept = 0, linetype="dashed", 
                 color = "black", size=0.3) 
    
    # Save to the outbeedir
    if (most_likely_pred %in% c('1010','0101','0110','1001')){
      ggsave(filename=paste(outdir, 'beeplots/',most_likely_pred,'/', pgi2$chrom[1], ':',pgi2$start[1], '-',pgi2$end[1], '_',round(pgi2$len[1]/1000, 3),'k.png',sep=''), width=30, height=12, units='cm', device='png')
    } else {
      if (test_for_invrejection(call_llhs_f, pgi2$chrom[1], pgi2$start[1])){
        # If TRUE was returned, it means the inversion can be rejected,
        # i.e. HOM and HET are both less likely than reference. 
        ggsave(filename=paste(outdir, 'beeplots/','reject','/', pgi2$chrom[1], ':',pgi2$start[1], '-',pgi2$end[1], '_',round(pgi2$len[1]/1000, 3),'k.png',sep=''), width=30, height=12, units='cm', device='png')
      } else {
        # These are all segments that have no strongest HOM or HET fighter, but also HOM and HET are not categorically out of question yet. This will be a mixture of light green, red and grey.
        ggsave(filename=paste(outdir, 'beeplots/','else','/', pgi2$chrom[1], ':',pgi2$start[1], '-',pgi2$end[1], '_',round(pgi2$len[1]/1000, 3),'k.png',sep=''), width=30, height=12, units='cm', device='png')
      }
    }
    print(paste0('Saved number ', n))
  }
}


#' Function B.2: return a final overview table
#'
#' @param call_llhs
#' @author Wolfram Hoeps
#' @export
make_table_finaledition <- function(call_llhs_f, group_f, sname_f){
  

  maxN2 <- function(x, N=2){
    len <- length(x)
    if(N>len){
      warning('N greater than length(x).  Setting N=length(x)')
      N <- length(x)
    }
    sort(x,partial=len-N+1)[len-N+1]
  }
  # take the segs_llhs (which are call_llhs), and add important information:
  # maxval, secval, len, maxname, interesting_minval and sum_hetinv_max
  call_llhs_f2 = attach_max_sec_name_to_call_llhs(call_llhs_f)
  
  # also the length. it's a bit a messy topic.
  call_llhs_f2$len = (call_llhs_f2$end - call_llhs_f2$start)/1000
  
  # add label
  call_llhs_f2$orig_label = group_f
  call_llhs_f2$sample = sname_f
  
  # simply extract the interesting information
  cols_of_interest = c('chrom', 'start', 'end', 'len', 'sample', 'orig_label', 'maxname','secname','maxval', 'secval', 'sum_0101', 'sum_0110', 'sum_1001', 'sum_1010', 'sum_hetinv_max')
  overview = call_llhs_f2[,cols_of_interest]
  
  # sort
  #overview = overview[with(overview, order(maxname, -maxval)), ]
  
  # get hardprior info that are easy to digest #
  dfhard = call_llhs_f2[,c('sum_0101', 'sum_0110', 'sum_1001', 'sum_1010')]
  overview[,'second_hard'] = substr(colnames(dfhard)[apply(dfhard, 1, function(x)  which(x == sort(x, decreasing = TRUE)[2])[1])],5,8)
  overview[,'pred_hard'] = substr(colnames(dfhard)[apply(dfhard,1,which.max)],5,8)

  dfhard$max2 = apply(dfhard , 1, maxN2)
  dfhard$max = apply(dfhard , 1, max)

  overview$confidence_hard_over_second = dfhard$max - dfhard$max2
  
  # get also semihard including dups and stuff #
  dfsoft = call_llhs_f2[,c('sum_1010',  'sum_0000',  'sum_0010',  'sum_1000',
                       'sum_0101',  'sum_0110',  'sum_1001',  'sum_1011',
                       'sum_1110',  'sum_2020',  'sum_1020',  'sum_2010')]
  
  overview[,'pred_soft'] = substr(colnames(dfsoft)[apply(dfsoft,1,which.max)],5,8)
  
  dfsoft$max2 = apply(dfsoft , 1, maxN2)
  dfsoft$max = apply(dfsoft , 1, max)
  
  overview$confidence_soft_over_hard = dfsoft$max - dfhard$max
  
  # also add prediction of no-prior.
  overview$pred_nobias = substr(call_llhs_f2$maxname,5,8)
  overview$confidence_nobias_over_hard = call_llhs_f2$maxval - apply(dfhard , 1, max)
  # Make it more understandable
  overview_return = overview[,c('chrom', 'start', 'end', 'len', 'sample', 'orig_label',
                                'pred_hard', 'pred_soft', 'pred_nobias', 
                                'confidence_hard_over_second', 'confidence_soft_over_hard',
                                'confidence_nobias_over_hard', 'second_hard')]
  #'prediction_nobias', 'confidence_nobias')]
  
  #overview_return$confidence_nobias = round(overview_return$confidence_nobias, 2)
  overview_return$confidence_hard_over_second = round(overview_return$confidence_hard_over_second, 2)
  overview_return$confidence_soft_over_hard = round(overview_return$confidence_soft_over_hard, 2)
  overview_return$confidence_nobias_over_hard = round(overview_return$confidence_nobias_over_hard, 2)
  
  #simplify names
  overview_return[overview_return=='1010'] = '0|0'
  overview_return[overview_return=='0101'] = '1|1'
  overview_return[overview_return=='0110'] = '1|0'
  overview_return[overview_return=='1001'] = '0|1'
  
  overview_return[overview_return=='HOM'] = '1|1'
  overview_return[overview_return=='HET'] = '1|0'
  overview_return[overview_return=='REF'] = '0|0'
  
  ## i like this for debugging sometimes.
  #library(pheatmap)
  #pheatmap(overview_return[,c('confidence_hard_over_second','confidence_soft_over_hard','confidence_nobias_over_hard')], cluster_cols = F)
  return(overview_return)
}




#' Function B: return an overview table
#'
#' @param call_llhs
#' @author Wolfram Hoeps
#' @export
make_table <- function(call_llhs_f, group_f, sname_f){
  
  maxN2 <- function(x, N=2){
    len <- length(x)
    if(N>len){
      warning('N greater than length(x).  Setting N=length(x)')
      N <- length(x)
    }
    sort(x,partial=len-N+1)[len-N+1]
  }
  # take the segs_llhs (which are call_llhs), and add important information:
  # maxval, secval, len, maxname, interesting_minval and sum_hetinv_max
  call_llhs_f2 = attach_max_sec_name_to_call_llhs(call_llhs_f)
  
  # also the length. it's a bit a messy topic.
  call_llhs_f2$len = (call_llhs_f2$end - call_llhs_f2$start)/1000

  # simply extract the interesting information
  cols_of_interest = c('chrom', 'start', 'end', 'len', 'maxname','secname','maxval', 'secval', 'sum_0101', 'sum_0110', 'sum_1001', 'sum_hetinv_max')
  overview = call_llhs_f2[,cols_of_interest]
  
  # Add label for easier postprocessing
  peter_a_mode = F
  hufsah_mode = T
  if (peter_a_mode){
  overview$orig_label = paste0(substring(group_f, 1,1), '|', substring(group_f, 3,3))
  overview$source = substring(group, 5,10)
  
  overview$orig_label[((overview$orig_label) == '0|1')] = '1|0,0|1'
  overview$orig_label[(as.character(overview$orig_label) == '1|0')] = '1|0,0|1'
  }
  
  if (hufsah_mode){
    colnames(overview) = c('chrom','start','end','len','maxname','secondhighestname','max_log_llh_ratio','secondhighest_log_llh_ratio','HOM_log_llh_ratio','0110_log_llh_ratio','1001_log_llh_ratio','HET_log_llh_ratio_max_of_0110_1001')
    return(overview)
    }
  
  overview$sample = sname_f
  
  # round some numbers

  
  # sort
  overview = overview[with(overview, order(maxname, -maxval)), ]

  # get hardprior info that are easy to digest
  overview$ref = 0
  df = overview[,c('sum_0101', 'sum_0110', 'sum_1001', 'ref')]
  overview[,'prediction'] = colnames(df)[apply(df,1,which.max)]
  overview$prediction[overview$prediction == 'sum_0101'] = '1|1'
  overview$prediction[overview$prediction == 'sum_0110'] = '1|0,0|1'
  overview$prediction[overview$prediction == 'sum_1001'] = '1|0,0|1'
  overview$prediction[overview$prediction == 'ref'] = '0|0'
  overview$confidence = apply(df,1,max) 
  
  df$max2 = apply(df , 1, maxN2)
  df$max = apply(df , 1, max)
  
  overview$confidence = df$max - df$max2
  
  # also add prediction of no-prior.
  overview$prediction_nobias = substr(overview$maxname,5,8)
  overview$confidence_nobias = overview$maxval - apply(df , 1, max)
  # Make it more understandable
  overview_return = overview[,c('sample', 'chrom', 'start', 'end', 'len',
                                'source', 'orig_label', 'prediction', 'confidence')]#, 
                                #'prediction_nobias', 'confidence_nobias')]
  
  #overview_return$confidence_nobias = round(overview_return$confidence_nobias, 2)
  overview_return$confidence = round(overview_return$confidence, 2)
  
  return(overview_return)
}

#' Little helperfunction to decide if an inversion can be automatically rejected. This is the case if
#' the LLHs for HOM and both HETs are below reference. It wouldn't have been necessary to write its own
#' function for this, but well here it is now.ok then 
#'
#' @param cl call_llhs
#' @param chrom the chromosome of the inversion of interest
#' @param start the start position of the inversion of interest
#' @author Wolfram Hoeps
#' @export
test_for_invrejection <- function(cl, chrom, start){
  # Read the values for inv_hom, inv_h1, inv_h2
  vals_to_check = cl[(cl$chrom==chrom) & (cl$start==start),][c('sum_0101', 'sum_0110', 'sum_1001')]
  # Are all of them < 0? Then we can reject.
  can_be_rejected = all(vals_to_check < 0)
  # True: reject. False: keep
  return(can_be_rejected)
}
