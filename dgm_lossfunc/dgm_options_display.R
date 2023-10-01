require(ggplot2)
require(ggthemes)
require(data.table)
require(dplyr)
require(interp)
require(gridExtra)
require(ggpubr)


setwd("D:/Code/dgm_lossfunc")

eucall_data = fread("dgm_lossfunctional_EuCallss.csv")


unique(eucall_data$Sample)


get_diff_abs_display <- function(data, interval=1.0) {
  
  #establish the min and max of scale 
  grandmin <- 0.0
  grandmax <- ceiling(max(data$diff_abs))
  #define the number of breaks.  In this case 8 +1 
  mybreaks <- seq(grandmin, grandmax, length.out = grandmax/interval)
  #Function to return the desired number of colors
  
  mycolors<- function(x) {
    colors<-colorRampPalette(c("darkblue","dodgerblue","green","yellow","orange","darkred" ))( grandmax/interval )
    colors[1:x]
  }
  
  #Function to create labels for legend
  breaklabel <- function(x, breaks = mybreaks){
    labels<- paste0(sapply( breaks[1:(length(breaks)-1)], function(x) {  as.character(round(x,1))   }), "-", 
                    sapply( breaks[2:length(breaks)], function(x) {  as.character(round(x,1))   }) )
    labels[1:x]
  }
  
  return(list(grandmin, grandmax, mybreaks, mycolors, breaklabel))
}

diff_abs_list = get_diff_abs_display(eucall_data)
diff_abs_breaks = diff_abs_list[[3]]
diff_abs_colorf = diff_abs_list[[4]]
diff_abs_lblsf = diff_abs_list[[5]]

# ggplot(sample_n(eucall_data[(Exact>0.0) & (Sample=="SL_U")],1000), aes(x=Time, y=S1, z=perc_abs)) + geom_contour()
# ggplot(data, aes(x, y, z = z1)) +
#   geom_contour_filled(breaks= mybreaks, show.legend = TRUE) +
#   scale_fill_manual(palette=mycolors, values=breaklabel(8), name="Value", drop=FALSE) +
#   theme(legend.position = "right")

contour_log10_perc_abs <- function(d1, title_str, palette = "Spectral") {
  d1$perc_abs = log10(d1$perc_abs)
  grid <- with(d1, interp::interp(Time, S1, perc_abs))
  griddf <- subset(data.frame(Time = rep(grid$x, nrow(grid$z)),
                              S1 = rep(grid$y, each = ncol(grid$z)),
                              perc_abs = as.numeric(grid$z)),
                   !is.na(perc_abs))
  p<- ggplot(griddf, aes(Time, S1, z = perc_abs)) +
    geom_contour_filled(colour = "white", show.legend = T) + 
    scale_fill_brewer(palette = palette, direction = -1) + 
    theme_fivethirtyeight() + ggtitle(title_str)
  
  return(p)
}

# d1 = eucall_data[(Exact>0.0) & (Sample=="SL_U")]
# title_str = "European Call, ABS difference - Adam on L2 Loss, Uniform Sampling"
# palette=diff_abs_colorf(length(diff_abs_breaks))
# breaks=diff_abs_breaks
# breaklabels =  diff_abs_lblsf(length(diff_abs_breaks)-1)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

contour_diff_abs <- function(d1, title_str, breaks, palette, breaklabels, show.legend = FALSE) {
  
  grid <- with(d1, interp::interp(Time, S1, diff_abs))
  griddf <- subset(data.frame(Time = rep(grid$x, nrow(grid$z)),
                              S1 = rep(grid$y, each = ncol(grid$z)),
                              diff_abs = as.numeric(grid$z)),
                   !is.na(diff_abs))
  
  p<- ggplot(griddf, aes(Time, S1, z = diff_abs)) +
    geom_contour_filled(colour = "white", breaks=breaks, show.legend = show.legend) + 
    scale_fill_manual(palette=palette, values=breaklabels, name="diff_abs", drop=FALSE) + 
    # scale_fill_brewer(palette = palette, direction = -1) + 
    theme_fivethirtyeight() + ggtitle(title_str)
  return(p)
}

  
SL_U_Acc_perc_abs = contour_log10_perc_abs(eucall_data[(Exact>0.0) & (Sample=="SL_U_NST")], "European Call, log(10) percentage difference - Adam on L2 Loss, Uniform Sampling")
SL_LN_Acc_perc_abs = contour_log10_perc_abs(eucall_data[(Exact>0.0) & (Sample=="SL_LN_NST")], "European Call, log(10) percentage difference - Adam on L2 Loss, LN Sampling")
KL_U_Acc_perc_abs = contour_log10_perc_abs(eucall_data[(Exact>0.0) & (Sample=="KL_U_NST")], "European Call, log(10) percentage difference - Adam on KL(L2) Loss, Uniform Sampling")
KL_LN_Acc_perc_abs = contour_log10_perc_abs(eucall_data[(Exact>0.0) & (Sample=="KL_LN_NST")], "European Call, log(10) percentage difference - Adam on KL(L2) Loss, LN Sampling")


SL_U_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="SL_U_NST")], 
                                     "Adam on L2 Loss, U Sampling", 
                                     palette=diff_abs_colorf,
                                     breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1), show.legend=F)
SL_LN_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="SL_LN_NST")], "Adam on L2 Loss, LN Sampling",
                                      palette=diff_abs_colorf,
                                      breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1), show.legend=F)
KL_U_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="KL_U_NST")], "Adam on KL(L2) Loss, U Sampling",
                                     palette=diff_abs_colorf,
                                     breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1), show.legend=F)
KL_LN_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="KL_LN_NST")], "Adam on KL(L2) Loss, LN Sampling",
                                      palette=diff_abs_colorf,
                                      breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1), show.legend=F)
HL_U_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="HL_U_NST")], 
                                     "Adam on H Loss, U Sampling", 
                                     palette=diff_abs_colorf,
                                     breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1), show.legend=F)
HL_LN_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="HL_LN_NST")], "Adam on H Loss, LN Sampling",
                                      palette=diff_abs_colorf,
                                      breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1),
                                      show.legend=F)

KLH_U_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="KLH_U_NST")], "Adam on KL(H) Loss, U Sampling",
                                      palette=diff_abs_colorf,
                                      breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1),
                                      show.legend = F
)
KLH_LN_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="KLH_LN_NST")], "Adam on KL(H) Loss, LN Sampling",
                                       palette=diff_abs_colorf,
                                       breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1),
                                       show.legend = F)

SL_U_ST_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="SL_U_ST")], "Adam on SL Loss, U ST Sampling",
                                        palette=diff_abs_colorf,
                                        breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1),
                                        show.legend = F)

KL_U_ST_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="KL_U_ST")], "Adam on KL Loss, U ST Sampling",
                                        palette=diff_abs_colorf,
                                        breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1),
                                        show.legend = F)

SLH_U_ST_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="SLH_U_ST")], "Adam on H Loss, U ST Sampling",
                                         palette=diff_abs_colorf,
                                         breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1),
                                         show.legend = F)

KLH_U_ST_Acc_diff_abs = contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="KLH_U_ST")], "Adam on KL(H) Loss, U ST Sampling",
                                         palette=diff_abs_colorf,
                                         breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1),
                                         show.legend = F)


diff_abs_legend = get_legend(contour_diff_abs(eucall_data[(Exact>0.0) & (Sample=="SL_U_NST")], 
                            "European Call, ABS difference - Adam on L2 Loss, U Sampling", 
                            palette=diff_abs_colorf,
                            breaks=diff_abs_breaks, breaklabels = diff_abs_lblsf(length(diff_abs_breaks)-1), show.legend=T))

blankPlot <- ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing()

ggarrange(
  SL_U_Acc_diff_abs, SL_LN_Acc_diff_abs, 
  common.legend = TRUE, legend = "bottom"
)

ggarrange(
  SL_U_Acc_diff_abs, KL_U_Acc_diff_abs, 
  common.legend = TRUE, legend = "bottom"
)

ggarrange(
  ggarrange(SL_U_Acc_diff_abs, SL_LN_Acc_diff_abs, ncol = 2),                # First row with line plot
  # Second row with box and dot plots
  ggarrange(KL_U_Acc_diff_abs, KL_LN_Acc_diff_abs, ncol = 2), 
  nrow = 2, 
  common.legend = TRUE
) 

ggarrange(
  ggarrange(SL_U_Acc_diff_abs, SL_LN_Acc_diff_abs, ncol = 2),                # First row with line plot
  # Second row with box and dot plots
  ggarrange(HL_U_Acc_diff_abs, HL_LN_Acc_diff_abs, ncol = 2), 
  nrow = 2, 
  common.legend = TRUE
) 

ggarrange(
  ggarrange(HL_U_Acc_diff_abs, HL_LN_Acc_diff_abs, ncol = 2), 
  ggarrange(KLH_U_Acc_diff_abs, KLH_LN_Acc_diff_abs, ncol = 2),                # First row with line plot
  nrow = 2, 
  common.legend = FALSE
) 

ggarrange(
  ggarrange(SL_U_Acc_diff_abs, SL_U_ST_Acc_diff_abs, KL_U_ST_Acc_diff_abs, ncol = 3), 
  ggarrange(HL_U_Acc_diff_abs, SLH_U_ST_Acc_diff_abs , KLH_U_ST_Acc_diff_abs, ncol = 3),                # First row with line plot
  nrow = 2, 
  common.legend = FALSE
) 


SL_U_ST_Acc_diff_abs
KL_U_ST_Acc_diff_abs
SLH_U_ST_Acc_diff_abs 
KLH_U_ST_Acc_diff_abs 


grid.arrange(diff_abs_legend, blankPlot,  SL_U_Acc_diff_abs, SL_LN_Acc_diff_abs,
             ncol=2, nrow = 2, 
             widths = c(2.7, 2.7), heights = c(0.5, 2.5))


grid.arrange(SL_U_Acc_diff_abs, SL_LN_Acc_diff_abs, diff_abs_legend, ncol=3, widths=c(2.3, 2.3, 0.4))

grid.arrange(SL_U_Acc_diff_abs,SL_LN_Acc_diff_abs,KL_U_Acc_diff_abs,KL_LN_Acc_diff_abs, nrow=1, ncol=4)

grid.arrange(SL_U_Acc_diff_abs,SL_LN_Acc_diff_abs,KL_U_Acc_diff_abs, nrow=1, ncol=3)
grid.arrange(SL_U_Acc_diff_abs,SL_LN_Acc_diff_abs,nrow=1,ncol=2)


# +  geom_point(data = d1)




