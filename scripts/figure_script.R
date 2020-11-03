
### load R packages####

library(tidyverse)
library(ggtree)
library(cowplot)
library(ggpubr)
library(ggalt)    
library(ggthemes) 
library(ggrepel)

#####

##########################################
#           Figure 1                     #
#      Swept haplotype, tree             #
##########################################

## Figure 1A #####

data_S1A <- read.csv("processed_data/S1A_File.csv", stringsAsFactors=FALSE)

fig1a <-  ggplot(data_S1A,
                 aes(xmin = start/1E6, xmax = stop/1E6,
                     ymin = -(plotpoint_hi + 0.5), ymax = -(plotpoint_hi - 0.5),
                     fill = swept_haplotype)) +
  geom_rect() +
  scale_fill_manual(values = c("plum4","gold2")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Genomic position (Mb)") +
  theme_bw() +
  facet_grid(.~chromosome, scales="free", space="free") +
  theme() +
  theme(axis.text =  element_text(size=12,  color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size=12,  color = "black"),
        legend.title =  element_text(size=12,  color = "black"),
        axis.title =  element_text(size=12,  color = "black"),
        strip.text = element_text(size=12, vjust = 1,  color = "black"),
        strip.background = element_blank(), 
        legend.position = "none",
        text=element_text(family="Helvetica"),
        plot.margin = unit(c(3, 0, 0.5, 0), "mm"))



## Figure 1B #####

load("processed_data/S1B_File.RData")

fig1b <- ggtree(data_S1B, layout="rectangular", branch.length="rate", size = 0.3,color="gray51") +
  geom_tippoint( aes(fill=label,  color=label),  
                 shape=21, 
                 size= 1) +
  scale_y_reverse() + 
  theme(legend.position =  c(0.6,0.7),
        legend.text = element_text(size=12,  color = "black"),
        legend.title =  element_blank(),
        text=element_text(family="Helvetica"),
        plot.margin = unit(c(0, 0, 0, 0), "mm")) + 
  scale_fill_manual(values = c("plum4", "deepskyblue4",   "lightskyblue2", 
                               "darkorange1",  "gold2")) +
  scale_color_manual(values = c("plum4",  "deepskyblue4",   "lightskyblue2", 
                                "darkorange1",  "gold2")) 


# cowplot
fig_1AB <- cowplot::plot_grid(fig1a, fig1b, 
                              labels = c('A', 'B'), 
                              rel_widths = c(6,1),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "lr",
                              nrow = 1)

ggsave(fig_1AB, filename = paste("figures/Fig_1.png",sep = ""), units = "mm",height = 225, width = 170)






##########################################
#           Figure 2                     #
#      Swept haplotype, tree             #
##########################################

## Figure 2A ####

data_S2A <- read.csv("processed_data/S2A_File.csv", stringsAsFactors=FALSE)


fig2a <- ggplot(data_S2A,aes(x=fct_reorder(strain, mean_b),y=mean_b)) + 
  geom_bar(stat='identity',fill="gray80",color="black") + 
  geom_errorbar(aes(ymin=mean_b-se_b, ymax=mean_b+se_b),
                cex = 0.2,width=.2,                   
                position=position_dodge(.6)) + theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12,  color = "black"),
        axis.text.y = element_text(size=12,  color = "black"),
        legend.text = element_text(size=12,  color = "black"),
        legend.title =  element_text(size=12,  color = "black"),
        panel.grid = ggplot2::element_blank(),
        legend.position=c(0.25,0.85)) +
  xlab("Strain") + 
  ylab("Brood size")  + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 370)) 



## Figure 2B total ######

box_swept <- data_S2A %>%
  dplyr::mutate(day =  "Brood size") 

stat_box_data <- function(y, upper_limit = max(box_swept$mean_b) * 1.15) {
  return( 
    data.frame(
      y = 0.98 * upper_limit,
      label = paste(length(y), '\n'
      )
    )
  )
}

fig2b_l <- ggplot(box_swept,aes(x=genotype,y=mean_b)) + 
  geom_jitter(shape=21,position=position_jitter(0.2),size=0.5,fill="black") +
  geom_boxplot(aes(fill=genotype),outlier.shape = NA,alpha=0.5) + 
  scale_fill_manual(values=c("gold2","plum4"),labels = c("swept", "divergent")) + 
  theme_bw() + facet_grid(.~day) +
  theme(axis.text.y =  element_text(size=12,  color = "black"),
        axis.text.x =  element_blank(),
        axis.title.x =  element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=12,  color = "black"),
        legend.title =   element_blank(),
        axis.title =  element_text(size=12,  color = "black"),
        strip.text = element_text(size=12, vjust = 1,  color = "black"),
        strip.background = element_blank(), 
        panel.grid = ggplot2::element_blank(),
        legend.position="bottom",  
        legend.spacing.x = unit(0.1, 'cm'),
        text=element_text(family="Helvetica"))+  
  ylab("Brood size")  + 
  xlab("Genotype")  + 
  scale_y_continuous(limits =c(100,400),breaks=seq(100,400,50)) +
  stat_compare_means(comparisons = list(c("swept","unswept")),label.y = c(350),label = "p.signif")+
  labs(color="Genotype",fill="Genotype")  +
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.4, size = 12*5/14
  ) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) 


## Figure 2B by day ######

data_S2B <- read.csv("processed_data/S2B_File.csv", stringsAsFactors=FALSE)

fig2b_r <- ggplot(data_S2B,aes(x=genotype,y=mean_b)) + 
  facet_grid(.~day) + 
  geom_jitter(shape=21,position=position_jitter(0.2),size=0.5,fill="black") +
  geom_boxplot(aes(fill=genotype),outlier.shape = NA,alpha=0.5) +
  scale_fill_manual(values=c("gold2","plum4")) + 
  theme_bw() +
  theme(axis.text.x =  element_blank(),
        axis.title =  element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y =  element_text(size=12,  color = "black"),
        legend.text = element_text(size=12,  color = "black"),
        legend.title =  element_text(size=12,  color = "black"),
        strip.text = element_text(size=12, vjust = 1,  color = "black"),
        strip.background = element_blank(), 
        legend.position="none",
        panel.grid = ggplot2::element_blank(),
        legend.spacing.x = unit(1, 'mm'),
        legend.background = element_rect(fill=alpha('white', 0)),
        text=element_text(family="Helvetica")) +  
  ylab("Brood size")  + 
  scale_y_continuous(limits =c(0,210),breaks=seq(0,210,50), oob = scales::squish) +  
  stat_compare_means(comparisons = list(c("swept","unswept")),label.y = c(201),label = "p.signif") +
  labs(color="Genotype",fill="Genotype") 


### cow fig2 



fig2b <- cowplot::plot_grid(fig2b_l,fig2b_r,
                            nrow = 1,
                            ncol = 2,
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "tb",
                            rel_widths = c(1.5,5))



fig_2 <- cowplot::plot_grid(fig2a, fig2b,
                            labels = c("A","B"), 
                            nrow = 2,
                            ncol = 1,
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr")

ggsave(fig_2, filename = paste( "figures/Fig_2.png",sep=""), units = "mm",height = 190, width = 170)



##########################################
#           Figure 3                     #
#      manhanton,pxg                     #
##########################################

## Figure 3A manhanton #####

data_S3A <- read.csv("processed_data/S3A_File.csv", stringsAsFactors=FALSE)


fig3a <- ggplot2::ggplot(data_S3A,aes(x = POS/1e6, y = log10p)) +
  ggplot2::scale_color_manual(values = c("0" = "black", 
                                         "1" = "red",
                                         "2" = "hotpink3")) +
  ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6,   
                                  xmax = endPOS/1e6,    
                                  ymin = 0, 
                                  ymax = Inf, 
                                  fill = "plum4"), 
                     color = "plum4",fill = "cyan",linetype = 2, 
                     alpha=.3) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                      color = "gray", 
                      alpha = .75,  
                      size = 1) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                      color = "gray", 
                      alpha = .75,  
                      size = 1,
                      linetype = 2) +
  ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG)),size=0.5 ) +
  ggplot2::facet_grid( . ~ CHROM, scales = "free" , space = "free") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12,  color = "black"),
                 axis.text.y = ggplot2::element_text(size = 12,  color = "black"),
                 axis.title.x = ggplot2::element_text(size = 12,  color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 12,  color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 12, color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 12, color = "black"), 
                 plot.title = ggplot2::element_text(size = 12, vjust = 1), 
                 panel.grid = ggplot2::element_blank(),
                 legend.position = "none",
                 strip.background = element_blank()) +
  ggplot2::labs(x = "Genomic position (Mb)",
                y = expression(-log[10](italic(p))))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 9.5)) 





## Figure 3B pxg#####

data_S3B <- read.csv("processed_data/S3B_File.csv", stringsAsFactors=FALSE)

fig3b <- ggplot(data_S3B, aes(x=factor(as.character(allele),labels = c("REF","ALT")),y = value)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21,
              position=position_jitter(0.2),
              size=1.5,
              aes(fill=chr_geno),
              alpha=0.8) +                      
  scale_fill_manual(values=c("gold2","plum4"),labels = c("swept", "divergent")) + 
  theme_bw() + 
  facet_grid(.~facet_marker) +
  theme(axis.text =  element_text(size=12,  color = "black"),
        plot.title = element_text(hjust = 0.1, size=12),
        legend.text = element_text(size=12,  color = "black"),
        legend.title =  element_blank() ,
        axis.title =  element_text(size=12,  color = "black"),
        axis.title.x = element_blank() ,
        strip.text = element_text(size=12, vjust = 1,  color = "black"),
        strip.background = element_blank(), 
        panel.grid = ggplot2::element_blank(),
        text=element_text(family="Helvetica")) +  
  ylab("Brood size")  +
  xlab("Genotype") +
  labs(fill="Swept ratio") 

### cow fig3 #


fig_3 <- cowplot::plot_grid(fig3a,fig3b,
                            labels = c("A","B"), 
                            nrow = 2,
                            ncol = 1,
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "l")


ggsave(fig_3, filename = paste( "figures/Fig_3.png",sep=""), units = "mm",height = 130, width = 170)




##########################################
#           Figure 4                     #
#      Geographical comparison           #
##########################################


data_S4 <- read.csv("processed_data/S4_File.csv", stringsAsFactors=FALSE)

box_swept_loc <- data_S4 

box_swept_loc$loc = factor(box_swept_loc$Location, levels = c("Hawaii","North America","Europe"))

my_comparisons <- list( c("Hawaii", "North America"), c("Hawaii", "Europe"))

cc<-compare_means(mean_b ~ loc,  data = box_swept_loc)


geo_box <- ggplot(box_swept_loc,aes(x=loc,y=mean_b)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21,
              position=position_jitter(0.2),
              size=2.5,
              aes(fill=genotype),
              alpha=0.8) +
  stat_compare_means(comparisons = my_comparisons, 
                     ref.group = "Hawaii",
                     label = "p.signif",
                     label.y = c(320,350)) +                      
  scale_fill_manual(values=c("gold2" ,"plum4"),labels = c("swept", "divergent")) + 
  theme_bw() +
  theme(axis.text =  element_text(size=12,  color = "black"),
        plot.title = element_text(hjust = 0.1, size=12,  color = "black"),
        legend.text = element_text(size=12,  color = "black"),
        axis.title =  element_text(size=12,  color = "black"),
        strip.text = element_text(size=12, vjust = 1,  color = "black"),
        strip.background = element_blank(), 
        legend.title = element_blank(),
        legend.spacing.x = unit(5, 'mm'),
        panel.grid = ggplot2::element_blank(),
        text=element_text(family="Helvetica"))+  
  ylab("Brood size")  + xlab("Sampling location")+
  labs(color="Sampling location",fill="Swept ratio") +
  ylim(100,370)



ggsave(geo_box, filename = paste( "figures/Fig_4.png",sep=""), units = "mm",height = 90, width =140)



##########################################
#           Figure S1                    #
#      Geographical distribution         #
##########################################

data_SS1 <- read.csv("processed_data/SS1_File.csv", stringsAsFactors=FALSE)


# modify settings
background_color <- "white"
axes_text_size <- "black"
axis_color <- "gray51"
axes_text_size <- 10


# load world map and select Hawaii > Big island
world <- map_data("world") %>% dplyr::filter(!region=="Antarctica")

strain_to_plot <- data_SS1 %>%
  dplyr::mutate(strain2=" ")


plt_wolrd <- ggplot()+ 
  geom_map(data=world, map=world,
           aes(x=long, y=lat, map_id=region),
           color=axis_color, fill=background_color, size=0.5,alpha=0.5) +
  geom_point(data = strain_to_plot, 
             aes(x=as.numeric(longitude), y=as.numeric(latitude)), 
             color = "black",
             shape = 21,   
             size = 0.1,
             alpha=0.1) + 
  geom_label_repel(data=strain_to_plot, 
                   aes(x=as.numeric(longitude), 
                       y=as.numeric(latitude),
                       label = strain2,
                       fill = factor(n_swept_chr)),
                   box.padding = 0.005,
                   label.padding = 0.15, point.padding = 0.1,
                   segment.alpha =0.5,
                   segment.size=0.3,
                   size=0.1) + 
  scale_fill_manual(values = c("plum4", "deepskyblue4",   "lightskyblue2", 
                               "darkorange1",  "gold2")) + 
  theme_map() +
  theme(legend.position ="bottom",
        legend.text = element_text(size=12,  color = "black"),
        legend.title =  element_text(size=12,  color = "black")) +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(fill='Number of swept chromosomes')

ggsave(plt_wolrd, filename = paste("figures/Fig_S1.png",sep = ""), units = "mm",height = 100, width = 170)



##########################################
#           Figure S2                    #
#   broods colored by swept chromosome   #
##########################################


data_SS2 <- read.csv("processed_data/SS2_File.csv", stringsAsFactors=FALSE)


stat_box_data <- function(y, upper_limit = max(data_SS2$mean_b) * 1.15) {
  return( 
    data.frame(
      y = 0.98 * upper_limit,
      label = paste(length(y), '\n'
      )
    )
  )
}

figS2b <- ggplot(data_SS2,aes(x=chr_geno,y=mean_b)) + 
  geom_jitter(shape=21,position=position_jitter(0.2),size=0.5,fill="black") +
  geom_boxplot(aes(fill=chr_geno),outlier.shape = NA,alpha=0.5) + 
  scale_fill_manual(values=c("gold2","plum4"),labels = c("swept", "divergent")) + 
  theme_bw() + facet_grid(.~swept_region) +
  theme(axis.text.y =  element_text(size=12,  color = "black"),
        axis.text.x =  element_blank(),
        axis.title.x =  element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size=12,  color = "black"),
        legend.title =   element_blank(),
        axis.title =  element_text(size=12,  color = "black"),
        strip.text = element_text(size=12, vjust = 1,  color = "black"),
        strip.background = element_blank(), 
        panel.grid = ggplot2::element_blank(),
        legend.position="right",  
        text=element_text(family="Helvetica"))+  
  ylab("Brood size")  + 
  xlab("Genotype")  + 
  scale_y_continuous(limits =c(100,400),breaks=seq(100,400,50)) +
  stat_compare_means(comparisons = list(c("swept","unswept")),label.y = c(350),label = "p.signif")+
  labs(color="Genotype",fill="Genotype")  +
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.4, size = 12*5/14
  ) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


ggsave(figS2b, filename = paste("figures/Fig_S2.png",sep = ""), units = "mm",height = 80, width = 170)




##########################################
#           Figure S3                    #
#       GWAS EIGEN & LD                  #
##########################################


## Figure S3A GWAS EIGEN manhanton #####

data_SS3A <- read.csv("processed_data/SS3A_File.csv", stringsAsFactors=FALSE)

figS3a <- ggplot2::ggplot(data_SS3A,aes(x = POS/1e6, y = log10p)) +
  ggplot2::scale_color_manual(values = c("0" = "black", 
                                         "1" = "red",
                                         "2" = "hotpink3")) +
  ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6,   
                                  xmax = endPOS/1e6,    
                                  ymin = 0, 
                                  ymax = Inf, 
                                  fill = "plum4"), 
                     color = "plum4",fill = "cyan",linetype = 2, 
                     alpha=.3) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                      color = "gray", 
                      alpha = .75,  
                      size = 1) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN_CUTOFF),
                      color = "gray", 
                      alpha = .75,  
                      size = 1,
                      linetype = 2) +
  ggplot2::geom_point( ggplot2::aes(color= factor(EIGEN_SIG)),size=0.5 ) +
  ggplot2::facet_grid( . ~ CHROM, scales = "free" , space = "free") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12,  color = "black"),
                 axis.text.y = ggplot2::element_text(size = 12,  color = "black"),
                 axis.title.x = ggplot2::element_text(size = 12,  color = "black", vjust = -0.3), 
                 axis.title.y = ggplot2::element_text(size = 12,  color = "black"), 
                 strip.text.x = ggplot2::element_text(size = 12, color = "black"), 
                 strip.text.y = ggplot2::element_text(size = 12, color = "black"), 
                 plot.title = ggplot2::element_text(size = 12, vjust = 1), 
                 panel.grid = ggplot2::element_blank(),
                 legend.position = "none",
                 strip.background = element_blank()) +
  ggplot2::labs(x = "Genomic position (Mb)",
                y = expression(-log[10](italic(p))))+ 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 9.5)) 


## Figure S3B GWAS EIGEN QTL LD #####

data_SS3B <- read.csv("processed_data/SS3B_File.csv", stringsAsFactors=FALSE)


figS3b <- ggplot2::ggplot(data_SS3B) +
  ggplot2::aes(x = factor(SNP1, levels = unique(SNP1), ordered = T), y = factor(SNP2, levels = unique(SNP1), ordered = T)) +
  ggplot2::geom_tile(ggplot2::aes(fill = r2)) +
  ggplot2::geom_text(ggplot2::aes(label = signif(r2, 3)),  size = 12*5/14) +
  ggplot2::theme(text = ggplot2::element_text(size = 12, color = "black"),
                 axis.text.x = ggplot2::element_text(size = 8, color = "black"),
                 axis.text.y = ggplot2::element_text(size = 12, color = "black"),
                 axis.title.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 legend.text = element_text(size=12,  color = "black"),
                 legend.title =   element_text(size=12,  color = "black"),
                 legend.position = c(0.1,0.6)
  ) +
  scale_x_discrete(labels = function(x) { gsub("_", ":", x) }, expand = c(0, 0)) + 
  scale_y_discrete(position = "right", labels = function(x) { gsub("_", ":", x) }, expand = c(0, 0)) + 
  scale_fill_continuous(high = "red", low = "white", na.value = "white") +
  labs(fill=expression(R^{2}))




### cow figS3 #


fig_S3 <- cowplot::plot_grid(figS3a,figS3b,
                             labels = c("A","B"), 
                             nrow = 2,
                             ncol = 1,
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             axis = "l")


ggsave(fig_S3, filename = paste( "figures/Fig_S3.png",sep=""), units = "mm",height = 130, width = 175)

