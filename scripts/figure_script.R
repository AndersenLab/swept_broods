
base_dir <- "~/swept_broods/"

setwd(base_dir)


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

data_S1A <- read.csv("processed_data/FileS1_haplotypes.csv", stringsAsFactors=FALSE)

fig1a <-  ggplot(data_S1A,
                 aes(xmin = start/1E6, xmax = stop/1E6,
                     ymin = -(plotpoint_hi + 0.5), ymax = -(plotpoint_hi - 0.5),
                     fill = swept_haplotype)) +
  geom_rect() +
  scale_fill_manual(values = c("gray","red")) +
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

load("processed_data/FileS2_tree.RData")

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

data_S2A <- read.csv("processed_data/FileS3_lifetimeFertility.csv", stringsAsFactors=FALSE)

bar_S2A <- data_S2A %>% dplyr::mutate(Strain=ifelse(strain %in% c("N2","CB4856"),strain,"Other strains"))


fig2a <- ggplot(bar_S2A,aes(x=fct_reorder(strain, mean_b),y=mean_b,fill=Strain)) + 
  geom_bar(stat='identity',color="gray51",size=.3) + 
  scale_fill_manual(values=c("N2"="orange","CB4856"="blue","Other strains"="gray88")) + 
  geom_errorbar(aes(ymin=mean_b-se_b, ymax=mean_b+se_b),
                cex = 0.2,width=.2,                   
                position=position_dodge(.6)) + theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12,  color = "black"),
        axis.text.y = element_text(size=12,  color = "black"),
        legend.text = element_text(size=12,  color = "black"),
        legend.title =  element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.position=c(0.3,0.9)) +
  xlab("Strain") + 
  ylab("Lifetime fertility")  + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 370)) + guides(fill = guide_legend(nrow = 1))


## Figure 2B total ######

box_swept <- data_S2A %>%
  dplyr::mutate(day =  "Lifetime") 


box_swept$Genotypes <- factor(box_swept$genotype,levels = c("swept","divergent"))


stat_box_data <- function(y, upper_limit = max(box_swept$mean_b) * 1.15) {
  return( 
    data.frame(
      y = 0.98 * upper_limit,
      label = paste(length(y), '\n'
      )
    )
  )
}

fig2b_l <- ggplot(box_swept,aes(x=Genotypes,y=mean_b)) + 
  geom_jitter(shape=21,position=position_jitter(0.2),size=0.5,fill="black") +
  geom_boxplot(aes(fill=Genotypes),outlier.shape = NA,alpha=0.5) + 
  scale_fill_manual(values=c("swept"="gold2","divergent"="plum4")) + 
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
  ylab("Fertility")  + 
  xlab("Genotype")  + 
  scale_y_continuous(limits =c(100,400),breaks=seq(100,400,50)) +
  stat_compare_means(comparisons = list(c("swept","divergent")),label.y = c(350),label = "p.signif")+
  labs(color="Genotype",fill="Genotype")  +
  stat_summary(
    fun.data = stat_box_data, 
    geom = "text", 
    hjust = 0.5,
    vjust = 0.4, size = 12*5/14
  ) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


## Figure 2B by day ######

data_S2B <- read.csv("processed_data/FileS4_dailyFertility.csv", stringsAsFactors=FALSE)

data_S2B$Genotypes <- factor(data_S2B$genotype,levels = c("swept","divergent"))

data_S2B$day[data_S2B$day=="Day 5"]  <- "Days 5-7" 


fig2b_r <- ggplot(data_S2B,aes(x=Genotypes,y=mean_b)) + 
  facet_grid(.~day) + 
  geom_jitter(shape=21,position=position_jitter(0.2),size=0.5,fill="black") +
  geom_boxplot(aes(fill=Genotypes),outlier.shape = NA,alpha=0.5) +
  scale_fill_manual(values=c("swept"="gold2","divergent"="plum4")) + 
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
  stat_compare_means(comparisons = list(c("swept","divergent")),label.y = c(201),label = "p.signif") +
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
                            axis = "lr",
                            rel_heights =  c(1,2))

ggsave(fig_2, filename = paste( "figures/Fig_2.png",sep=""), units = "mm",height = 140, width = 170)



##########################################
#           Figure 3                     #
#      manhanton,pxg                     #
##########################################

## Figure 3A manhanton #####

data_S3A <- read.csv("processed_data/FileS5_GWA121.csv", stringsAsFactors=FALSE)

fig3a <- manhaplot(data_S3A) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 10.4)) 
  


## Figure 3B pxg#####

data_S3B <- read.csv("processed_data/FileS6_pxg121.csv", stringsAsFactors=FALSE)

fig3b <- pxg_plot(data_S3B)

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


data_S4 <- read.csv("processed_data/FileS7_geo.csv", stringsAsFactors=FALSE)

data_S4$Genotypes <- factor(data_S4$genotype,levels = c("swept","divergent"))


box_swept_loc <- data_S4 

box_swept_loc$loc = factor(box_swept_loc$Location, levels = c("Hawaii","North America","Europe"))

my_comparisons <- list( c("Hawaii", "North America"), c("Hawaii", "Europe"))

cc<-compare_means(mean_b ~ loc,  data = box_swept_loc)


geo_box <- ggplot(box_swept_loc,aes(x=loc,y=mean_b)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21,
              position=position_jitter(0.2),
              size=2.5,
              aes(fill=Genotypes),
              alpha=0.8) +
  stat_compare_means(comparisons = my_comparisons, 
                     ref.group = "Hawaii",
                     label = "p.signif",
                     label.y = c(320,350)) +                      
  scale_fill_manual(values=c("swept"="gold2","divergent"="plum4")) + 
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
  ylab("Lifetime fertility")  + xlab("Sampling location")+
  labs(color="Sampling location",fill="Swept ratio") +
  ylim(100,370)



ggsave(geo_box, filename = paste( "figures/Fig_4.png",sep=""), units = "mm",height = 90, width =140)



##########################################
#           Figure S1                    #
#      Geographical distribution         #
##########################################

data_SS1 <- read.csv("processed_data/FileS8_distribution.csv", stringsAsFactors=FALSE)


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


data_SS2 <- read.csv("processed_data/FileS9_swpchr.csv", stringsAsFactors=FALSE)

data_SS2$Genotypes <- factor(data_SS2$chr_geno,levels = c("swept","divergent"))

stat_box_data <- function(y, upper_limit = max(data_SS2$mean_b) * 1.15) {
  return( 
    data.frame(
      y = 0.98 * upper_limit,
      label = paste(length(y), '\n'
      )
    )
  )
}

figS2b <- ggplot(data_SS2,aes(x=Genotypes,y=mean_b)) + 
  geom_jitter(shape=21,position=position_jitter(0.2),size=0.5,fill="black") +
  geom_boxplot(aes(fill=Genotypes),outlier.shape = NA,alpha=0.5) + 
  scale_fill_manual(values=c("swept"="gold2","divergent"="plum4")) + 
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
  ylab("Lifetime fertility")  + 
  xlab("Genotype")  + 
  scale_y_continuous(limits =c(100,400),breaks=seq(100,400,50)) +
  stat_compare_means(comparisons = list(c("swept","divergent")),label.y = c(350),label = "p.signif")+
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
#               LD                       #
##########################################


## Figure S3B GWAS EIGEN QTL LD #####

data_SS3 <- read.csv("processed_data/FileS10_LD121.csv", stringsAsFactors=FALSE)


figS3 <- ggplot2::ggplot(data_SS3) +
  ggplot2::aes(x = factor(SNP1, levels = unique(SNP1), ordered = T), y = factor(SNP2, levels = unique(SNP1), ordered = T)) +
  ggplot2::geom_tile(ggplot2::aes(fill = r2)) +
  ggplot2::geom_text(ggplot2::aes(label = signif(r2, 3)),  size = 12*5/14) +
  ggplot2::theme(text = ggplot2::element_text(size = 12, color = "black"),
                 axis.text.x = ggplot2::element_text(size = 12, color = "black"),
                 axis.text.y = ggplot2::element_text(size = 12, color = "black"),
                 axis.title.x = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 legend.text = element_text(size=12,  color = "black"),
                 legend.title =   element_text(size=12,  color = "black"),
                 legend.position = c(0.1,0.7)
  ) +
  scale_x_discrete(labels = function(x) { gsub("_", ":", x) }, expand = c(0, 0)) + 
  scale_y_discrete(position = "right", labels = function(x) { gsub("_", ":", x) }, expand = c(0, 0)) + 
  scale_fill_continuous(high = "red", low = "white", na.value = "white") +
  labs(fill=expression(r^{2}))


  
ggsave(figS3, filename = paste( "figures/Fig_S3.png",sep=""), units = "mm",height = 80, width = 100)


##########################################
#           Figure S4                    #
#      QTL interval haplotype            #
##########################################


## Figure S4 QTL interval haplotype #####

SS4_File <- read.csv("processed_data/FileS11_QTLhaplotype121.csv", stringsAsFactors=FALSE)


####ref_swept ##
figSS4_ref <- iv_hap(subset(SS4_File,ALL=="REF")) +
  ylab("REF") +
  facet_wrap_custom(.~facet_marker, scales="free",nrow = 1,
                    scale_overrides = list(
                      scale_override(1, scale_x_continuous(breaks = c(13512687/1E6, 13917228/1E6, 14195881/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(2, scale_x_continuous(breaks = c(4171/1E6, 543326/1E6, 758213/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(3, scale_x_continuous(breaks = c(7887491/1E6, 14534671/1E6, 16016868/1E6),labels = scales::number_format(accuracy = 0.1)))
                    ))


####alt_swept ##
figSS4_alt <- iv_hap(subset(SS4_File,ALL=="ALT")) +
  ylab("ALT")  +
  facet_wrap_custom(.~facet_marker, scales="free",nrow = 1,
                    scale_overrides = list(
                      scale_override(1, scale_x_continuous(breaks = c(13512687/1E6, 13917228/1E6, 14195881/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(2, scale_x_continuous(breaks = c(4171/1E6, 543326/1E6, 758213/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(3, scale_x_continuous(breaks = c(7887491/1E6, 14534671/1E6, 16016868/1E6),labels = scales::number_format(accuracy = 0.1)))
                    )) +
  theme(axis.text =  element_text(size=12,  color = "black"),
        strip.text = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=12,  color = "black"),
        legend.title = element_blank(),
        axis.title.x =  element_text(size=12, vjust = 1,  color = "black"))


#### cow s4 ##

figS4_legend <- get_legend(figSS4_alt)


figS4_main  <- cowplot::plot_grid(figSS4_ref ,
                                  figSS4_alt + theme(legend.position="none"),
                                  nrow = 2,ncol = 1,
                                  axis = "lr")



figS4 <- cowplot::plot_grid(figS4_main ,figS4_legend,
                            nrow = 2,ncol = 1,
                            rel_heights = c(10,1.5))



ggsave(figS4, filename = paste( "figures/Fig_S4.png",sep=""), units = "mm",height = 150, width = 170)





##########################################
#           Figure S5                    #
#          norm.n GWA dmso                #
##########################################


## Figure S5 A manhanton #####

data_SS5A <- read.csv("processed_data/FileS12_GWA236.csv", stringsAsFactors=FALSE)

figS5_A <- manhaplot(data_SS5A) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 


## Figure S5 B pxg #####


data_SS5B <- read.csv("processed_data/FileS13_pxg236.csv", stringsAsFactors=FALSE)


figS5_B <- pxg_plot(data_SS5B) +
  theme(legend.position = "none")  +  
  ylab("Fertility") 


## Figure S5 C interval haplotype #####

data_SS5C <- read.csv("processed_data/FileS14_QTLhaplotype236.csv", stringsAsFactors=FALSE)

figS5C_ref <- iv_hap(subset(data_SS5C,ALL=="REF")) +
  ylab("REF") +
  theme( plot.margin = unit(c(0, 2, 0, 10), "mm"))+
  facet_wrap_custom(.~facet_marker, scales="free",nrow = 1, 
                    scale_overrides = list(
                     # scale_override(1, scale_x_continuous(breaks = c(9529207/1E6, 10256873/1E6, 10691570/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(1, scale_x_continuous(breaks = c(3852481/1E6, 4831537/1E6, 5424864/1E6),labels = scales::number_format(accuracy = 0.1)))
                    ))

figS5C_alt <- iv_hap(subset(data_SS5C,ALL=="ALT")) +
  ylab("ALT")  +
  theme( strip.text = element_blank(),
         axis.text =  element_text(size=12,  color = "black"),
         plot.margin = unit(c(0, 2, 0, 10), "mm"),
         axis.title.x =  element_text(size=12, vjust = 1,  color = "black"))+
  facet_wrap_custom(.~facet_marker, scales="free",nrow = 1, 
                    scale_overrides = list(
                   #   scale_override(1, scale_x_continuous(breaks = c(9529207/1E6, 10256873/1E6, 10691570/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(1, scale_x_continuous(breaks = c(3852481/1E6, 4831537/1E6, 5424864/1E6),labels = scales::number_format(accuracy = 0.1)))
                    ))


#### cow s5 ##

figS5_C <- cowplot::plot_grid(figS5C_ref ,figS5C_alt,
                              nrow = 2,ncol = 1,
                              axis = "lr",
                              rel_heights = c(4,1))

figS5_BC <- cowplot::plot_grid(figS5_B,figS5_C,
                                nrow = 1,ncol = 2,
                                labels = c("B","C"), 
                                label_size = 12, 
                                label_fontfamily="Helvetica",
                                axis = "lrtb")



figS5 <- cowplot::plot_grid(figS5_A ,figS5_BC,
                            nrow = 2,ncol = 1,
                            labels = c("A",""), 
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr",
                            rel_heights = c(1,1.1))




ggsave(figS5, filename = paste( "figures/Fig_S5.png",sep=""), units = "mm",height = 120, width = 170)





##########################################
#           Figure S6                    #
#      norm.n linkage-mapping 1% H2O     #
##########################################


load("processed_data/FileS15_lk1h2o.RData")

figS6_A <- lk_lod_plot(lkmap_1h2o,cis_1h2o)



data_SS6B <- read.csv("processed_data/FileS16_pxg1h2o.csv", stringsAsFactors=FALSE)


data_SS6B$peakmarker <- factor(data_SS6B$marker, levels = c(" Parental", "II:3646789"))


figS6_B <- lk_pxg_plot(data_SS6B) + 
  ggplot2::facet_wrap(~peakmarker, ncol = 2, scales = "free_x") + 
  ggplot2::labs(x = "", y = paste("Fertility","(1% H2O)",sep = " "))



figS6 <- cowplot::plot_grid(figS6_A,figS6_B,
                            labels = c("A","B"), 
                            nrow = 2,
                            ncol = 1,
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr")


ggsave(figS6, filename = paste( "figures/Fig_S6.png",sep=""), units = "mm",height = 120, width = 170)




##########################################
#           Figure S7                    #
#      norm.n linkage-mapping 1% DMSO    #
##########################################


load("processed_data/FileS17_lk1dmso.RData")

figS7_A <- lk_lod_plot(lkmap_1dmso,cis_1dmso)



data_SS7B <- read.csv("processed_data/FileS18_pxg1dmso.csv", stringsAsFactors=FALSE)


data_SS7B$peakmarker <- factor(data_SS7B$marker, levels = c(" Parental",  "IV:9031007" ,"V:11883496"  ))



figS7_B <- lk_pxg_plot(data_SS7B) + 
  ggplot2::facet_wrap(~peakmarker, ncol = 3, scales = "free_x") + 
  ggplot2::labs(x = "", y = paste("Fertility","(1% DMSO)",sep = " "))



figS7 <- cowplot::plot_grid(figS7_A,figS7_B,
                            labels = c("A","B"), 
                            nrow = 2,
                            ncol = 1,
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr")


ggsave(figS7, filename = paste( "figures/Fig_S7.png",sep=""), units = "mm",height = 120, width = 170)








##########################################
#           Figure S8                    #
#      norm.n linkage-mapping 0.5% DMSO  #
##########################################


load("processed_data/FileS19_lk05dmso.RData")

figS8_A <- lk_lod_plot(lkmap_05dmso,cis_05dmso)



data_SS8B <- read.csv("processed_data/FileS20_pxg05dmso.csv", stringsAsFactors=FALSE)


data_SS8B$peakmarker <- factor(data_SS8B$marker, levels = c(" Parental", "II:8964146" , "IV:8044494" , "IV:16522256" ,"V:11883496"  ))



figS8_B <- lk_pxg_plot(data_SS8B) + 
  ggplot2::facet_wrap(~peakmarker, ncol = 5, scales = "free_x") + 
  ggplot2::labs(x = "", y = paste("Fertility","(0.5% DMSO)",sep = " "))



figS8 <- cowplot::plot_grid(figS8_A,figS8_B,
                            labels = c("A","B"), 
                            nrow = 2,
                            ncol = 1,
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr")


ggsave(figS8, filename = paste( "figures/Fig_S8.png",sep=""), units = "mm",height = 120, width = 170)






##########################################
#           Figure 5                     #
#       Mapping summary                  #
##########################################



newrows_size <- data.frame(trait=NA,
                           CHROM=c("I","II","III","IV","V","X"),
                           POS=1,
                           startPOS=0,
                           endPOS = c(14972282, 15173999, 13829314, 17450860, 20914693, 17748731),
                           lod=NA)


fertility_peak_sum_gwa <- bind_rows(data_S3A,data_SS5A) %>% 
  dplyr::distinct(CHROM,POS,startPOS,endPOS,log10p) %>% 
  na.omit() %>% 
  dplyr::mutate(trait=ifelse(CHROM=="X","1% DMSO","Agar plate"))  %>% 
  dplyr::mutate(type="GWA mapping")


gwa_sum <- ggplot(fertility_peak_sum_gwa)+
  aes(x=POS/1E6, y=trait) +
  geom_segment(aes(x = startPOS/1e6, y = trait, xend = endPOS/1e6, yend = trait, color = log10p), size = 2, alpha = 1) +
  geom_segment(data=newrows_size,aes(x = 0, y = trait, xend = endPOS/1e6, yend = trait), size = 2.5, alpha = 0) +
  geom_point(aes(fill=log10p),colour = "black",size = 2, alpha = 1, shape = 25) +
  facet_grid(type ~ CHROM, scales = "free", space = "free") +
  scale_fill_gradient(high = "#D7263D", low = "#0072B2",
                      name = expression(-log[10](italic(p))))+
  scale_color_gradient(high = "#D7263D", low = "#0072B2",
                       name = expression(-log[10](italic(p)))) +
  theme_bw(15) +
  theme(axis.text.y =  element_text(size=12,  color = "black"),
        axis.text.x =  element_blank(), 
        axis.ticks.x =  element_blank(), 
        plot.title = element_text(hjust = 0.1, size=12,  color = "black"),
        legend.text = element_text(size=12,  color = "black"),
        axis.title =  element_blank(), 
        strip.text = element_text(size=12, vjust = 1,  color = "black"),
        strip.background = element_blank(), 
        legend.title = element_text(size=12,  color = "black"),
        legend.spacing.x = unit(5, 'mm'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.grid = ggplot2::element_blank(),
         panel.spacing.x=unit(0.05, "lines"),
         legend.margin=margin(0,0,0,-10),
         text=element_text(family="Helvetica")) 


fertility_peak_sum_lk <- bind_rows(cis_1h2o,cis_05dmso,cis_1dmso) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(trait=ifelse(trait=="X1percwater.norm.n","1% water",
                             ifelse(trait=="halfpercDMSO.norm.n","0.5% DMSO","1% DMSO"))) %>% 
  dplyr::select(trait,CHROM=chr,POS=pos,startPOS=ci_l_pos,endPOS=ci_r_pos,lod)  %>% 
  dplyr::mutate(type="Linkage mapping")




lkm_sum <- ggplot(fertility_peak_sum_lk)+
  aes(x=POS/1E6, y=trait) +
  geom_segment(aes(x = startPOS/1e6, y = trait, xend = endPOS/1e6, yend = trait, color = lod), size = 2, alpha = 1) +
  geom_segment(data=newrows_size,aes(x = 0, y = trait, xend = endPOS/1e6, yend = trait), size = 2.5, alpha = 0) +
  geom_point(aes(fill=lod),colour = "black",size = 2, alpha = 1, shape = 25) +
  facet_grid(type ~ CHROM, scales = "free", space = "free") +
  scale_fill_gradient(high = "#D7263D", low = "#0072B2",
                      name = "LOD") +
  scale_color_gradient(high = "#D7263D", low = "#0072B2",
                       name = "LOD") +
  xlab("Genomic position (Mb)") + 
  ylab("") +
  theme_bw(15)  +
  theme(axis.text.y =  element_text(size=12,  color = "black"),
        axis.text.x =  element_blank(),
        plot.title = element_text(hjust = 0.1, size=12,  color = "black"),
        legend.text = element_text(size=12,  color = "black"),
        axis.title =  element_text(size=12,  color = "black"),
        strip.text.x = element_blank(),
        strip.text = element_text(size=12, vjust = 1,  color = "black"),
        strip.background = element_blank(), 
        legend.title = element_text(size=12,  color = "black"),
        legend.spacing.x = unit(5, 'mm'),
        panel.spacing.x=unit(0.05, "lines"),
        panel.grid = ggplot2::element_blank(),
        legend.margin=margin(0,10,0,-10),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        text=element_text(family="Helvetica")) 




fig5 <- cowplot::plot_grid(gwa_sum,
                           lkm_sum,
                         #  labels = c('A', 'B'), 
                           label_size = 10, 
                           label_fontfamily="Helvetica",
                           rel_heights = c(0.93,1),
                           axis = "lr",
                           align = "v",
                           nrow = 2)


ggsave(fig5, filename = paste( "figures/Fig_5.png",sep=""), units = "mm",height = 100, width =170)

