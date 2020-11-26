
base_dir <- "swept_broods/"

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

data_S1A <- read.csv("processed_data/S1A_File.csv", stringsAsFactors=FALSE)

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

bar_S2A <- data_S2A %>% dplyr::mutate(Strain=ifelse(strain %in% c("N2","CB4856"),strain,"Other strains"))


fig2a <- ggplot(bar_S2A,aes(x=fct_reorder(strain, mean_b),y=mean_b,fill=Strain)) + 
  geom_bar(stat='identity',color="black") + 
  scale_fill_manual(values=c("N2"="orange","CB4856"="blue","Other strains"="gray80")) + 
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
  ylab("Brood size")  + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 370)) + guides(fill = guide_legend(nrow = 1))


## Figure 2B total ######

box_swept <- data_S2A %>%
  dplyr::mutate(day =  "Brood size") 


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
  ylab("Number of progeny")  + 
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

data_S2B <- read.csv("processed_data/S2B_File.csv", stringsAsFactors=FALSE)

data_S2B$Genotypes <- factor(data_S2B$genotype,levels = c("swept","divergent"))


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

data_S3A <- read.csv("processed_data/S3A_File.csv", stringsAsFactors=FALSE)

fig3a <- manhaplot(data_S3A)+ 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 10.4)) 
  


## Figure 3B pxg#####

data_S3B <- read.csv("processed_data/S3B_File.csv", stringsAsFactors=FALSE)

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


data_S4 <- read.csv("processed_data/S4_File.csv", stringsAsFactors=FALSE)

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
  ylab("Brood size")  + 
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

data_SS3 <- read.csv("processed_data/SS3_File.csv", stringsAsFactors=FALSE)


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

SS4_File <- read.csv("~/OneDrive - Northwestern University/projects/manuscripts/swept_broods121/processed_data/SS4_File.csv", stringsAsFactors=FALSE)


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
#      swept strans GWA                  #
##########################################


## Figure S5 A manhanton #####

data_SS5A <- read.csv("processed_data/SS5A_File.csv", stringsAsFactors=FALSE)

figS5_A <- manhaplot(data_SS5A) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6.25)) 


## Figure S5 B pxg #####


data_SS5B <- read.csv("processed_data/SS5B_File.csv", stringsAsFactors=FALSE)


figS5_B <- pxg_plot(data_SS5B) +
  theme(legend.position = "none")


## Figure S5 C interval haplotype #####

data_SS5C <- read.csv("processed_data/SS5C_File.csv", stringsAsFactors=FALSE)

figS5C_ref <- iv_hap(subset(data_SS5C,ALL=="REF")) +
  ylab("REF") +
  theme( plot.margin = unit(c(0, 2, 0, 10), "mm"))+
  facet_wrap_custom(.~facet_marker, scales="free",nrow = 1, 
                    scale_overrides = list(
                      scale_override(1, scale_x_continuous(breaks = c(2184029/1E6, 2622125/1E6, 3280545/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(2, scale_x_continuous(breaks = c(2876548/1E6, 3107587/1E6, 3717385/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(3, scale_x_continuous(breaks = c(15288009/1E6, 15883817/1E6, 16167354/1E6),labels = scales::number_format(accuracy = 0.1)))
                    ))

figS5C_alt <- iv_hap(subset(data_SS5C,ALL=="ALT")) +
  ylab("ALT")  +
  theme( strip.text = element_blank(),
         axis.text =  element_text(size=12,  color = "black"),
         plot.margin = unit(c(0, 2, 0, 10), "mm"),
         axis.title.x =  element_text(size=12, vjust = 1,  color = "black"))+
  facet_wrap_custom(.~facet_marker, scales="free",nrow = 1, 
                    scale_overrides = list(
                      scale_override(1, scale_x_continuous(breaks = c(2184029/1E6, 2622125/1E6, 3280545/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(2, scale_x_continuous(breaks = c(2876548/1E6, 3107587/1E6, 3717385/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(3, scale_x_continuous(breaks = c(15288009/1E6, 15883817/1E6, 16167354/1E6),labels = scales::number_format(accuracy = 0.1)))
                    ))


#### cow s5 ##

figS5_C <- cowplot::plot_grid(figS5C_ref ,figS5C_alt,
                             nrow = 2,ncol = 1,
                             axis = "lr",
                             rel_heights = c(4.3,1))




figS5 <- cowplot::plot_grid(figS5_A ,figS5_B,figS5_C,
                                    nrow = 3,ncol = 1,
                                    labels = c("A","B","C"), 
                                    label_size = 12, 
                                    label_fontfamily="Helvetica",
                                    axis = "lr",
                                    rel_heights = c(1,0.8,1.6))




ggsave(figS5, filename = paste( "figures/Fig_S5.png",sep=""), units = "mm",height = 180, width = 170)




##########################################
#           Figure S6                    #
#      divergent strans GWA              #
##########################################


## Figure S6 manhanton #####

data_SS6 <- read.csv("processed_data/SS6_File.csv", stringsAsFactors=FALSE)

figS6 <- manhaplot(data_SS6) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5.8)) 

ggsave(figS6, filename = paste( "figures/Fig_S6.png",sep=""), units = "mm",height = 60, width = 170)





##########################################
#           Figure S7                    #
#          norm.n GWA                   #
##########################################


## Figure S7 A manhanton #####

data_SS7A <- read.csv("processed_data/SS7A_File.csv", stringsAsFactors=FALSE)

figS7_A <- manhaplot(data_SS7A) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 


## Figure S7 B pxg #####


data_SS7B <- read.csv("processed_data/SS7B_File.csv", stringsAsFactors=FALSE)


figS7_B <- pxg_plot(data_SS7B) +
  theme(legend.position = "none") 


## Figure S7 C interval haplotype #####

data_SS7C <- read.csv("processed_data/SS7C_File.csv", stringsAsFactors=FALSE)

figS7C_ref <- iv_hap(subset(data_SS7C,ALL=="REF")) +
  ylab("REF") +
  theme( plot.margin = unit(c(0, 2, 0, 10), "mm"))+
  facet_wrap_custom(.~facet_marker, scales="free",nrow = 1, 
                    scale_overrides = list(
                      scale_override(1, scale_x_continuous(breaks = c(9529207/1E6, 10256873/1E6, 10691570/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(2, scale_x_continuous(breaks = c(3673426/1E6, 4831537/1E6, 6132617/1E6),labels = scales::number_format(accuracy = 0.1)))
                    ))

figS7C_alt <- iv_hap(subset(data_SS7C,ALL=="ALT")) +
  ylab("ALT")  +
  theme( strip.text = element_blank(),
         axis.text =  element_text(size=12,  color = "black"),
         plot.margin = unit(c(0, 2, 0, 10), "mm"),
         axis.title.x =  element_text(size=12, vjust = 1,  color = "black"))+
  facet_wrap_custom(.~facet_marker, scales="free",nrow = 1, 
                    scale_overrides = list(
                      scale_override(1, scale_x_continuous(breaks = c(9529207/1E6, 10256873/1E6, 10691570/1E6),labels = scales::number_format(accuracy = 0.1))),
                      scale_override(2, scale_x_continuous(breaks = c(3673426/1E6, 4831537/1E6, 6132617/1E6),labels = scales::number_format(accuracy = 0.1)))
                    ))


#### cow s7 ##

figS7_C <- cowplot::plot_grid(figS7C_ref ,figS7C_alt,
                              nrow = 2,ncol = 1,
                              axis = "lr",
                              rel_heights = c(4,1))




figS7 <- cowplot::plot_grid(figS7_A ,figS7_B,figS7_C,
                            nrow = 3,ncol = 1,
                            labels = c("A","B","C"), 
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr",
                            rel_heights = c(1,1,1.5))




ggsave(figS7, filename = paste( "figures/Fig_S7.png",sep=""), units = "mm",height = 180, width = 170)


##########################################
#           Figure S8                    #
#      norm.n linkage-mapping 1% H2O     #
##########################################


load("processed_data/SS8A_File.RData")

figS8_A <- lk_lod_plot(lkmap_SS8,cis_SS8)



data_SS8B <- read.csv("processed_data/SS8B_File.csv", stringsAsFactors=FALSE)


data_SS8B$peakmarker <- factor(data_SS8B$marker, levels = c(" Parental", "II:3646789"))


figS8_B <- lk_pxg_plot(data_SS8B) + 
  ggplot2::facet_wrap(~peakmarker, ncol = 2, scales = "free_x") 


figS8 <- cowplot::plot_grid(figS8_A,figS8_B,
                             labels = c("A","B"), 
                            nrow = 2,
                            ncol = 1,
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr")


ggsave(figS8, filename = paste( "figures/Fig_S8.png",sep=""), units = "mm",height = 120, width = 170)



##########################################
#           Figure S9                    #
#      norm.n linkage-mapping 0.5% DMSO  #
##########################################


load("processed_data/SS9A_File.RData")

figS9_A <- lk_lod_plot(lkmap_SS9,cis_SS9)



data_SS9B <- read.csv("processed_data/SS9B_File.csv", stringsAsFactors=FALSE)


data_SS9B$peakmarker <- factor(data_SS9B$marker, levels = c(" Parental", "II:8964146" , "IV:8044494" , "IV:16522256" ,"V:11883496"  ))



figS9_B <- lk_pxg_plot(data_SS9B) + 
  ggplot2::facet_wrap(~peakmarker, ncol = 5, scales = "free_x") + 
  ggplot2::labs(x = "", y = paste("norm.n","(0.5% DMSO)",sep = " "))



figS9 <- cowplot::plot_grid(figS9_A,figS9_B,
                            labels = c("A","B"), 
                            nrow = 2,
                            ncol = 1,
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr")


ggsave(figS9, filename = paste( "figures/Fig_S9.png",sep=""), units = "mm",height = 120, width = 170)




##########################################
#           Figure S10                    #
#      norm.n linkage-mapping 1% DMSO  #
##########################################


load("processed_data/SS10A_File.RData")

figS10_A <- lk_lod_plot(lkmap_SS10,cis_SS10)



data_SS10B <- read.csv("processed_data/SS10B_File.csv", stringsAsFactors=FALSE)


data_SS10B$peakmarker <- factor(data_SS10B$marker, levels = c(" Parental",  "IV:9031007" ,"V:11883496"  ))



figS10_B <- lk_pxg_plot(data_SS10B) + 
  ggplot2::facet_wrap(~peakmarker, ncol = 3, scales = "free_x") + 
  ggplot2::labs(x = "", y = paste("norm.n","(1% DMSO)",sep = " "))



figS10 <- cowplot::plot_grid(figS10_A,figS10_B,
                            labels = c("A","B"), 
                            nrow = 2,
                            ncol = 1,
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr")


ggsave(figS10, filename = paste( "figures/Fig_S10.png",sep=""), units = "mm",height = 120, width = 170)
