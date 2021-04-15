



##### function for custom facet #### 
#https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/

scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)


facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}





### linkage pxg plot ####
pxgplot_par_kt <- function (cross, map, parpheno, tit = "", ylab = "norm.n",
                            textsize = 12, titlesize = 12, pointsize = 0.5) {
  peaks <- map %>% 
    dplyr::group_by(iteration) %>% 
    dplyr::filter(!is.na(var_exp)) %>% 
    dplyr::do(head(., n = 1))
  
  if (nrow(peaks) == 0) {
    stop("No QTL identified")
  }
  
  uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
  # colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
  
  pheno <- cross$pheno %>% 
    dplyr::select(map$trait[1])
  
  geno <- data.frame(linkagemapping:::extract_genotype(cross)) %>% 
    dplyr::select(which(colnames(.) %in% uniquemarkers)) %>% 
    data.frame(., pheno)
  
  colnames(geno)[1:(ncol(geno) - 1)] <- sapply(colnames(geno)[1:(ncol(geno) - 1)], 
                                               function(marker) {paste(unlist(peaks[peaks$marker == gsub("\\.", "-", marker),
                                                                                    c("chr", "pos")]), collapse = ":")})
  colnames(geno)[ncol(geno)] <- "pheno"
  
  split <- tidyr::gather(geno, marker, genotype, -pheno) %>%
    tidyr::drop_na(genotype)
  
  split$genotype <- sapply(split$genotype, function(x) {
    if (x == -1) {
      "N2-RIAILs"
    }
    else {
      "CB-RIAILs"
    }
  })
  
  # add parent phenotype
  parpheno <- parpheno %>% 
    dplyr::mutate(marker = " Parental") %>%
    dplyr::mutate(genotype = strain) %>%
    dplyr::select(pheno = phenotype, marker, genotype)
  
  split <- split %>%
    dplyr::bind_rows(parpheno)
  
  split$genotype <- factor(split$genotype, levels = c("N2", "CB4856", "N2-RIAILs", "CB-RIAILs"))
  
  ggplot2::ggplot(split) + 
    ggplot2::geom_jitter(ggplot2::aes(x = genotype, y = pheno), alpha = 1, size = pointsize, width = 0.1) + 
    ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype, alpha = 0.5), outlier.shape = NA) + 
    ggplot2::scale_fill_manual(values = c(`N2-RIAILs` = "orange", `CB-RIAILs` = "blue", "N2" = "orange", "CB4856" = "blue")) + 
    ggplot2::facet_wrap(~marker, ncol = 5, scales = "free_x") + 
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = textsize, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = textsize, color = "black"), 
                   axis.title.x = ggplot2::element_text(size = titlesize, color = "black", vjust = -0.3), 
                   axis.title.y = ggplot2::element_text(size = titlesize,  color = "black"), 
                   strip.text.x = ggplot2::element_text(size = textsize,  color = "black"),
                   strip.text.y = ggplot2::element_text(size = textsize,  color = "black"),
                   plot.title = ggplot2::element_blank(), 
                   strip.background = element_blank(),
                   text=element_text(family="Helvetica"),
                   legend.position = "none", 
                   panel.grid = ggplot2::element_blank()) + 
    ggplot2::labs(x = "", y = ylab) + 
    scale_x_discrete(labels=c("CB-RIAILs" = "CB\nRIAILs", "N2-RIAILs" = "N2\nRIAILs"))
  
}



lk_pxg_plot <- function(lk_pxg_data){
  

  lk_pxg_data$genotype <- factor(lk_pxg_data$genotype, levels = c("QX1430", "CB4856", "QX-RIAILs", "CB-RIAILs"))
  
  lk_pxg_plt <- ggplot2::ggplot(lk_pxg_data) + 
    ggplot2::geom_jitter(ggplot2::aes(x = genotype, y = pheno), alpha = 1, size = 0.5, width = 0.1) + 
    ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype, alpha = 0.5), outlier.shape = NA) + 
    ggplot2::scale_fill_manual(values = c(`QX-RIAILs` = "orange", `CB-RIAILs` = "blue", "QX1430" = "orange", "CB4856" = "blue")) + 
    ggplot2::facet_wrap(~marker,  scales = "free_x") + 
    ggplot2::theme_bw() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12, color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 12, color = "black"), 
                   axis.title.x = ggplot2::element_text(size = 12, color = "black", vjust = -0.3), 
                   axis.title.y = ggplot2::element_text(size = 12,  color = "black"), 
                   strip.text.x = ggplot2::element_text(size = 12,  color = "black"),
                   strip.text.y = ggplot2::element_text(size = 12,  color = "black"),
                   plot.title = ggplot2::element_blank(), 
                   strip.background = element_blank(),
                   text=element_text(family="Helvetica"),
                   legend.position = "none", 
                   panel.grid = ggplot2::element_blank()) + 
    ggplot2::labs(x = "", y = "norm.n") + 
    scale_x_discrete(labels=c("CB-RIAILs" = "CB\nAllele", "QX-RIAILs" = "QX\nAllele"))
  
  return(lk_pxg_plt)
  
}



### linkage LOD plot ####
maxlodplot_kt <- function (map, textsize = 12, titlesize = 12, linesize = 0.5) {
  
  map<-mappingdf
  map1 <- map %>% dplyr::group_by(marker) %>% dplyr::filter(lod ==max(lod))
  cis <- map %>% dplyr::group_by(marker) %>% dplyr::mutate(maxlod = max(lod)) %>%
    dplyr::group_by(iteration) %>% dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n = 1))
  if (nrow(cis) == 0) {
    plot <- ggplot2::ggplot(map1, ggplot2::aes(x = genotype,y = pheno)) + ggplot2::geom_blank()
    return(plot)
  }
  map1 <- linkagemapping:::cidefiner(cis, map1)
  plot <- ggplot2::ggplot(map1) + 
    ggplot2::aes(x = pos/1e+06,y = lod)
  if (nrow(cis) != 0) {
    plot <- plot + ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e+06,ymin = 0, ymax = ci_lod), fill = "cyan",alpha=.3) +
      ggplot2::geom_point(data = cis, ggplot2::aes(x = pos/1e+06,y = (1.05 * maxlod)), fill = "red", shape = 25,
                          size = 3.2, show.legend = FALSE) + 
      ggplot2::geom_text(data = cis, ggplot2::aes(x = pos/1e+06, y = (1.25 * maxlod), label = paste0(100 *round(var_exp, digits = 3), "%")), 
                         colour = "black",size = 8*5 /14, hjust = "inward")
  }
  plot <- plot + ggplot2::geom_line(size = linesize, alpha = 0.85) +
    ggplot2::facet_grid(. ~ chr, scales = "free", space = "free") +
    ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD") +
    ggplot2::scale_colour_discrete(name = "Mapping\nIteration") +
    ggplot2::ggtitle(map1$trait[1]) + 
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = textsize, color = "black"),
      axis.text.y = ggplot2::element_text(size = textsize, color = "black"),
      axis.title.x = ggplot2::element_text(size = titlesize, color="black"),
      axis.title.y = ggplot2::element_text(size = titlesize, color="black"),
      strip.text.x = ggplot2::element_text(size = textsize,  color = "black"),
      strip.text.y = ggplot2::element_text(size = textsize, color = "black"),
      plot.title = ggplot2::element_blank(),
      strip.background = element_blank(),
      panel.grid = ggplot2::element_blank(),
      text=element_text(family="Helvetica")) + 
    scale_y_continuous(expand = expand_scale(mult = c(0, .05))) 
  return(plot)
}


lk_lod_plot <- function(lkmap,cis){
  
  lk_lod_plt <- ggplot2::ggplot(lkmap) + 
    ggplot2::aes(x = pos/1e+06,y = lod)+ 
    ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e+06,ymin = 0, ymax = ci_lod), fill = "cyan",alpha=.3) +
    ggplot2::geom_point(data = cis, ggplot2::aes(x = pos/1e+06,y = (1.05 * maxlod)), fill = "red", shape = 25,
                        size = 3.2, show.legend = FALSE) + 
    ggplot2::geom_text(data = cis, ggplot2::aes(x = pos/1e+06, y = (1.25 * maxlod), label = paste0(100 *round(var_exp, digits = 3), "%")), 
                       colour = "black",size = 8*5 /14, hjust = "inward") + 
    ggplot2::geom_line(size = 0.5, alpha = 0.85) +
    ggplot2::facet_grid(. ~ chr, scales = "free", space = "free") +
    ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD") +
    ggplot2::scale_colour_discrete(name = "Mapping\nIteration") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 12, color = "black"),
      axis.text.y = ggplot2::element_text(size = 12, color = "black"),
      axis.title.x = ggplot2::element_text(size = 12, color="black"),
      axis.title.y = ggplot2::element_text(size = 12, color="black"),
      strip.text.x = ggplot2::element_text(size = 12,  color = "black"),
      strip.text.y = ggplot2::element_text(size = 12, color = "black"),
      plot.title = ggplot2::element_blank(),
      strip.background = element_blank(),
      panel.grid = ggplot2::element_blank(),
      text=element_text(family="Helvetica")) + 
    scale_y_continuous(expand = expand_scale(mult = c(0, .05))) 
  
  return(lk_lod_plt)
  
}


### manhtton ####
manhaplot <- function (data_manh){  
  
#  data_manh <- data_SS11A
  
  plt_manh <-  ggplot2::ggplot(data_manh) +
  ggplot2::aes(x = POS/1e6, y = log10p) +
  ggplot2::scale_color_manual(values = c("0" = "black", 
                                         "1" = "red",
                                         "2" = "hotpink3")) +
 # ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6,    # this is the plot boundary for LD and gene plots
 #                                 xmax = endPOS/1e6,    # this is the plot boundary for LD and gene plots
 #                                 ymin = 0, 
  #                                ymax = Inf, 
  #                                fill = "plum4"), 
 #                    color = "plum4",fill = "cyan",linetype = 2, 
   #                  alpha=.3) +
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
                y = expression(-log[10](italic(p))))
  
  return(plt_manh)
  
}


### GWA pxg ####
pxg_plot <- function (data_pxg){  
  
  data_pxg$Genotypes <- factor(data_pxg$genotype,levels = c("swept","divergent"))
  
  
  plt_pxg <-  ggplot(data_pxg, aes(x=factor(as.character(allele),labels = c("REF","ALT")),y = value)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(shape=21,
                position=position_jitter(0.2),
                size=1.5,
                aes(fill=Genotypes),
                alpha=0.8) +                      
    scale_fill_manual(values=c("swept"="gold2","divergent"="plum4")) + 
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
    ylab("Lifetime fecundity")  +
    xlab("Genotype") +
    labs(fill="Swept ratio") 
  
  
  return(plt_pxg)
  
}


### interval haplotype ####

iv_hap <- function(data_inhap){
  

  data_inhap$Genotypes <- factor(data_inhap$genotype,levels = c("swept","divergent"))
  
  
  plt_iv_hap <- ggplot(data_inhap,
         aes(xmin = start_s/1E6, xmax = stop_s/1E6,
             ymin = -(plotpoint_hi + 0.5), ymax = -(plotpoint_hi - 0.5))) +
    geom_rect(show.legend = FALSE,
              aes(fill = swept_haplotype)) +
    scale_fill_manual(values = c("FALSE"="gray","TRUE"="red")) +
    scale_y_continuous(breaks = unique(data_inhap$plotpoint_hi),
                       expand = c(0, 0)) +
    geom_point(aes(x=(peakPOS/1E6),y=-plotpoint_hi,color=Genotypes),size=0.5) +
    scale_color_manual(values=c("swept"="gold2","divergent"="plum4")) + 
    xlab("Genomic position (Mb)") +
    theme_bw() +
    facet_wrap_custom(.~facet_marker, scales="free",nrow = 1, 
                      scale_overrides = list(
                        scale_override(1, scale_x_continuous(breaks = c(2183288/1E6, 2700622/1E6, 3050428/1E6),labels = scales::number_format(accuracy = 0.1))),
                        scale_override(2, scale_x_continuous(breaks = c(12433268/1E6, 14492150/1E6, 15150384/1E6),labels = scales::number_format(accuracy = 0.1)))
                      )) +
    theme(axis.text =  element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x =  element_blank(),
          axis.title.y =  element_text(size=12, vjust = 1,  color = "black"),
          strip.text = element_text(size=12, vjust = 1,  color = "black"),
          strip.background = element_blank(), 
          legend.position = "none",
          panel.grid = ggplot2::element_blank(),
          text=element_text(family="Helvetica"), #,
          plot.margin = unit(c(2, 2, 0, 1), "mm"))#,
        #  panel.spacing = unit(0.3, "lines")) 
  
  return(plt_iv_hap)
  
}

