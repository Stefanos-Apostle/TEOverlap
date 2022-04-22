

## Utils for distance from SNP -> nearest TE 


find_nearest_TE <- function(SNP, TE_age_df, direction = "upstream") {
  "
  Find distance to the nearest TE from SNP position
  
  Inputs:
    SNP = Single row data frame of a single SNP (one row of LD_df)
    TE_age_df <- Data.frame of specific TE elements with their ranges and associated age (MYA)
    direction = whether to find nearest upstream or downstream TE
    
  Output:
    Distance (bp) to nearest TE
  "
  
  chr_df <- TE_age_df[which(TE_age_df$chromosome == SNP$chromosome), ]
  cc_df <- chr_df
  cc_df$start[which(cc_df$strand == "-")] <- chr_df$end[which(chr_df$strand == "-")]
  cc_df$end[which(cc_df$strand == "-")] <- chr_df$start[which(chr_df$strand == "-")]
  
  if (direction == "upstream") {
    up_d <- SNP$position - cc_df$end
    up_d <- up_d[-which(up_d < 0)]
    md <- min(up_d)
    
    #n_up <- cc_df[which(abs(SNP$position - cc_df$end) == up_md), ]
  }else if (direction == "downstream") {
    down_d <- SNP$position - cc_df$start
    down_d <- down_d[-which(down_d > 0)]
    md <- min(abs(down_d))
    
    #n_down <- cc_df[which(abs(SNP$position - cc_df$start) == down_md), ]
  } else{
    stop("Direction must be one of the following; c('upstream', 'downstream')")
  }
  
  return(md)
}




plot_nearest_density <- function(LD_df, TE_age_df, directional = F, label = "") {
  "
  Function to plot the nearest upstream and downstream TE from SNPs by distance for a set of SNPs.
  
  Input:
    LD_df = GWAS clump data frame
    TE_age_df = Data.frame of specific TE elements with their ranges and associated age (MYA)
    directional = Boolean whether to separate upstream and downstream nearest TEs
    label = any character array to be included in the title
    
  Output:
    GGPlot of up and down-stream density plots to nearest TE from SNPs
  "
  
  up_dist <- c()
  down_dist <- c()
  for (i in c(1:nrow(LD_df))) {
    up_d <- find_nearest_TE(LD_df[i, ],  TE_age_df, direction = "upstream")
    up_dist <- c(up_dist, up_d)
    
    down_d <- find_nearest_TE(LD_df[i, ],  TE_age_df, direction = "downstream")
    down_dist <- c(down_dist, down_d)
  }
  
  up_df <- data.frame(dist = up_dist, Direction = "Upstream")
  down_df <- data.frame(dist = down_dist, Direction = "Downstream")
  
  comb_df <- rbind(up_df, down_df)
  comb_df$Direction <- factor(comb_df$Direction, levels = c("Upstream", "Downstream"))
  
  fig_theme <- theme_classic() +
    theme(text = element_text(face = "bold", size = 15), plot.title = element_text(hjust = 0.5), line = element_line(size = 1))
  
  
  if (directional){
    kst <- ks.test(comb_df$dist[which(comb_df$Direction == "Upstream")], 
                   comb_df$dist[which(comb_df$Direction == "Downstream")]
    )
    
    p <- ggplot(data = comb_df, aes(x = dist, color = Direction)) +
      scale_color_brewer(palette = "Dark2")+
      ggtitle(paste("Distance from SNP to Nearest TE",label, paste("K-S pval =", round(kst$p.value, digits = 3)), sep = "\n"))
  }else {
    p <- ggplot(data = comb_df, aes(x = dist))+
      ggtitle(paste("Distance from SNP to Nearest TE", label, sep = "\n"))
  }
  
  p <- p +
    geom_density(size = 1) +
    geom_rug() +
    xlab("Distance (bp)") +
    fig_theme
  
  return(p)
}



LD_nearest_dists <- function(LD_df, TE_age_df, directional = F) {
  "
  Plot the nearest distance to TE for each LD clump
  
  Inputs:
    LD_df = GWAS clump data frame
    TE_age_df = Data.frame of specific TE elements with their ranges and associated age (MYA)
    directional = Boolean whether to separate upstream and downstream nearest TEs
    
  Outputs:
    GGplot for each LD clump SNP distances to nearest TE
  "
  plot_list <- list()
  for (r in c(1:nrow(LD_df))) {
    SNP <- LD_df[r, ]
    db <- get_variants(genomic_range = list(chromosome = SNP$chromosome, start = SNP$clumpStart, end = SNP$clumpEnd))
    SNP_df <- db@variants
    
    colnames(SNP_df)[grep("position", colnames(SNP_df))] <- "position"
    colnames(SNP_df)[grep("chromosome", colnames(SNP_df))] <- "chromosome"
    
    ## subsetting TE_age_df to only allow TEs in the LD clump
    LD_TEs_range <- TE_age_df[which(TE_age_df$chromosome == SNP$chromosome), ]
    LD_TEs_range <- LD_TEs_range[which(LD_TEs_range$start < SNP$clumpEnd & LD_TEs_range$end > SNP$clumpStart), ]
    
    p <- plot_nearest_density(LD_df = SNP_df, TE_age_df = LD_TEs_range, directional = directional, label = SNP$refsnp_id)
    plot_list[[r]] <- p 
  }
  return (plot_list)
}



















































