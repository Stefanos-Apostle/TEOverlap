

LDClump_window <- function(LD_df, win) {
  "
  Function to alter the window of the LD Clump to a specific range

  Inputs:
    LD_df <- GWAS clump data frame
    win <- BP range to extend position left and right

  Output:
    Edited LD_df with new clumpStart and clumpEnd coordinates based on window
  "
  LD_df$clumpStart <- LD_df$position - win
  LD_df$clumpEnd <- LD_df$position + win

  return(LD_df)
}



read_in_TEage <- function(file) {
  "
  Reads in .bed file of TE ages and returns as a data.frame with foramtted column names and strsplit chromosome names.
  "

  ERV_age_df <- read_delim(file, delim = "\t", col_names = F)

  if (ncol(ERV_age_df) == 6) {
    colnames(ERV_age_df) <- c("chromosome", "start", "end", "name", "age", "strand")
  }else {
    warning("File must contain only 6 colums in the following order; c('chromosome', 'start', 'end', 'name', 'age', 'strand')")
  }

  if (class(ERV_age_df$chromosome) == "character") {
    chr_num <- unlist(lapply(ERV_age_df$chromosome, FUN = function(x){strsplit(x, "chr")[[1]][2]}))

    if (length(unique(ERV_age_df$chromosome)) == length(unique(chr_num))) {
      ERV_age_df$chromosome <- chr_num
    }else{
      warning("Renaming chromosome number did not work properly. Please make sure your chromosomes column are in the following formal; 'chr#'")
    }
  }else{
    warning("Chromosome naming must be in the following format; ie 'chr#")
  }

  return(ERV_age_df)
}




TE_overlap <- function(TE, TE_age_df, LD_df) {
  "
  Finds TEs that overlap a specific SNP position.

  Inputs:
    TE <- Name of TE from variants DF colums 4; ie L1MdFanc_I/LINE/L1
    TE_age_df <- Data.frame of specific TE elements with their ranges and associated age (MYA)
    LD_df <- GWAS clumps data frame

  Output:
    Returns a data.frame showing overlap of SNPs in TEs.
  "

  TE_oi <- TE_age_df[which(TE_age_df$name == TE),]




  TE_LD_df <- data.frame()
  for (i in c(1:nrow(TE_oi))) {
    rdf <- TE_oi[i, ]
    chr_LD_df <- LD_df[which(LD_df$chromosome == rdf$chromosome & LD_df$position > rdf$start & LD_df$position < rdf$end), ]

    if (nrow(chr_LD_df) > 0) {
      chr_LD_df$TE <- rdf$name
      chr_LD_df$age <- rdf$age
      if (nrow(TE_LD_df) == 0) {
        TE_LD_df <- chr_LD_df
      }else{
        TE_LD_df <- rbind(TE_LD_df, chr_LD_df)
      }
    }
  }

  return(TE_LD_df)
}


runOverlap <- function(TE_age_df, LD_df) {
  "
  Runs TE_overlap for each TE from age_df to find if it overlaps with a single SNP position

  Inputs:
    TE_age_df <- data.frame of TE and assocaited ages from laoded in .bed file
    LD_df <- data.frame GWAS clumps and SNP positions

  Output:
    Returns data.frame of SNPs which are in TEs and the associated TE name and age

  "
  start.time <- Sys.time()

  TE_names <- unique(TE_age_df$name)

  TEi = lapply(TE_names, function(x){TE_overlap(x, TE_age_df, LD_df)})
  overlap_df <- do.call(rbind, TEi)

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste("Time Taken:", time.taken))

  return(overlap_df)
}





TE_overlap_clumps <- function(TE_age_df, LD_df) {
  "
  Finds TEs that overlap our LD clump windows.

  Inputs:
    TE_age_df <- Data.frame of specific TE elements with their ranges and associated age (MYA)
    LD_df <- GWAS clumps data frame

  Output:
    Returns a data.frame showing overlap of SNPs in TEs.
  "

  ## add columns for length of TE, length of clump, length of overlap, and % coverage

  #TE_oi <- TE_age_df[which(TE_age_df$name == TE),]

  print("Running TE - Clump Overlap")
  TE_LD_df <- data.frame()
  for (i in c(1:nrow(LD_df))) {
    rdf <- LD_df[i, ]
    chr_TE_df <- TE_age_df[which(TE_age_df$chromosome == rdf$chromosome), ]
    chr_TE_df <- chr_TE_df[-which(chr_TE_df$start > rdf$clumpEnd | chr_TE_df$end < rdf$clumpStart), ]

    if (nrow(chr_TE_df) > 0) {

      chr_TE_df <- cbind(chr_TE_df, rdf[,-which(colnames(rdf) %in% c("alignment", "chromosome"))])

      if (nrow(TE_LD_df) == 0) {
        TE_LD_df <- chr_TE_df
      }else{
        TE_LD_df <- rbind(TE_LD_df, chr_TE_df)
      }
    }

    log <- (i / nrow(LD_df)) *100
    if (log %in% seq(0, 100, by = 10)) {
      print(paste(log, "% Complete...", sep = ""))
    }
  }
  #print("Calculating % TE Coverage...")
  #TE_LD_df$perc_cover <- calc_overlap_percent(TE_LD_df)
  print("Done!")
  return(TE_LD_df)
}


calc_axis_set <- function(alt) {
  "
  Function to format GWAS_df correctly

  Input:
    alt = GWAS_df either as whole or subsetted by subfamily

  Output:
    Formatted data frame for coroboration with other functions
  "
  alt$chromosome[which(alt$chromosome == "X")] <- 20
  alt$chromosome[which(alt$chromosome == "Y")] <- 21
  alt$chromosome[which(alt$chromosome == "M")] <- 22
  alt$chromosome <- as.numeric(alt$chromosome)

  axis_set <- data.frame()
  for (i in unique(alt$chromosome)) {
    cdf = alt[which(alt$chromosome == i), ]
    cdf <- cdf[order(cdf$position, decreasing = F), ]
    if (nrow(axis_set) == 0) {
      axis_set <- cdf
    }else{
      axis_set <- rbind(axis_set, cdf)
    }
  }

  genes <- c()
  for (i in axis_set$nearest) {
    a <- strsplit(i, ",")[[1]]
    genes <- c(genes, a[length(a)])
  }
  axis_set$gene <- genes
  return(axis_set)
}



plot_manhattan <- function(GWAS_df, subfamily = "All", label_cutoff = 75, size = "Pvalue"){
  "
  Manhattan plot of the SNPs that contains overlapping TEs

  Inputs:
    GWAS_df = Output data frame of overlap functions
    subfamily = Option to subset data frame to plot only a specific subfamily of TEs; default is to plot All TEs
    label_cutoff = -log10(pValue) cutoff to plot labels on the manhattan plot
    size = Selecting the size variable for the plot; either by Pvalue or Age

  Output:
    Ggplot2 formatted manhattan plot
  "


  if (subfamily == "All") {
    alt <- GWAS_df
  }else{
    alt <- GWAS_df[which(GWAS_df$subfamily == subfamily), ]
  }

  axis_set <- calc_axis_set(alt)

  if (size == "Pvalue") {
    p <- ggplot(axis_set, aes(x = position, y = -log10(pValue),size = -log10(pValue)))
  }else if (size == "Age") {
    p <- ggplot(axis_set, aes(x = position, y = -log10(pValue),size = (1/age)))
  }


  p +
    facet_grid(. ~ chromosome, scales = "free_x", switch = "x")  +
    geom_jitter(aes( colour = age))+
    theme_classic() +
    geom_text_repel(size = 2, aes(label = ifelse((-log10(pValue) > label_cutoff), gene, "")), point.padding = 1) +
    labs(color = "Age (MYA)") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    #scale_size_continuous(range = c(1,2)) +
    ggtitle(paste("SNPs in ", subfamily, sep = "")) +
    theme(text = element_text(face = "bold", size = 15), plot.title = element_text(hjust = 0.5), line = element_line(size = 1)) +
    coord_cartesian(clip = 'off') +
    scale_color_viridis_c()
}



corr_plot <- function(GWAS_df, subfamily = "All", Age = "raw", TE_age_df = "not required for raw age") {
  "
  Function to plot a scatter plot of the SNPs by age and pValue as well as calculate correlation between the two variables.

  Inputs:
    GWAS_df = Output data frame of overlap functions
    subfamily = Option to subset data frame to plot only a specific subfamily of TEs; default is to plot All TEs
    Age = Option to plot raw age values or relative subfamily age by z score to compare cross-subfamily overlapping TEs
    TE_age_df <- Data.frame of specific TE elements with their ranges and associated age (MYA)

  Output:
    Scatter plot with correlation results
  "
  if (subfamily == "All") {
    alt <- GWAS_df
    if (Age == "relative") {
      alt <- calc_subfamily_agez(alt, TE_age_df)
    }else if (Age != "raw") {
      stop("Age must be one of the following; c('raw', 'relative')")
    }
  }else{
    alt <- GWAS_df[which(GWAS_df$subfamily == subfamily), ]
  }

  if (Age == "relative" && TE_age_df == "not required for raw age") {
    stop("If running relative age, TE_age_df must be supplied as input")
  }

  axis_set <- calc_axis_set(alt)

  rank_df <- data.frame(Row = order(-log10(axis_set$pValue), decreasing = T),
                        Age = if (Age == "raw") {
                          axis_set$age[order(-log10(axis_set$pValue), decreasing = T)]
                        }else if (Age == "relative") {
                          axis_set$subf_z[order(-log10(axis_set$pValue), decreasing = T)]
                        },
                        nlog10p = -log10(axis_set$pValue[order(-log10(axis_set$pValue), decreasing = T)]),
                        gene = axis_set$gene[order(-log10(axis_set$pValue), decreasing = T)],
                        Rank = c(1:nrow(axis_set)))

  lm <- lm(Age ~ Rank, rank_df)

  ggscatter(rank_df, x = "nlog10p", y = "Age",
            add = "reg.line", conf.int = F,
            cor.coef = T) +
    #annotate("text", x = 45, y = 115, label = paste("y = ", round(lm$coefficients[2], digits = 4), "x + ", round(lm$coefficients[1], digits = 4), sep = "")) +
    ggtitle(paste("SNPs of ", subfamily, sep = "")) +
    ylab(ifelse(Age == 'raw', 'Age (MYA)', 'Age (z-score by subfamily)')) +
    xlab("-log10pValue") +
    theme(text = element_text(face = "bold", size = 15), plot.title = element_text(hjust = 0.5), line = element_line(size = 1))

}





nTE_genes <- function(GWAS_df, labels_repelled = F) {
  "
  Function to plot the frequencies of overlapping TEs with their associated genes.

  Input:
    GWAS_df = Output data frame of overlap functions
    labels = Selecting whether or it labels should be above bars or not

  Output:
    Ggplot2 formatted bar plot
  "

  axis_set <- calc_axis_set(GWAS_df)

  dup <- c()
  age <- c()
  for (i in unique(axis_set$gene)) {
    if (length(which(axis_set$gene == i)) > 1){
      dup <- c(dup, i)
      age <- c(age,mean(axis_set$age[which(axis_set$gene == i)]))
    }
  }
  dup_df <- data.frame(gene = dup,
                       age = age)

  asg <-data.frame()
  for (i in unique(axis_set$gene)) {
    tdf <- data.frame(Gene = i,
                      Count = length(which(axis_set$gene == i)),
                      max_nlogp = max(-log10(axis_set$pValue[which(axis_set$gene == i)])))
    if (nrow(asg) == 0) {
      asg <- tdf
    }else {
      asg <- rbind(asg, tdf)
    }
  }

  p <- ggbarplot(asg, x = "Gene", y = "Count", color = "Gene")+
    theme(text = element_text(face = "bold", size = 15), plot.title = element_text(hjust = 0.5), line = element_line(size = 1))  +
    theme(legend.position = "none") +
    ggtitle("Count of Regions in TEs Closest Gene")

  if (labels_repelled == T) {
    p <- p +
      theme(axis.text.x = element_blank()) +
      geom_text_repel(label = ifelse(asg$Count > 1, asg$Gene, ""), color = ifelse(asg$max_nlogp > 100, "red", "black"))
  } else{
    p <- p +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  return(p)
}


perc_lowdecile <- function(GWAS_df, d = 1) {
  "
  Function to calcualte the % enrichment of TEs in the bottom d decile(s) of age.

  Input:
    GWAS_df = Output data frame of overlap functions
    d = Number of bottom deciles to calculate enrichment

  Output:
    Data frame of each subfamily describing enrichment in the bottom decile of age.
  "
  GWAS_df <- subfamily_from_TE(GWAS_df)

  xr <- data.frame()
  for (i in unique(GWAS_df$subfamily)) {
    subr <- GWAS_df[which(GWAS_df$subfamily == i), ]
    mi <- min(subr$age)
    ma <- max(subr$age)
    low_dec <- (((ma - mi) / 10)*d)  + mi
    x <- length(which(subr$age < low_dec)) / length(subr$age)
    #names(x) <- i
    tdf <- data.frame(Subfamily = i,
                      PercDecile = x*100,
                      NSubFam = length(subr$age)
    )
    if (nrow(xr) == 0) {
      xr <- tdf
    }else {
      xr <- rbind(xr, tdf)
    }
  }

  return(xr)
}



prop_overlap <- function (LD_df, GWAS_df) {
  "\n  Function to calculate and plot the proportion of LD Clumps with overlapping TEs\n\n  Inputs:\n    LD_df <- GWAS clumps data frame\n    GWAS_df <- Output data frame of overlap functions\n\n  Output:\n    GGplot2 formatted pie chart of percent LD clumps with/without overlapping TEs\n  "

  if ("varId" %in% colnames(LD_df)) {
    cn <- "varID"
  }else if ("refsnp_id" %in% colnames(LD_df)) {
    cn <- "refsnp_id"
  } else {
    stop("VarID or refsnp_id must be present in the colnames of LD_df")
  }

  ol1 <- LD_df[, which(colnames(LD_df) == cn)]
  ol2 <- GWAS_df[, which(colnames(GWAS_df) == cn)]

  snp_df <- data.frame(No_TE_Overlap = (length(ol1) -
                                          length(which(ol1 %in% ol2)))/length(ol1),
                       TE_Overlap = length(which(ol1 %in% ol2))/length(ol1))
  snp_df <- as.data.frame(t(snp_df))
  snp_df$column <- c("No TE Overlap", "TE Overlapping SNPs")
  snp_df$label <- paste(round(snp_df$V1 * 100, digits = 1),
                        "%", sep = "")
  p <- ggpie(snp_df, x = "V1", label = "label", fill = "column",
             color = "white", lab.pos = "in") + theme(legend.position = "right") +
    theme(text = element_text(face = "bold", size = 15),
          plot.title = element_text(hjust = 0.5), line = element_line(size = 1)) +
    ggtitle("Proportion of Regions \n with Overlapping TEs")
  return(p)
}


prop_subfamilies <- function(GWAS_df){
  "
  Function to calcualte and plot the proportion of each TE subfamily in the TEs that overlap LD clumps

  Input:
    GWAS_df <- Output data frame of overlap functions

  Output:
    GGPlot2 formatted pie chart showing percent of each subfamily represented in overlapping TEs
  "

  GWAS_df$subfamily <- subfamily_from_TE(GWAS_df)

  pie_df <- data.frame()
  for (t in unique(GWAS_df$subfamily)){
    tdf <- data.frame(Subfamily = t,
                      Percent = (length(which(GWAS_df$subfamily == t)) / length(GWAS_df$subfamily))*100)
    if (nrow(pie_df) == 0) {
      pie_df <- tdf
    }else{
      pie_df <- rbind(pie_df, tdf)
    }
  }
  pie_df$label = ifelse(pie_df$Percent > 2, round(pie_df$Percent, digits = 1), "")
  #pie_df$label = ifelse(pie_df$Percent > 2, paste(round(pie_df$Percent, digits = 1), pie_df$Subfamily, sep = " \n "), "")

  pie_plot <- ggpie(pie_df, "Percent", label = "label", fill = "Subfamily",
                    color = "white", lab.pos = "out", lab.font = 3, lab.adjust = 5) +
    ggtitle("Percent of TE Subfamilies in TE_BMI_SNPs") +
    theme(legend.position = "right") +
    theme(text = element_text(face = "bold", size = 15), plot.title = element_text(hjust = 0.5), line = element_line(size = 1))

  p1 <- ggpie(pie_df, "Percent", label = "label", fill = "Subfamily",
        color = "white", lab.pos = "out", lab.font = 1, lab.adjust = 5) +
    ggtitle("Percept Proportion of TE Subfamilies") +
    theme(legend.position = "none") +
    theme(text = element_text(face = "bold", size = 15), plot.title = element_text(hjust = 0.5), line = element_line(size = 1))


  leg <- get_legend(pie_plot)
  p2 <- as_ggplot(leg)

  p3 <- ggarrange(p1, p2, ncol = 2)
  return(p3)
}




calc_subfamily_agez <- function(GWAS_df, TE_age_df){
  "
  Function to calculate the relative age of each overlapping TE by subfamily to compare all TEs

  Input:
    GWAS_df = Output data frame of overlap functions
    TE_age_df <- Data.frame of specific TE elements with their ranges and associated age (MYA)

  Output:
    Updated GWAS_df with column subf_z for realtive ages (z score by subfamily)
  "

  GWAS_df <- subfamily_from_TE(GWAS_df)

  #calculating mean and sd for each subfamily
  sub_stats <- data.frame()
  for (s in unique(GWAS_df$subfamily)) {
    m <- mean(TE_age_df$age[grep(s, TE_age_df$name)])
    sd <- sd(TE_age_df$age[grep(s, TE_age_df$name)])
    tdf <- data.frame(subfamily = s,
                      mean = m,
                      sd = sd)
    if (nrow(sub_stats) == 0) {
      sub_stats <- tdf
    }else{
      sub_stats <- rbind(sub_stats, tdf)
    }
  }

  # calculating z score of subfamily z-score for each overlapping TE
  agez <- c()
  for (r in c(1:nrow(GWAS_df))) {
    sf <- GWAS_df$name[r]

    for (k in sub_stats$subfamily) {
      gs <- grep(k, sf)
      if (length(gs) == 1) {
        ss <- which(sub_stats$subfamily == k)
      }
    }
    #ss <- grep(sub_stats$subfamily, sf)
    ta <- (GWAS_df$age[r] - sub_stats$mean[ss]) / sub_stats$sd[ss]
    agez <- c(agez, ta)
  }

  GWAS_df$subf_z <- agez
  return(GWAS_df)
}


subfamily_from_TE <- function(GWAS_df) {
  "
  Function to generate column for TE subfamilies from TE names

  Input:
    GWAS_df = Output data frame of overlap functions

  Output:
    Updated GWAS_df with subfamilies in column name 'subfamily'
  "
  TE_family <- c()
  for (i in GWAS_df$name) {
    if("?" %in% strsplit(i, "")[[1]]) {
      i <- strsplit(i, "")[[1]]
      i <- i[-which(i == "?")]
      i <- paste(i, collapse = "")
    }
    f = strsplit(i, "/")[[1]]
    if (length(f) == 2) {
      TE_family <- c(TE_family, f[2])
    }else{
      TE_family <- c(TE_family, paste(f[2], f[3], sep = "/"))
    }
  }
  GWAS_df$subfamily <- TE_family
  return(GWAS_df)
}



distance_plot <- function(SNP, TE_age_df, window = 500, color_by_subfamily = T, plot_expected = F) {
  "
  Function to plot the density of TEs from the SNP

  Inputs:
    SNP = Single row data frame of a single SNP (one row of LD_df)
    TE_age_df <- Data.frame of specific TE elements with their ranges and associated age (MYA)
    window = +- base pair window from SNP position to plot density of TEs
    color_by_subfamily = Whether or not the plot should be colored by TE subfamilies

  Output:
    Ggplot formatted plot of density of TEs around the SNP position

  "
  TE_age_df$position <- ((TE_age_df$end - TE_age_df$start)/2 + as.numeric(TE_age_df$start))
  chr_df <- TE_age_df[which(TE_age_df$chromosome == SNP$chromosome), ]

  # plotting window to search for TEs within this distance from the SNP
  dist_df <- chr_df[which(chr_df$position < (SNP$position + window) & chr_df$position > (SNP$position - window)) , ]

  if (color_by_subfamily == T) {dist_df <- subfamily_from_TE(dist_df)}

  if (color_by_subfamily == F) {
    p <- ggplot(dist_df, aes(x = position, color = "black")) +
      scale_color_manual(name = "Density", values = c("black" = "black", "red" = "red"), labels = c("Observed", "Expected"))
  }else{
    p <- ggplot(dist_df, aes(x = position, color = subfamily)) +
      labs(color = "Subfamily") +
      ggtitle(SNP$refsnp_id)
  }

  if (plot_expected == T) {
    exp_df <- calc_expected_positions(SNP, TE_age_df, window)

    kst <- ks.test(dist_df$position, exp_df$exp_position)
    p <- p +
      geom_density(data = exp_df, aes(x = exp_position, color = "red")) +
      ggtitle(paste(SNP$refsnp_id, paste("K-S pval = ", round(kst$p.value, digits = 3)), sep = "\n"))
  }

  p + geom_density() +
    geom_rug(size = 1) +
    geom_vline(xintercept = SNP$position, color = "black", linetype = "dashed") +
    xlim(c(SNP$position-window, SNP$position+window)) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          legend.position = "right") +
    xlab(paste("Chromosome", dist_df$chromosome))

}



calc_expected_positions <- function(SNP, TE_age_df, window, perc_genome_coverage = .45) {
  "
  Function to generate the expected TE insertions around a SNP

  Inputs:
    SNP = Single row data frame of a single SNP (one row of LD_df)
    TE_age_df <- Data.frame of specific TE elements with their ranges and associated age (MYA)
    window = +- base pair window from SNP position to plot density of TEs
    perc_genome_coverage = Average coverage of TEs in the genome from literature (45%)

  Output:
    Data frame of expected TE positions around the SNP to use for density plot
  "
  avg_TE_length <- mean(TE_age_df$end - TE_age_df$start)
  num_expected_TE_insertions <- round(((window*2)*perc_genome_coverage)/ avg_TE_length)
  sep <- (window*2)/num_expected_TE_insertions
  exp_pos <- unlist(lapply(c(1:num_expected_TE_insertions), FUN = function(x) {sep/2 + sep*(x-1)}))
  exp_pos <- exp_pos + (SNP$position-window)
  ep_df <- data.frame(exp_position = exp_pos)
  return(ep_df)
}















