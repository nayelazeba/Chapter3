library(DESeq2)
library(phyloseq)
library(tidyverse)
library(vegan)
library(readxl)
library(RColorBrewer)
library("ggpubr")
library(scales)

####### PRIMARY FUNCTIONS #######

# Reads in PS file, creates phyloseq object, appends data from excel file 
# inputs : 
#   rds_filename - PS file 
#   excel_filename - excel file containing soil data (can be null)
# returns : 
#   ps - with or without soil data depending on if excel is supplied.
get_ps <- function(rds_filename, excel_filename = NULL) {
  ps = read_in_rds(rds_filename)
  if (!is.null(excel_filename)) {
    soil_data = read_in_xl(excel_filename)
    samdat = sample_data(soil_data)
    sample_data(ps)=samdat  
  }
  return(ps)
}

# Generates ordinate plot based on supplied PS object
# inputs : 
#   ps - PS object (any normalization is OK, pruned or not is OK)
#   sample_type - Just for plot title, bacterial / fungal
#   shape_v, color_v, size_v - variables for ggplot (size is optional)
# returns : 
#   plot - plot object, run in order to gen plot. 

####### NMDS ################
get_nmds_vegan <- function(ps) {
  s = data.frame(sample_data(ps))
  psotu2veg <- function(ps) {
    OTU <- otu_table(ps)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  ps.otu.veg = psotu2veg(ps)
  dist <- vegdist(ps.otu.veg, method  = "bray")
  ps_nmds = metaMDS(dist, distance = "bray", k=3, trymax = 1000)
  return(ps_nmds)
}

plot_nmds <- function(ps_nmds, ps, temp, sample_type, color_v, shape_v, size_v=NULL, hjust=NULL, vjust=NULL) {
  color_title = title_case(color_v)
  shape_title = title_case(shape_v)
  points_nmds = data.frame(ps_nmds$points)
  plotNMDS = merge(points_nmds, data.frame(sample_data(ps)), by="row.names")
  plot_title = paste("NMDS plot of sampled soil", sample_type, "communities")
  color_v = sym(color_v)
  shape_v = sym(shape_v)
  stress_lab  = paste("stress =", round(ps_nmds$stress,2))
  annotations <- data.frame(
    xpos = Inf,
    ypos = Inf,
    annotateText = stress_lab,
    hjustvar = 1,
    vjustvar = 1)
  p_nmds = ggplot(plotNMDS, 
                  aes(color=as.factor(!!color_v), 
                      shape=as.factor(!!shape_v),
                      x= MDS1,
                      y= MDS2)) + 
    #ggtitle(plot_title) +
    theme_bw(base_size = 24) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"), 
          axis.ticks = element_blank(),
          axis.text = element_text(size = 18, colour = "black"),
          axis.title = element_text(size = 20, face = "bold"),
          plot.title = element_blank(),
          legend.title=element_text(size=20, face = "bold"), 
          legend.text=element_text(size=18)) +
    #annotate("text", x = -Inf, y = Inf, hjust=hjust, vjust=vjust, label = stress_lab , color="black", size=7)+
    #geom_text(data = data.frame(), aes(x=Inf, y=Inf, hjust="right",
    #                                   vjust="top", label=stress_lab,
    #                                   colour= "black", size = 10)) +
    labs(shape=shape_title, colour=color_title)
  
  if(temp == 350) {
    p_nmds = p_nmds + scale_color_manual(name = "Treatment", values = c("darkorange3", "black"), labels = c("soil", "soil + 350 PyOM"))+
      scale_shape_manual(name = "Time", values = c(0,15, 16, 17, 18), labels = c("Day 0", "Day 2", "Day 7", "Day 18", "Day 26")) 
  } else {
    p_nmds = p_nmds +  scale_color_manual(name = "Treatment", values = c("darkorange3", "black"), labels = c("soil", "soil + 550 PyOM")) +
      scale_shape_manual(name = "Time", values = c(0,15, 16, 17), labels = c("Day 0", "Day 4", "Day 10", "Day 22"))
  }

  if(is.null(size_v)) {
    p_nmds = p_nmds + geom_point(size=3, stroke = 1)
  } else {
    size_title = title_case(size_v)
    size_v = sym(size_v)
    p_nmds = p_nmds + geom_point(aes(size=as.factor(!!size_v))) + 
      labs(size=size_title)
  }
  return(p_nmds)
}


####### HELPER FUNCTIONS #######
read_in_rds <- function(filename) { 
  ps = readRDS(filename)
  return(ps)
}

read_in_xl <- function(filename) {
  soil_data = data.frame(read_excel(filename))
  rownames(soil_data) <- soil_data$Sample_ID
  soil_data$Sample_ID <- NULL
  return(soil_data)
}

prune_ps <- function(ps, ps_type) {
  if(ps_type == "bacteria") {
    ps <- prune_samples(sample_data(ps)$Sample_type != "pos_control" & sample_data(ps)$Sample_type != "blank", ps)
  } else if (ps_type == "fungi") {
    ps <- prune_samples(sample_data(ps)$Sample_type != "pos_control" & sample_data(ps)$Sample_type != "blank", ps)
  }
  return(ps)
}

normalize_ps <- function(ps, norm_type) {
  if(norm_type == "relative") {
    ps = transform_sample_counts(ps, function(x) x / sum(x))
  } else if (norm_type == "log") {
    ps = transform_sample_counts(ps, function(x) log(1 + x))
  } else if (norm_type == "sqrt") {
    ps = transform_sample_counts(ps, function(x) (x)^0.5)
  } 
  return(ps)
}

title_case <- function(y) {
  y <- sub("_", " ", y)
  return(y)
  #c <- strsplit(y, " ")[[1]]
  #paste(toupper(substring(c, 1,1)), substring(c, 2),
  #      sep="", collapse=" ")
}
