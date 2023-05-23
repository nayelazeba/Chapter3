library(data.table)
library(tidyverse)
library(phyloseq)
library(vegan)
#########################################################
############ 2020721_NZ complete LCMS run ###############
########### load data + inital processing  ##############
#########################################################
# Peak Height data exported from MS-DIAL as a .txt
# .txt file edited in Excel: 
# --removed first three rows of administrative junk from the MS
# --changed the name of column1 to "AlignmentID"
# --Most Importantly: changed all peak height values to scientific notation with 6 decimal places
# re-saved in Excel as a .csv file, and loaded here:

# msdat <- fread("./20220721_NZ_PeakHeights_FinalAnalysis.csv") #8061 features
# dim(msdat)
# #check the class of each column:
# msdat[ ,lapply(.SD, class)] #all peak height data columns should be numeric!

## Basic anatomy of the MS-DIAL output:
## rows = features or ions
## columns 1-32 = info about each feature
## columns 33-n = peak height values for each each feature in each sample
## columns at the end of the table with generalized sample names are averages followed by SD's
# peak height values are analogous to relative abundance in DNA/RNA sequencing data
# a peak height value >100000 is generally considered "real" 
# LC/MS data is inherently SUPER NOISEY because the MS is incredibly sensitive.
#
##
#
## INITITAL DATA PROCESSING

#create a reference table that lists column names & numbers
#msdat.colnames <- data.table(colnames(msdat)) 

#keep only the features that have a Peak Height value >100000 in any sample (exluding PooledQCs)
#msdat.sub1 <- msdat[apply(msdat[, c(33:227)], MARGIN = 1, function(x) any(x > 100000)), ] #6478 features

#keep only the features that have a Peak Height value <100000 in all Blank samples ("MeOH+Nonact")
#msdat.sub2 <- msdat.sub1[apply(msdat.sub1[, c(213:227)], MARGIN = 1, function(x) all(x < 100000)), ] #3546 features

# remove all features that with a Peak Height value >100000 in the Extract Controls 
#msdat.sub3 <- msdat.sub2[apply(msdat.sub2[, c(204:212)], MARGIN = 1, function(x) all(x < 100000)), ] #2624 features

#Export:
#write.csv(msdat.sub3, "/Users/monikafischer/Desktop/NayelaLCMS/20220721_NZ_FINAL.msdat.sub3.csv")

##################################################################
######### START HERE ########################################
## Read in the file
msdat.sub3 <- read.csv("20220721_NZ_FINAL.msdat.sub3.csv", row.names = 1, header = TRUE)
## Selecting only the columns with sample abundances :
msdat.PeakHeight <- msdat.sub3[,c(1,33:203)]
names(msdat.PeakHeight) <- gsub(".", "-", names(msdat.PeakHeight), fixed=TRUE)


## Run this code to extract features observed in soil-PyOM samples 
## Checking for duplicated columns
msdat.PeakHeight.nodupes <- data.frame(msdat.PeakHeight)[!duplicated(colnames(msdat.PeakHeight))]
## Removing features/rows with abundance >0 in soil samples
## Keeping features/rows that are present in soil-PyOM samples only
msdat.PeakHeight.res <- msdat.PeakHeight.nodupes %>%
  mutate(row_sums = rowSums(select(., matches("DOE.*S[0-9]")))) %>%
  filter(row_sums == 0) %>%
  dplyr::select(matches("DOE.*SB[0-9]|AlignmentId|row_sums"))

## Counting each feature/row present in soil-PyOM samples
msdat.PeakHeight.res <- data.frame(msdat.PeakHeight.res) %>%
  select(-"row_sums") %>%
  mutate(row_sums  = rowSums(select(., matches("DOE.*SB[0-9]"))> 0)) %>%
  arrange(-row_sums)
  
names(msdat.PeakHeight.res) <- gsub(".", "-", names(msdat.PeakHeight.res), fixed=TRUE)

## Merging the soil-PyOM data with the features data table 
msdat.features = msdat.sub3[, 1:32]

msdat.PeakHeight.merged <- msdat.PeakHeight.res %>%
  left_join(msdat.features, by = "AlignmentID")
#Export
#write.csv(msdat.PeakHeight.merged, "msdat.pyom.peaks.csv")

## For plotting ordinations and peaks
## Create an "OTUtable" of sample abundances:
msdat.PeakHeight.rn <- as.matrix(column_to_rownames(msdat.PeakHeight, var="AlignmentID"))
## Use msdat.PeakHeight df here if not correcting for soil only features.
## Use msdat.PeakHeight.res otherwise

#check for NAs, and replace any NAs with a zero:
sum(rowSums(is.na(msdat.PeakHeight.rn)) > 0) #no NAs! good!
#msdat.PeakHeight.rn[is.na(msdat.PeakHeight.rn)] <- 0

# Load metadata table:
msdat.SampleMetadata <- fread("./20220721_NZpilot_Metadata.csv")
#Comparing rownames and colnames before merging
msdat.SampleMetadata <- column_to_rownames(msdat.SampleMetadata, var="SampleID")
rownames(msdat.SampleMetadata)
colnames(msdat.PeakHeight.rn)
rownames(msdat.SampleMetadata) == colnames(msdat.PeakHeight.rn)

#combine tables into a phyloseq object:
msdat.ps <- phyloseq(otu_table(msdat.PeakHeight.rn, taxa_are_rows=TRUE), #matrix (rownames = AlignmentID, colnames = SampleNames, values = abundance)
                     sample_data(msdat.SampleMetadata)) #data.frame, (rownames = SampleNames, colnames & values = additional info for each SampleName)
# msdat.ps.350.S <- subset_samples(msdat.ps, PyOM =="350" & Treatment == "S")
# msdat.ps.350.SB <- subset_samples(msdat.ps, PyOM =="350" & Treatment == "SB350")
# msdat.ps.550.S <- subset_samples(msdat.ps, PyOM =="550" & Treatment == "S")
# msdat.ps.550.SB <- subset_samples(msdat.ps, PyOM =="550" & Treatment == "SB550")

#subset for samples in the big-kahuna/complete/full run:
msdat.ps.sub <- prune_samples(sample_data(msdat.ps)$Pilot == "full", msdat.ps)
### use this to run only timepoints t2, t3 and t4 for 350 samples
#subset for only the pilot samples that were part of the big-kahuna/complete/full run:
#msdat.ps.sub <- prune_samples(sample_data(msdat.ps)$Pilot == "full_PilotSamp", msdat.ps)
#subset for the pilot run:
#msdat.ps.sub <- prune_samples(sample_data(msdat.ps)$Pilot == "pilot", msdat.ps)
# subset for T1&T5 samples from the Pilot Run + T2,T3,T4 & 550 samples from the Complete Run:
#msdat.ps.sub <- prune_samples(sample_data(msdat.ps)$PilotSwap == "yes", msdat.ps)
# subset for full run (350:T2,T3,T4; 550:all and pilot rerun 350T1 and 350T5)
#msdat.ps.sub <- prune_samples(sample_data(msdat.ps.350)$Pilot != "pilot", msdat.ps.350)

#Remove T2 from 550 PyOM samples to stay consistent with no. of time points in the plots
msdat.ps.sub <- subset_samples(msdat.ps.sub, !(PyOM == "550" & Timepoint == "T2"))

#subset samples for temperature and treatment
msdat.ps.350 <- subset_samples(msdat.ps.sub, PyOM =="350")
msdat.ps.550 <- subset_samples(msdat.ps.sub, PyOM =="550")

## Plotting
#NMDS plots
# Defining functions to get ordinations using vegan
get_nmds <- function(ps) {
  s = data.frame(sample_data(ps))
  psotu2veg <- function(ps) {
    OTU <- otu_table(ps)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  ps.otu.veg = psotu2veg(ps)
  ps_nmds = metaMDS(ps.otu.veg, distance = "bray", k=2, trymax = 1000)
  return(ps_nmds)
}

plot_nmds <- function(ps_nmds, ps, temp, color_v, shape_v, size_v=NULL, hjust=NULL, vjust=NULL) {
  color_title = color_v
  shape_title = shape_v
  points_nmds = data.frame(ps_nmds$points)
  plotNMDS = merge(points_nmds, data.frame(sample_data(ps)), by="row.names")
  #plot_title = paste("NMDS plot of sampled soil", sample_type, "communities")
  color_v = sym(color_v)
  shape_v = sym(shape_v)
  stress_lab  = paste("stress =", round(ps_nmds$stress,2))
  dim_lab = paste("k =", ps_nmds$ndim)
  p_nmds = ggplot(plotNMDS, 
                  aes(color=as.factor(!!color_v), 
                      shape=as.factor(!!shape_v),
                      x= MDS1,
                      y= MDS2)) + 
    #ggtitle(plot_title) +
    theme_bw(base_size = 18) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"), 
          axis.ticks = element_blank(),
          axis.text = element_text(colour = "black"),
          axis.title = element_text(face = "bold"),
          plot.title = element_blank(),
          legend.title=element_text(face = "bold"), 
          legend.text=element_text(),
          legend.position = "right") +
    annotate("text", x = -Inf, y = Inf, hjust=hjust, vjust=vjust, label = c(stress_lab, dim_lab)  , color="black", size=7)+
    labs(shape=shape_title, colour=color_title)
  
  if(temp == 350) {
    p_nmds = p_nmds + scale_color_manual(name = "Treatment", values = c("steelblue", "gray54"), labels = c("soil", "soil + 350 PyOM"))+
      scale_shape_manual(name = "Timepoint", values = c(15, 16, 17,18), labels = c("Day 2", "Day 7", "Day 18")) 
    ###Day 0 and day 26 not included
  } else {
    p_nmds = p_nmds +  scale_color_manual(name = "Treatment", values = c("steelblue", "gray54"), labels = c("soil", "soil + 550 PyOM")) +
      scale_shape_manual(name = "Timepoint", values = c(15, 16, 17), labels = c("Day 4", "Day 10", "Day 22"))
    ####Day 0 (shape 0) not included
  }
  
  if(is.null(size_v)) {
    p_nmds = p_nmds + geom_point(size=4, stroke = 1)
  } else {
    size_title = title_case(size_v)
    size_v = sym(size_v)
    p_nmds = p_nmds + geom_point(aes(size=as.factor(!!size_v))) + 
      labs(size=size_title)
  }
  return(p_nmds)
}

# Runs adonis test on supplied ps with supplied formula
# inputs : 
#   ps - Must be a PS where the sample_data contains the columns you're referencing
#        in the formula. 
#   formula - Cannot be a string, must be a formula object (use as.formula(string))
# returns : 
#   None (prints results of adonis test to terminal)
run_adonis2 <- function(ps, formula) {
  ps_bray <- vegdist(t(as.data.frame(as.matrix(otu_table(ps)))))
  adonis2(formula, data = data.frame(sample_data(ps)))
}

#Plotting ordination using the function (Fig. 4A and 4B)
ps_nmds = get_nmds(msdat.ps.550)
plot = plot_nmds(ps_nmds, msdat.ps.550, 550, "Treatment", "Timepoint", hjust=-0.8, vjust=1.4)
#plot = plot + xlim(-1.0, 2.0)
plot
#ggsave("./nmds_lcms_550B_leg.png", dpi=300, device="png", width=6, height=6)

#PERMANOVA
ps_bray <- vegdist(t(as.data.frame(as.matrix(otu_table(msdat.ps.550)))))
formula = as.formula("ps_bray ~ Timepoint+Treatment")
run_adonis2(msdat.ps.550, formula)

#Plotting PyOM - specific peaks
features <- c("2384", "2634", "2758", "2126") ### these are features that are present in almost all PyOM-soil samples (both 350 and 550)
msdat.ps.melt <- psmelt(msdat.ps.sub)

#Creating confidence interval functions 
lower_se <- function(mean, se){
  res <- mean - se
  return(res)
}
upper_se <- function(mean, se){
  res <- mean + se
  return(res)
}

msdat.ps.peaks <- msdat.ps.melt %>%
  filter(OTU %in% features & Treatment != "S") ###soil samples have 0 abundance values, remove them

#Changing time point labels from "T1", "T2" etc to days
msdat.ps.peaks <- msdat.ps.peaks %>%
  mutate(time_days  = case_when(PyOM  == "350" & Timepoint == "T2" ~ 2,
                                PyOM  == "350" & Timepoint == "T3" ~ 7,
                                PyOM  == "350" & Timepoint == "T4" ~ 18,
                                PyOM  == "550" & Timepoint == "T3" ~ 4,
                                PyOM  == "550" & Timepoint == "T4" ~ 10,
                                PyOM  == "550" & Timepoint == "T5" ~ 22))

# msdat.ps.peaks <- msdat.ps.melt %>%
#   filter(OTU %in% features) %>%
#   filter(Treatment != "S") %>%
#   group_by(OTU, Treatment, Timepoint, Sample)%>%
#   summarize(Abundance = sum(Abundance),
#             sd_abund = sd(Abundance),
#             n_abund = n(), 
#             se_abund = sd_abund/sqrt(n_abund),
#             lower_se_abund = mean_abund - se_abund,
#             upper_se_abund = mean_abund + se_abund)


#Plotting Fig. S4
# Create a named vector with the new facet labels
unique_peaks <- unique(msdat.ps.peaks$OTU)
peak_names <- setNames(paste("Feature", seq_along(unique_peaks)), unique_peaks)

# Add a new column with the new labels
msdat.ps.peaks$peak_label <- as.factor(peak_names[as.character(msdat.ps.peaks$OTU)])

# Plot
p = ggplot(msdat.ps.peaks,
           aes(x = as.factor(time_days), y = Abundance , color = as.factor(Treatment)))
p = p+scale_color_manual(name = "Treatment", values = c("black", "red"), labels = c("soil + 350 PyOM", "soil + 550 PyOM"))
p = p+theme_bw(base_size = 18)
p = p+theme(axis.ticks = element_line(),
            axis.text = element_text(colour = "black"),
            axis.title = element_text(face = "bold"),
            plot.title = element_blank(),
            legend.position = "top",
            legend.title=element_blank(), 
            legend.text=element_text(),
            strip.text = element_text())
p = p+geom_boxplot()
p = p+labs(x = "Incubation time (days)", y = "Relative abundance\n")
p = p+facet_wrap(~ peak_label, scales = "free")
#p = p+facet_wrap(~ OTU, scales = "free", labeller = as_labeller(peak_names)) ##use this to change facet labels
p = p + scale_y_continuous(labels = scales::scientific)
p
#ggsave("./pyom_peaks2.png", dpi = 300, device = "png", width = 8, height = 7)



 