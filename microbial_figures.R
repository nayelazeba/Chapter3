source("nz_ps_lib.R")
# Read the ps object
#filename = "./ps.doe.043022.ITS" #For ITS dataset
filename = "./ps.doe.5.16S" #For 16S dataset
ps = get_ps(filename)

# Removing chloroplasts and mitochondria for bacteria
ps <- ps %>%
  subset_taxa(
    Domain == "Bacteria"&
      Family != "Mitochondria" &
      Class != "Chloroplast"
  )

# Filter the phyloseq object to keep only the blank samples
blank_samples <- subset_samples(ps, Sample_type == "blank")

# Calculate the mean number of reads in the blank samples
max_reads_blank <- max(sample_sums(blank_samples))
max_reads_blank

# Pruning to remove blank and low abundance samples
# For bacteria
ps = prune_ps(ps,"bacteria")
ps = prune_samples(sample_sums(ps) > max_reads_blank, ps) #cut off based on the the mean no. of reads observed in the blank samples
#ps_cut = prune_samples(sample_sums(ps)< max_reads_blank, ps)

# For fungi
ps = prune_ps(ps,"fungi")
ps = prune_samples(sample_sums(ps) > max_reads_blank, ps) #cut off based on the the mean no. of reads observed in the blank samples

# Subsetting 350 samples
ps_350 = subset_samples(ps, Temperature == 350)
# Converting sample_type variable to factor and setting the reference level
sample_data(ps_350)$Sample_type <- factor(sample_data(ps_350)$Sample_type, levels = c("control","sample"))
sample_data(ps_350)$Sample_type <- relevel(sample_data(ps_350)$Sample_type, ref = "control")
# Converting time variable to factor and setting the reference level
sample_data(ps_350)$Timepoint <- factor(sample_data(ps_350)$Timepoint, levels = c("t1","t2", "t3", "t4", "t5"))
sample_data(ps_350)$Timepoint <- relevel(sample_data(ps_350)$Timepoint, ref = "t1")
# Normalizing
ps_350_relative <- normalize_ps(ps_350, "relative")

# Subsetting 550 samples
ps_550 = subset_samples(ps, Temperature == 550 & Timepoint!= "t1")
# Converting sample_type variable to factor and setting the reference level
sample_data(ps_550)$Sample_type <- factor(sample_data(ps_550)$Sample_type, levels = c("control","sample"))
sample_data(ps_550)$Sample_type <- relevel(sample_data(ps_550)$Sample_type, ref = "control")
# Converting time variable to factor and setting the reference level
sample_data(ps_550)$Timepoint <- factor(sample_data(ps_550)$Timepoint, levels = c("t2", "t3", "t4", "t5"))
sample_data(ps_550)$Timepoint <- relevel(sample_data(ps_550)$Timepoint, ref = "t2")
# Normalizing
ps_550_relative <- normalize_ps(ps_550, "relative")

# Plotting NMDS figures (all panels in Fig. 3)
source("nz_ps_lib.R")
ps_nmds = get_nmds_vegan(ps_550_relative) 
plot = plot_nmds(ps_nmds, ps_550_relative, 550, "bacterial", "Sample_type", "Timepoint", hjust=-0.8, vjust=1.3)
plot  = plot + theme(legend.position = "none")
plot
#ggsave("./plots/nmds_fungi_550.png", dpi=300, device="png", width=4, height=4)

################ PERMANOVA ##############################
run_adonis2 <- function(ps) {
  s = data.frame(sample_data(ps))
  psotu2veg <- function(ps) {
    OTU <- otu_table(ps)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  ps.otu.veg = psotu2veg(ps)
  ps_bray <- vegdist(ps.otu.veg, method  = "bray")
  adonis2(ps_bray ~ Timepoint+Sample_type, data = s)
}
run_adonis2(ps_350_relative)

#####################DESeq2####################################
# Checking out the sample data for the phyloseq object
colnames(sample_data(ps_350))

#Modeling the difference between sample and control treatments at time 0, 
#the difference over time, and any sample type specific differences over time
formula = as.formula("~Sample_type + Timepoint + Sample_type:Timepoint") 

# Convert to DeSeq2 object
diagdds = phyloseq_to_deseq2(ps_350, formula)

#Function to calculate geometric means to include 0 counts 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

##https://support.bioconductor.org/p/76548/
##https://support.bioconductor.org/p/62246/#62250
##https://github.com/joey711/phyloseq/issues/283

#Calculate geometric means
geoMeans = apply(counts(diagdds), 1, gm_mean)

#Estimate scaling factor
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)

#Running deseq
#Performing likelihood ratio test
#Using LRT with reduced formula 
#OTUs with small p values from this test are those which at one or more 
#time points after time 0 showed a sample type specific effect. 
#Does not include OTUS that moved up or down over time in the same way in both sample and control.
diagdds = DESeq(diagdds, fitType="local", quiet = FALSE, test="LRT", reduced = ~ Timepoint + Sample_type)

#Extracting results 
#List coefficients
resultsNames(diagdds)
res <- results(diagdds)
#Extracting results at a given time point
#res5 <- results(diagdds, name="Sample_typesample.Timepointt5", test="Wald")
#View results
head(res)
summary(res)
#res <- res5
# Adding in taxonomy data to deseq results.
restab = cbind(as(res, "data.frame"), as(tax_table(ps_350)[rownames(res), ], "matrix"))
restab = rownames_to_column(restab, var= "OTU")

#Filtering to remove NA padj rows
restab = restab[!is.na(restab$padj),]

#Selecting significant responders
alpha = 0.05
sigtab = restab[(restab$padj < alpha), ]

#Setting levels for padj values
#restab$sig = ifelse(restab$padj<0.05,1,0)
#restab$sig = as.factor(restab$sig)

# Creating a fasta file of the incorporators for BLASTing
#library(stringr)
# read in fasta file and extract sequences
# fasta <- readLines("./qiime_output_16S/dna-sequences.fasta")
# otu_ids <- character()
# otu_seqs <- character()
# 
# for (i in 1:length(fasta)) {
#   if (substr(fasta[i], 1, 1) == ">") {
#     otu_ids <- c(otu_ids, substr(fasta[i], 2, nchar(fasta[i])))
#     otu_seqs <- c(otu_seqs, fasta[i+1])
#   }
# }
# create named list of sequences
# fasta_list <- setNames(otu_seqs, otu_ids)

# select OTUs
# otus_to_select <- sigtab$OTU
# selected_seqs <- fasta_list[otus_to_select]
# 
# # write extracted sequences to fasta file
# writeLines(paste0(">", names(selected_seqs), "\n", selected_seqs), "350resp_otu_sequences.fasta")

#For making the figures
#Filtering out low counts 
abund_cutoff <- quantile(restab$baseMean, 0.25)
sigtab <- sigtab %>%
  filter(baseMean >= abund_cutoff)

# Replace genus names with "Unknown" if they contain "unknown" or "uncultured"
sigtab$Genus <- ifelse(grepl("(unknown|uncultured)",sigtab$Genus, ignore.case = TRUE), "Unknown", sigtab$Genus)
sigtab$Genus <- ifelse(grepl("(unknown|uncultured)",sigtab$Genus, ignore.case = TRUE), "Unknown", sigtab$Genus)
# Replace NA values in the "genus" column with "Unknown"
sigtab$Genus[is.na(sigtab$Genus)] <- "Unknown"
sigtab$Genus <- ifelse(sigtab$Genus == "", "Unknown", sigtab$Genus)
#write.csv(sigtab, file="bact_sigpos_350pyom_t5.csv")

#Plotting Fig. 5A and S2
options(repr.plot.width=4, repr.plot.height=1.5)
p = ggplot(sigtab,aes(x=Genus, y=log2FoldChange, fill=Phylum))
p = p + scale_fill_manual(values=c(brewer.pal(length(unique(sigtab$Phylum)),"Spectral")))
p = p + geom_jitter(shape=21, stroke = 1, width = 0.2,  size = 4)#stroke = 1, width = 0.2,  size = 4
p = p + scale_size(guide=FALSE)
#p = p + scale_alpha_manual(guide=FALSE, values=c(0.2,0.8)) #if setting alpha for p-values
p = p + theme_bw()
p = p + geom_hline(yintercept=0)
#p = p + ylim(c(0,1.5))
p = p + ylab(expression(atop(bold(Response~to~"350"~PyOM~addition),
                             ~(log[2]-fold~change~italic(vs)~soil))))
p = p + theme(axis.text.x = element_text(angle=65, hjust=1, vjust = 1, size=16, face="italic"),
              axis.text.y = element_text(size = 16),
              legend.position="top",
              legend.title = element_blank(),
              legend.text = element_text(size=16,face="italic"),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size=16))
p = p + guides(fill=guide_legend(title="Phylum",override.aes = list(size=4)),alpha=FALSE)
p
#ggsave("./manuscript_fig/deseq_bac350.png", dpi=300, device="png", width=9.5, height=4.7)

#####To make abundance figures  #############
#ps_550_relative = normalize_ps(ps_550, "relative")

#Covert to dataframe
ps_tab <- ps_350_relative %>%
  psmelt ()   

#Selecting OTUs from the deseq2 output that responded positively 
#Selecting positive responders
sigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
resp_otu <- paste(unique(sigtab$OTU))
ps_otutab <- ps_tab %>% 
  filter(OTU %in% resp_otu) %>%
  filter(!(Genus == "Unknown")) %>%
  filter(!(Genus == "uncultured")) %>%
  group_by(Genus, Family, Phylum, Sample, Timepoint, Sample_type) %>%
  summarise(Abundance=sum(Abundance),
            n = n())

# #Cleaning the genus name (if including family)
# ps_gentab = ps_otutab %>%
#   filter(!(Family == "")) %>%
#   mutate(Genus = ifelse(is.na(Genus) | grepl("uncultured", Genus, ignore.case = TRUE),
#                         paste("Unk.", Family),
#                         paste(Genus)))
#write.csv(ps_350_gentab, "resp_genera_350.csv")

#Extracting the most abundant phyla
RespOrder = ps_gentab %>%
  arrange(-n)
RespOrder = paste(RespOrder$Genus)
RespOrder = unique(RespOrder)
nGen = 6
ps_gentab = ps_gentab[ps_gentab$Genus %in% RespOrder[1:nGen],]

#Selecting bacterial genera for plotting for 350PyOM
# plotted_genera350 <- c("Noviherbaspirillum", "Gemmatimonas", "Unk. Comamonadaceae", "Pseudonocardia", "Saccharimonadales", "Massilia")
# ps_gentab <- ps_gentab[ps_gentab$Genus %in% plotted_genera350,]

#Plotting Fig. 5B and S6
p = ggplot(ps_gentab,
           aes(x = Timepoint, y = Abundance , color = Sample_type))
p = p+scale_color_manual(name = "Treatment", values = c("darkorange3", "black"), labels = c("soil", "soil + 350 PyOM"))
p = p+theme_bw(base_size = 18)
p = p+theme(axis.ticks = element_line(),
            axis.text = element_text(colour = "black"),
            axis.title = element_text(face = "bold"),
            plot.title = element_blank(),
            legend.position = "top",
            legend.title=element_blank(), 
            legend.text=element_text(),
            strip.text = element_text(face = "italic"))
p = p+scale_x_discrete(labels=c("t1" = "0", "t2" = "2",
                                "t3" = "7", "t4" = "18", "t5" = "26"))
#For 550 PyOM treatment
# p = p+scale_x_discrete(labels=c("t2" = "0",
#                                 "t3" = "4", "t4" = "10", "t5" = "22"))
p = p+geom_boxplot()
p = p+labs(x = "Incubation time (days)", y = "Relative abundance\n")
p = p+facet_wrap(~ Genus, scales = "free")
p = p+stat_compare_means(aes(group = Sample_type), label = "p.signif", hide.ns = TRUE, show.legend = FALSE, size = 6) 
p = p+scale_y_continuous(expand = c(0.3, 0), labels = label_number(accuracy = 0.001))
p
#ggsave("./manuscript_fig/bact350_Gen6B.png", dpi = 300, device = "png", width = 9.5, height = 5.5)
#ggsave("./manuscript_fig/bac550_gen.png", dpi = 300, device = "png", width = 4, height = 4)




