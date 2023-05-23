#Plotting supplementary figure S5 showing change in relative abundance adn C data simultaneously
source("../sequencing/nz_ps_lib.R")
library(cowplot)
library(egg)

#Read in the C mineralization data
df.rate <-read.csv("../picarro_data_350/rate_respired.csv", row.names = 1)
df.rate$temperature <- as.factor(df.rate$temperature)

#Creating standard error functions
lower_se <- function(mean, se){
  res <- mean - se
  return(res)
}
upper_se <- function(mean, se){
  res <- mean + se
  return(res)
}


df.means <- df.rate %>%
  group_by(treatment, temperature, cycle) %>%
  dplyr::summarize(
    mean_C_rate = mean(CO2C_mg_per_day),
    sd_C_rate = sd(CO2C_mg_per_day),
    time_days = mean(time_days),
    n_C_rate= n()) %>%
  dplyr::mutate(se_C_rate= sd_C_rate/sqrt(n_C_rate),
                lower_se_C_rate=lower_se(mean_C_rate, se_C_rate),
                upper_se_C_rate=upper_se(mean_C_rate, se_C_rate))

#Filter treatments
trt <- c("water-extractable PyOM-C", "non-water-extractable PyOM-C")

df.means <- df.means %>%
  dplyr::filter(treatment %in% trt) %>%
  dplyr::filter(temperature == 350)


df.means$treatment = relevel(as.factor(df.means$treatment), "water-extractable PyOM-C")

##Read the rel. abund. data for 350 PyOM responders
ps_350_gentab <- read.csv("../sequencing/resp_genera_350.csv", row.names = 1)
#Renaming time points
ps_350_gentab <- ps_350_gentab %>%
  mutate(time_days  = case_when(Timepoint == "t1" ~ 0,
                                Timepoint == "t2" ~ 2,
                                Timepoint == "t3" ~ 7,
                                Timepoint == "t4" ~ 18,
                                Timepoint == "t5" ~ 26))
#Filtering data to select significant responders
gemma <- ps_350_gentab %>%
  filter(Genus == "Gemmatimonas")

pseudo <- ps_350_gentab %>%
  filter(Genus == "Pseudonocardia")

#Plotting separate plots and combining them
p1 = ggplot(df.means, aes(x = as.factor(time_days), y = mean_C_rate)) +
  geom_point(aes(shape = treatment), size=3, stroke = 1.2) +
  geom_line(aes(group=treatment, linetype = treatment), size = 1) +
  geom_abline(slope=0,intercept=0) +
  geom_errorbar(aes(x= as.factor(time_days), ymin=lower_se_C_rate,ymax=upper_se_C_rate,), width = 0.3, size = 1) +
  scale_y_continuous(limits = c(-1, max(df.means$upper_se_C_rate, na.rm = TRUE)), 
                     expand = c(0,0)) +
  labs(y = expression(atop(bold(C~mineralized~per~day),
                           ~(mg~CO[2]-C~g^"-1"~added~C~day^"-1")))) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        panel.grid.major = element_line(),
        panel.grid.minor = element_line(),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  scale_shape_manual(values=c(1, 2)) +
  scale_linetype_manual(values  = c("solid", "dashed"))

p2 = ggplot(gemma, aes(x = as.factor(time_days), y = Abundance)) +
  geom_boxplot(aes(color = Sample_type)) +
  scale_y_continuous(limits = c(min(gemma$Abundance, na.rm = TRUE), max(gemma$Abundance, na.rm = TRUE)), 
                     expand = c(0.3, 0), 
                     labels = label_number(accuracy = 0.001)) +
  labs(x = "Time (days)", y = "Relative\n Abundance") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(),
        axis.title.y = element_text(face = "bold", size = 18),
        axis.text.y = element_text(size = 18),
        panel.grid.major = element_line(),
        panel.grid.minor = element_line(),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_color_manual(name = "Treatment", values = c("darkorange3", "black"), labels = c("soil", "soil + 350 PyOM")) +
  stat_compare_means(aes(group = Sample_type), label = "p.signif", hide.ns = TRUE, show.legend = FALSE, size = 6) 
                     
p3 = ggplot(pseudo, aes(x = as.factor(time_days), y = Abundance)) +
  geom_boxplot(aes(color = Sample_type)) +
  scale_y_continuous(limits = c(min(pseudo$Abundance, na.rm = TRUE), max(pseudo$Abundance, na.rm = TRUE)), 
                     expand = c(0.3, 0), 
                     labels = label_number(accuracy = 0.001)) +
  labs(x = "Time (days)", y = "Relative\n Abundance") +
  theme(axis.line.y = element_line(),
        axis.line.x = element_line(),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 18),
        panel.grid.major = element_line(),
        panel.grid.minor = element_line(),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_color_manual(name = "Treatment", values = c("darkorange3", "black"), labels = c("soil", "soil + 350 PyOM")) +
  stat_compare_means(aes(group = Sample_type), label = "p.signif", hide.ns = TRUE, show.legend = FALSE, size = 6) 



cowplot::plot_grid(p1, p2, p3, align = "v", ncol = 1, rel_heights = c(0.43, 0.28, 0.28))
egg::ggarrange(p1, p2, p3, heights = c(0.43, 0.28, 0.28))

#ggsave("./supp_fig.png", device = "png", dpi = 300, width  = 8, height = 9)
