#Load necessary packages
library (tidyverse)
library(wesanderson)

######Reading in dataframe ######
df.cml <- read.csv("./cml_respired.csv", row.names = 1, header = TRUE, as.is = TRUE, check.names = FALSE)
df.cml$temperature <- as.factor(df.cml$temperature)

#NOTE: Change this data frame to "rate_respired.csv" to plot time series with rates of C mineralized

#Creating standard error functions
lower_se <- function(mean, se){
  res <- mean - se
  return(res)
}
upper_se <- function(mean, se){
  res <- mean + se
  return(res)
}

df.c.means <- df.cml %>%
  group_by(treatment, temperature, cycle) %>%
  dplyr::summarize(
    mean_C_mg_cml = mean(C_mg_norm_cml),
    sd_C_mg_cml = sd(C_mg_norm_cml),
    mean_time_days = mean(time_days),
    n_C_mg_cml= n()) %>%
  dplyr::mutate(se_C_mg_cml= sd_C_mg_cml/sqrt(n_C_mg_cml),
                lower_se_C_mg_cml=lower_se(mean_C_mg_cml, se_C_mg_cml),
                upper_se_C_mg_cml=upper_se(mean_C_mg_cml, se_C_mg_cml))

######### FOR cumulative C figure ######################
cml_fig_trs <- c("water-extractable PyOM-C","SOC",
"non-water-extractable PyOM-C")

df.c.means.trt <- df.c.means %>%
  dplyr::filter(treatment %in% cml_fig_trs)

df.c.means.trt$treatment <- relevel(as.factor(df.c.means.trt$treatment), "water-extractable PyOM-C")

#Plotting cumulative C mineralized for all fractions (Fig. 2A and Fig. S1)
p = ggplot(df.c.means.trt)
# Add data
p = p + geom_errorbar(aes(x=mean_time_days, ymin=lower_se_C_mg_cml,ymax=upper_se_C_mg_cml,
                          colour = temperature), width = 0.3, size = 1)
p = p + geom_point(aes(x=mean_time_days,y=mean_C_mg_cml, shape=treatment, colour = temperature), size=6, stroke = 1.2)
p = p + geom_line(aes(x=mean_time_days,y=mean_C_mg_cml, colour = temperature, alpha = treatment), size = 1)
# Formatting
p = p + theme_bw(base_size = 24)
p = p + geom_abline(slope=0,intercept=0)
p = p + xlab("Incubation time (days)") + ylab(expression(atop(bold(Cumulative~C~mineralized), 
                                                              ~(mg~CO[2]-C))))
p = p + theme(legend.title = element_blank(),legend.text=element_text(), legend.position = "top")
p = p + theme(axis.line = element_line(colour = "black"), axis.text=element_text(),axis.title = element_text(face = "bold"))
p = p + theme(strip.text.x = element_text())
p = p + scale_colour_manual(values = wes_palette("Royal1"),labels = c("350 PyOM", "550 PyOM" ))
#p = p + scale_colour_manual(values = wes_palette("Royal1")[2],labels = c("550 PyOM"))
p = p + scale_shape_manual(values=c(1, 2, 0))
p = p + scale_alpha_manual(values=c(1, 0.9999, 0.9998))
# p = p + geom_vline(xintercept = c(0,2,7,18,26), linetype="dashed", 
#                 color = "black", size=0.5)
# p = p + geom_vline(xintercept = c(0,4,10,22), linetype="dashed",
#                 color = "black", size=0.5)
p
#ggsave("./plots/partitioned_cml_cor2_noleg.png", dpi = 300, device="png", width=8, height=8)

###For priming figure###########
##Separating soil means for PyOM amended and unamended soils
df.soil.c.means <- df.c.means %>%
  filter(treatment == "SOC")
df.soil.control.c.means <- df.c.means %>%
  filter(treatment == "SOC(no PyOM)")
df.soil.c.means['treatment'] = NULL
df.soil.control.c.means['treatment'] = NULL
df.soil.c.means.cols <- c("temperature", "cycle", "mean_soilC_mg_cml",
                       "sd_soilC_mg_cml", "mean_time_days", 
                       "n_soilC_mg_cml","se_soilC_mg_cml",
                       "lower_se_soilC_mg_cml", "upper_se_soilC_mg_cml")
colnames(df.soil.c.means) <- df.soil.c.means.cols

df.prim.c <- merge(df.soil.c.means, df.soil.control.c.means, by.x = c("temperature", "cycle"), by.y = c("temperature","cycle"))
df.prim.c$temperature <- as.factor(df.prim.c$temperature)
colnames(df.prim.c)
#Calculating cumulative priming and error propagation
df.prim.c = df.prim.c %>%
  dplyr::mutate(soilC_prim_mg_cml = mean_soilC_mg_cml - mean_C_mg_cml) %>%
  dplyr::mutate(error_prim_cml = sqrt((sd_soilC_mg_cml)^2 + (sd_C_mg_cml)^2)) %>%
  dplyr::mutate(frac_soilC_prim_cml = (soilC_prim_mg_cml/mean_C_mg_cml)) %>%
  dplyr::mutate(error_prim_frac = (sqrt((error_prim_cml/soilC_prim_mg_cml)^2 + 
                  (sd_C_mg_cml/mean_C_mg_cml)^2))*frac_soilC_prim_cml) %>%
  dplyr::mutate(perc_soilC_prim_cml = frac_soilC_prim_cml*100) %>%
  dplyr::mutate(error_prim_perc = error_prim_frac*100)
  
ggplot(df.prim.c) + geom_point(aes(x = mean_time_days.x, y = perc_soilC_prim_cml)) + facet_wrap(~temperature)

#Creating confidence interval functions 
# lower_ci <- function(mean, se, zvalue = 1.960){
#   res <- mean - (zvalue * se)
#   return(res)
# }
# upper_ci <- function(mean, se, zvalue = 1.960){
#   res <- mean + (zvalue * se)
#   return(res)
# }

df.prim.c <- df.prim.c %>%
  dplyr::mutate(se_prim = error_prim_perc/sqrt(n_soilC_mg_cml),
                lower_se_prim = lower_se(perc_soilC_prim_cml, se_prim),
                upper_se_prim = upper_se(perc_soilC_prim_cml, se_prim))
#Plotting Fig. 2B
p = ggplot(df.prim.c)
p = p + geom_point(aes(x = mean_time_days.x, y = perc_soilC_prim_cml, shape = temperature, color = temperature), size = 6, stroke  = 1.2)
p = p + geom_errorbar(aes(x = mean_time_days.x,
                          ymin = lower_se_prim,
                          ymax = upper_se_prim,
                          linetype = temperature), width = 0.2)
p = p + geom_line(aes(x = mean_time_days.x, y = perc_soilC_prim_cml, linetype = temperature), size = 1.4)
#Formatting
p = p + scale_shape_manual(values = c(15,0), labels = c("soil + 350 PyOM", "soil + 550 PyOM"))
p = p + scale_linetype_manual(values = c("solid", "dashed"), labels = c("soil + 350 PyOM", "soil + 550 PyOM"))
p = p + scale_colour_manual(values = c("black", "grey54"), labels = c("soil + 350 PyOM", "soil + 550 PyOM"))
p = p + geom_abline(slope=0,intercept=0, linetype = "dotted")
p = p + xlab("Incubation time (days)") + ylab(expression(atop(bold(Change~"in"~SOC~mineralization),                                                               ~("%"~of~no-PyOM~treatment))))
#p = p + ylim(-20,100)
p = p + theme_bw(base_size = 24)
p = p + theme(legend.title = element_blank(),legend.text= element_text(size = 24), legend.position = "top")
p = p + theme(axis.text=element_text(size=24),axis.title = element_text(size=26, face = "bold"), axis.line = element_line(colour = "black"))
p

#ggsave("./plots/soil_priming_cml_leg_2.png", dpi = 300, device="png", width=8, height=8)


#####################################
###Plotting priming mechanism correlations
###Reading dataset
df.rate <- read.csv("./rate_respired.csv", row.names = 1, header = TRUE)
#Calculating the means and standard errors for C values at each time point 
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
    mean_C_mg_per_day = mean(CO2C_mg_per_day),
    sd_C_mg_per_day = sd(CO2C_mg_per_day),
    mean_time_days = mean(time_days),
    n_C_mg_per_day= n()) %>%
  dplyr::mutate(se_C_mg_per_day= sd_C_mg_per_day/sqrt(n_C_mg_per_day),
                lower_se_C_mg_per_day=lower_se(mean_C_mg_per_day, se_C_mg_per_day),
                upper_se_C_mg_per_day=upper_se(mean_C_mg_per_day, se_C_mg_per_day))
##Separating soil means for PyOM amended and unamended treatments
df.soil.means <- df.means %>%
  filter(treatment == "SOC")
df.soil.control.means <- df.means %>%
  filter(treatment == "SOC(no PyOM)")
df.soil.means['treatment'] = NULL
df.soil.control.means['treatment'] = NULL
df.soil.means.cols <- c("temperature", "cycle", "mean_soilC_mg_per_day",
                       "sd_soilC_mg_per_day", "mean_time_days", 
                      "n_soilC_mg_per_day","se_soilC_mg_per_day",
                      "lower_se_soilC_mg_per_day", "upper_se_soilC_mg_per_day")
colnames(df.soil.means) <- df.soil.means.cols

###Calculating priming
df.prim <- merge(df.soil.means, df.soil.control.means, by.x = c("temperature", "cycle"), by.y = c("temperature","cycle"))
df.prim$temperature <- as.factor(df.prim$temperature)

df.prim <- df.prim %>%
  mutate(soilC_prim_mg = mean_soilC_mg_per_day - 
           mean_C_mg_per_day) %>%
  mutate(error_prim = sqrt((sd_soilC_mg_per_day)^2 + 
                             (sd_C_mg_per_day)^2)) %>%
  mutate(frac_soilC_prim = soilC_prim_mg/mean_C_mg_per_day) %>%
  mutate(error_prim_frac = (sqrt((error_prim/soilC_prim_mg)^2 + 
                                          (sd_C_mg_per_day/mean_C_mg_per_day)^2))*frac_soilC_prim) %>%
  mutate(perc_soilC_prim = frac_soilC_prim*100) %>%
  mutate(error_prim_perc = error_prim_frac*100)

df.prim <- df.prim %>%
  dplyr::mutate(se_prim = error_prim_perc/sqrt(n_soilC_mg_per_day),
                lower_se_prim = lower_se(perc_soilC_prim, se_prim),
                upper_se_prim = upper_se(perc_soilC_prim, se_prim))
  
  
  

ggplot(df.prim) + geom_point(aes(x = mean_time_days.x, y = perc_soilC_prim)) + facet_wrap(~temperature)


df.prim.rows <- c("temperature", "cycle", "perc_soilC_prim", "mean_time_days.x", "lower_se_prim","upper_se_prim")
df.prim.trimmed <- df.prim[df.prim.rows]

##Separating df.wx means dataset
df.wx.means <- df.means %>%
  filter(treatment == "water-extractable PyOM-C")
df.wx.means['treatment'] = NULL

df.wx.means.cols <- c("temperature", "cycle", "mean_wxC_mg_per_day",
                       "sd_wxC_mg_per_day", "mean_time_days", 
                       "n_wxC_mg_per_day","se_wxC_mg_per_day",
                       "lower_se_wxC_mg_per_day", "upper_se_wxC_mg_per_day")
colnames(df.wx.means) <- df.wx.means.cols

df.wx.rows <- c("temperature", "cycle", "mean_wxC_mg_per_day", "sd_wxC_mg_per_day", "n_wxC_mg_per_day")
df.wx.trimmed <- df.wx.means[df.wx.rows]

###Separating df.nonwx means
df.nonwx.means <- df.means %>%
  filter(treatment == "non-water-extractable PyOM-C")
df.nonwx.means['treatment'] = NULL

df.nonwx.means.cols <- c("temperature", "cycle", "mean_nonwxC_mg_per_day",
                     "sd_nonwxC_mg_per_day", "mean_time_days", 
                     "n_nonwxC_mg_per_day","se_nonwxC_mg_per_day",
                     "lower_se_nonwxC_mg_per_day", "upper_se_nonwxC_mg_per_day")
colnames(df.nonwx.means) <- df.nonwx.means.cols

df.nonwx.rows <- c("temperature", "cycle", "mean_nonwxC_mg_per_day", "sd_nonwxC_mg_per_day", "n_nonwxC_mg_per_day")
df.nonwx.trimmed <- df.nonwx.means[df.nonwx.rows]

#Merging dataframes for plotting
df.merged <- merge(df.prim.trimmed, df.wx.trimmed, by=c("temperature", "cycle"))
df.merged <- merge(df.merged, df.nonwx.trimmed, by=c("temperature", "cycle"))
#Adding up the C in different PyOM fractions and calculating errors
df.merged = df.merged %>%
  dplyr::filter(temperature==350) %>% #change to 550 for Fig. S3 
  dplyr::mutate(sum = mean_nonwxC_mg_per_day + mean_wxC_mg_per_day) %>%
  dplyr::mutate(error_pyc = sqrt((sd_nonwxC_mg_per_day)^2 + 
                               (sd_wxC_mg_per_day)^2))  %>%
  dplyr::mutate(se_pyc = error_pyc/sqrt(n_wxC_mg_per_day),
                 lower_se_pyc = lower_se(sum, se_pyc),
                 upper_se_pyc = upper_se(sum, se_pyc)) %>%
  dplyr::filter(mean_time_days.x <=2) # for first 48 hours 
  #dplyr::filter(mean_time_days.x > 2) %>% # for post-48 hours 
  #dplyr::filter(cycle != 37) #cycle !=37 for 350
  

colnames(df.merged)
mod_eqn = function(df){
  m1 = lm(perc_soilC_prim ~ sum, df)
  r2 = summary(m1)$r.squared
  m2 = aov(perc_soilC_prim ~ sum, df)
  p_val = summary(m2)[[1]][["Pr(>F)"]][[1]]
  if (p_val < 0.001) {
    eq <- substitute({italic(R)^{2}}~"="~r2~","~{italic(p)}~"<"~0.001,
                     list(r2 = format(r2, digits = 2)))
  } else {
    eq <- substitute({italic(R)^{2}}~"="~r2~","~{italic(p)}~"="~p_val,
                     list(r2 = format(r2, digits = 2), p_val = format(p_val, digits = 1)))
  }
  as.character(as.expression(eq))
}


#Plotting Fig. 6A and 6B
p = ggplot(df.merged)
p = p + geom_point(aes(x = sum, y = perc_soilC_prim), size = 4) #colour = as.numeric(mean_time_days.x)
# p = p + geom_smooth(method = lm, level = FALSE, colour = "black", 
#                     mapping = aes(x = sum, y = perc_soilC_prim))
p = p + geom_text(aes(label = mod_eqn(df.merged), x=Inf, y=Inf), parse=TRUE, hjust=1.3, vjust=1, size = 7)
p = p + geom_errorbarh(aes(xmax = upper_se_pyc, xmin = lower_se_pyc, y = perc_soilC_prim),
                       color = "gray72")
p = p + geom_errorbar(aes(ymin = lower_se_prim, ymax = upper_se_prim, x = sum),
                      color = "gray72")
#Formatting
p = p + theme_bw(base_size = 24)
p = p + theme(axis.text.x=element_text(), axis.text.y=element_text(), 
              axis.title = element_text(face = "bold"),
              axis.line = element_line(colour = "black"))
p = p + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
p = p + xlab(expression(atop(bold(PyOM-C~mineralization~per~day),
                             ~(mg~CO[2]-C~g^"-1"~added~C~day^"-1")))) 
p = p + ylab(expression(atop(bold(SOC~priming~per~day), 
                             ~("%"~of~no-PyOM~treatment))))
p = p + guides(colour = guide_legend(title="Time elapsed (days)"))
p
#ggsave("./plots/prim_vs_py350_48hr.png", dpi = 300, device="png", width=8, height=8)
#ggsave("./plots/soil_priming_550.png", dpi = 300, device="png", width=8, height=8)
