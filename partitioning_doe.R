#Load necessary packages
library (tidyverse)

#Parameters specific to incubation experimental design:
jarvol_L <- 0.473 # jar volume (L) for 16oz pint size Mason jars
pyom_350_13C_d13C <- 1513.12 #isotopic composition of 13C labeled 350 PyOM (d13C per mille vPDB)
pyom_550_13C_d13C <- 1596.51 #isotopic composition of 13C labeled 550 PyOM (d13C per mille vPDB)

total_wxC_g_350 = 0.0001 #mass of water-extractable 350 PyOM-C in each jar (g)
total_wxC_g_550 = 0.00007 #mass of water-extractable 550 PyOM-C in each jar (g)
total_nonwxC_g_350 = 0.06 #mass of non-water-extractable 350 PyOM-C in each jar (g)
total_nonwxC_g_550 = 0.06 #mass of non-water-extractable 550 PyOM-C in each jar (g)
total_soilC_g = 0.32 # mass of SOC in each jar (g)
total_pyomC_g_350 = total_nonwxC_g_350+ total_wxC_g_350 #mass of 350 PyOm-C in each jar (g)
total_pyomC_g_550 = total_nonwxC_g_550+ total_wxC_g_550 #mass of 550 PyOm-C in each jar (g)

#############################################################################
##############################################################################
#Reading data tables
df.flux <- read.csv("CO2Flux.csv", row.names = 1)
df.flux$temperature <- as.factor(df.flux$temperature)
df.flux$treatment <- as.factor(df.flux$treatment)

#Covert CO2 values from ppm to a mass basis. We make some simple assumptions that the 
#sample is at standard temperature and pressure. We also adjust time into time since 
#incubation started and convenient units
df.min <- df.flux %>%
  dplyr::group_by(sample_num) %>%
  dplyr::mutate(time_s = time-min(time)) %>%
  dplyr::mutate(time_hr = (time_s/3600)) %>%
  dplyr::mutate(time_days = time_hr/24) %>%
  dplyr::arrange(time_days) %>%
  dplyr::mutate(CO2C_mg = sample_CO2/1000000/22.4*jarvol_L*12.01*1000) %>% 
  dplyr::mutate(time_int = time_days - lag(time_days)) %>% #calculating the time interval between each measurement
  dplyr::filter(!is.na(sample_CO2)) %>%
  dplyr::filter(time_days<30) %>% #keeping data collected within 30 days of incubation
  dplyr::select(-time_s, -time_hr, -time)

#Treatment information for both 350 and 550 PyOM (Refer to Figure 1 in MS)
# 1. 12C_bk_soil: Soil + unlabeled PyOM
# 2. 13C_wx_soil: Soil + 13C labeled water-extractable PyOM
# 3. 13C_nox_soil: Soil + 13C labeled non-water-extractable PyOM
# 4. 13C_bk_soil: Soil + 13C labeled PyOM
# 5. soil
# 6. blank

#Using blank treatments to correct CO2 values in the dataset at each cycle
df.blank.means1 = df.min %>%
  filter(treatment == "blank") %>%
  group_by(cycle) %>%
  dplyr::summarize(mean_blank_CO2C_mg = mean(CO2C_mg))

mineralization = df.min %>%
  filter (treatment!= "blank") %>%
  dplyr::left_join(df.blank.means1) %>%
  mutate(CO2C_mg_corr = CO2C_mg - mean_blank_CO2C_mg)

#################################################################################
# Plotting C mineralized from blank treatments
#############################################################################
df.blank = df.min %>%
  filter(treatment == "blank" & sample_num!=2) %>% ###values from sample 2 are very high compared to the other three
  mutate(CO2C_mg_per_day = CO2C_mg/time_int)

ggplot(df.blank) + geom_point(aes(x = time_days, y = CO2C_mg_per_day))+ facet_wrap(~temperature)

## Cumulative
df.blank.c <- df.blank %>%
  dplyr::group_by(sample_num)%>%
  dplyr::arrange(cycle)%>%
  dplyr::mutate(CO2_C_mg_cml= cumsum(CO2C_mg))

ggplot(df.blank.c) + geom_point(aes(x = time_days, y = CO2_C_mg_cml)) + facet_wrap(~temperature)

###############################################################################
#Partitioning water-extractable PyOM-C mineralized
################################################################################
#We can determine the fraction of CO2 emitted from water-extractable PyOM (frac_13C_wx) by partitioning 
#the total CO2 emissions from the “Soil + 13C water-extractable PyOM” treatments (13C_wx_soil) into two sources: 
#CO2 emitted from the 13C labelled water-extractable PyOM source and CO2 emitted from the 
#unlabeled soil and non-water-extractable PyOM source 

#Equation for partitioning: sample_iCO2 = frac_13C_wx*pyom_13C_d13C + (1-frac_13C_wx)*mean_12C_bk_soil_iCO2.

#sample_iCO2 <- isotopic composititon of the total CO2 emited in 13_wx_soil treatments
#pyom_13C_d13C <- isotopic composition of 350 PyOM; proxy for isotopic composition of 
#CO2 emitted from the 13C labeled water-exctractable PyOM
#mean_12C_bk_soil_iCO2 <- isotopic composition of CO2 emitted from soil and unlabeled non-water-extractable PyOM

#We can calculate the mean_12C_bk_soil_iCO2 for each 
#cycle from “Soil + unlabeled PyOM” treatments (12C_bk_soil)
df.12C.bk = mineralization %>%
  filter(treatment == "12C_bk_soil") %>%
  dplyr::group_by(temperature,cycle) %>%
  dplyr::summarize(mean_12C_bk_soil_iCO2 = mean(sample_iCO2)) %>%
  dplyr::mutate(pyom_d13C = ifelse(temperature == 350, 
                                   pyom_350_13C_d13C,
                                   pyom_550_13C_d13C))

#Creating a new data frame for water-extractable PyOM-C partitioning 
df.wx = mineralization %>%
  dplyr::filter(treatment == "13C_wx_soil") %>%
  dplyr::filter(sample_num !=42) %>%
  dplyr::left_join(df.12C.bk) 

df.wx = df.wx %>%
  dplyr::mutate(frac_13C_wx = ((sample_iCO2 - mean_12C_bk_soil_iCO2) /(pyom_d13C - mean_12C_bk_soil_iCO2))) %>%
  dplyr::mutate(wx_C_mg = CO2C_mg_corr*frac_13C_wx) %>%
  #normalizing with the total amount of water-extractable PyOM-C in each jar
  dplyr::mutate(wx_C_mg_norm = ifelse(temperature == 350, wx_C_mg/total_wxC_g_350, wx_C_mg/total_wxC_g_550)) %>%
  dplyr::mutate(wx_C_mg_per_day = wx_C_mg_norm/time_int) %>%
  dplyr::select(c("temperature","cycle","sample_num", "replicate","wx_C_mg","wx_C_mg_norm", "wx_C_mg_per_day","time_days"))

ggplot(df.wx) + geom_point(aes(x = time_days, y = wx_C_mg_per_day)) + facet_wrap(~temperature)

#Cumulative 
df.wx.c <- df.wx %>%
  dplyr::group_by(sample_num)%>%
  dplyr::arrange(cycle)%>%
  dplyr::mutate(wx_C_mg_cml= cumsum(wx_C_mg)) %>%
  dplyr::mutate(wx_C_mg_norm_cml= cumsum(wx_C_mg_norm)) %>%
  dplyr::select(c("temperature","cycle","replicate","wx_C_mg_cml","wx_C_mg_norm_cml","time_days"))

ggplot(df.wx.c) + geom_point(aes(x = time_days, y = wx_C_mg_norm_cml))+facet_wrap(~temperature)

##################################################################################
# Partitioning non-water-extractable PyOM-C mineralized 
#################################################################################
#We can determine the fraction of CO2 emitted from non-water-extractable PyOM (frac_13C_nonwx) by partitioning 
#the total CO2 emissions from the “Soil + 13C non-water-extractable PyOM” treatments (13C_nox_soil) into two sources: 
#CO2 emitted from the 13C labelled non-water-extractable PyOM source and CO2 emitted from the 
#unlabeled soil and water-extractable PyOM source 

# Equation for partitioning: sample_iCO2 = frac_13C_nonwx*pyom_13C_d13C + (1-frac_13C_nonwx)*mean_12C_bk_soil_iCO2.

#sample_iCO2 <- isotopic composititon of the total CO2 emited in 13_nox_soil treatments
#pyom_13C_d13C <- isotopic composition of 13C labeleled PyOM; proxy for isotopic composition of 
#CO2 emitted from the 13C labeled non-water-exctractable PyOM
#mean_12C_bk_soil_iCO2 <- isotopic composition of CO2 emitted from soil and unlabeled water-extractable PyOM

#We use the values from “Soil + unlabeled PyOM” treatments (12C_bk_soil) for mean_12C_bk_soil_iCO2

#Creating a new data frame for non-water-extractable PyOM-C partitioning 
df.nonwx = mineralization %>%
  dplyr::filter(treatment == "13C_nox_soil") %>%
  dplyr::left_join(df.12C.bk)

df.nonwx = df.nonwx %>%
  dplyr::mutate(frac_13C_nonwx = ((sample_iCO2 - mean_12C_bk_soil_iCO2) /(pyom_d13C - mean_12C_bk_soil_iCO2))) %>%
  dplyr::mutate(nonwx_C_mg = CO2C_mg_corr*frac_13C_nonwx) %>%
  #normalizing with the total amount of non-water-extractable PyOM-C in each jar
  dplyr::mutate(nonwx_C_mg_norm = ifelse(temperature == 350, nonwx_C_mg/total_nonwxC_g_350, nonwx_C_mg/total_nonwxC_g_550)) %>%
  dplyr::mutate(nonwx_C_mg_per_day = nonwx_C_mg_norm/time_int ) %>%
  dplyr::select(c("temperature","cycle","sample_num", "replicate","nonwx_C_mg","nonwx_C_mg_norm","nonwx_C_mg_per_day","time_days"))


ggplot(df.nonwx) + geom_point(aes(x = time_days, y = nonwx_C_mg_per_day)) + facet_wrap(~temperature)

#Cumulative
df.nonwx.c <- df.nonwx %>%
  dplyr::group_by(sample_num)%>%
  dplyr::arrange(cycle)%>%
  dplyr::mutate(nonwx_C_mg_cml= cumsum(nonwx_C_mg))%>%
  dplyr::mutate(nonwx_C_mg_norm_cml= cumsum(nonwx_C_mg_norm))%>%
  dplyr::select(c("temperature","cycle","replicate","nonwx_C_mg_cml","nonwx_C_mg_norm_cml","time_days"))

ggplot(df.nonwx.c) + geom_point(aes(x = time_days, y = nonwx_C_mg_cml)) + facet_wrap(~temperature)

#################################################################################
# Partitioning soil C and bulk PyOM-C mineralized
#############################################################################
#We can determine the fraction of CO2 emitted from SOC (frac_soil) by partitioning the total CO2 emissions 
#from the “Soil + 13C PyOM” treatments (13C_bk_soil) into two sources: CO2 emitted from the 13C labelled PyOM 
#source and CO2 emitted from the unlabeled soil source

#Equation for partitioning:sample_iCO2 = frac_soil*soil_d13C + (1-frac_soil)*pyom_13C_d13C.
#sample_iCO2 <- isotopic composititon of the total CO2 emited in 13_bk_soil treatments
#pyom_13C_d13C <- isotopic composition of 13c labeled PyOM
#soil_d13C <- isotopic composition of soil

#Fraction of PyOM-C can be calculated as:
#frac_pyom = 1-frac_soil

#We will need an estimate of the soil_d13C. We can do this from the unamended "control soil" treatments (soil).
soil.control <- df.flux %>%
  dplyr::filter(treatment == "soil")

#ggplot(data = soil.control, aes(y=sample_iCO2,x=cycle,color=as.factor(sample_num)))+
#  geom_point() + facet_wrap(~temperature)

#Calculating soil_d13C for each cycle
soil.13C.cycle = soil.control %>%
  dplyr::filter(sample_num !=31) %>%
  dplyr::group_by(temperature,cycle) %>%
  dplyr::summarise(soil_d13C = mean(sample_iCO2), SD = sd(sample_iCO2))%>%
  dplyr::mutate(pyom_d13C = ifelse(temperature == 350,pyom_350_13C_d13C,pyom_550_13C_d13C ))

#Creating a new data frame for soil C and bulk PyOM-C partitioning 
df.soil = mineralization %>%
  dplyr::filter(treatment == "13C_bk_soil") %>%
  dplyr::left_join(soil.13C.cycle) %>%
  dplyr::mutate(frac_soil = ((sample_iCO2 - pyom_d13C) /(soil_d13C - pyom_d13C))) %>%
  dplyr::mutate(frac_pyom = 1-frac_soil) %>% #fraction coming from PyOM-C
  dplyr::mutate(soilC_mg = CO2C_mg_corr*frac_soil) %>%
  dplyr::mutate(pyom_C_mg = CO2C_mg_corr*frac_pyom) %>%
  dplyr::mutate(soilC_mg_norm = soilC_mg/total_soilC_g) %>% #normalized by soil C in each jar 
  dplyr::mutate(pyom_C_mg_norm = ifelse(temperature == 350, pyom_C_mg/total_pyomC_g_350, pyom_C_mg/total_pyomC_g_550)) %>% 
  dplyr::mutate(soilC_mg_per_day = soilC_mg_norm/time_int) %>%
  dplyr::mutate(pyom_C_mg_per_day = pyom_C_mg_norm/time_int) %>%
  dplyr::select(c("temperature","cycle","sample_num", "replicate",
                  "soilC_mg","pyom_C_mg", "soilC_mg_norm","pyom_C_mg_norm","soilC_mg_per_day","pyom_C_mg_per_day","time_days"))

ggplot(df.soil) + geom_point(aes(x = time_days, y = soilC_mg_per_day))+ facet_wrap(~temperature)
ggplot(df.soil) + geom_point(aes(x = time_days, y = pyom_C_mg_per_day))+ facet_wrap(~temperature)

#Cumulative
df.soil.c <- df.soil %>%
  dplyr::group_by(sample_num)%>%
  dplyr::arrange(cycle)%>%
  dplyr::mutate(soil_C_mg_cml= cumsum(soilC_mg))%>%
  dplyr::mutate(soil_C_mg_norm_cml= cumsum(soilC_mg_norm))%>%
  dplyr::select(c("temperature","cycle","replicate","soil_C_mg_cml","soil_C_mg_norm_cml","time_days"))

ggplot(df.soil.c) + geom_point(aes(x = time_days, y = soil_C_mg_cml)) + facet_wrap(~temperature)

# Calculating cumulative PyOM-C mineralized
df.pyom.c <- df.soil %>%
  dplyr::group_by(sample_num)%>%
  dplyr::arrange(cycle)%>%
  dplyr::mutate(pyom_C_mg_cml= cumsum(pyom_C_mg)) %>%
  dplyr::mutate(pyom_C_mg_norm_cml= cumsum(pyom_C_mg_norm)) %>%
  dplyr::select(c("temperature","cycle","replicate","pyom_C_mg_cml","pyom_C_mg_norm_cml","time_days"))

ggplot(df.pyom.c) + geom_point(aes(x = time_days, y = pyom_C_mg_cml)) + facet_wrap(~temperature)

#################################################################################
# Calculating rate and cumulative soil C mineralized from unamended soil treatments
#############################################################################
# We can do this from the control soil jars.
df.soil.control <- mineralization %>%
  dplyr::filter(treatment == "soil") %>%
  dplyr::mutate(CO2C_mg_norm = CO2C_mg_corr/total_soilC_g) %>%
  dplyr::mutate(CO2C_mg_per_day = CO2C_mg_norm/time_int) %>%
  dplyr::select(c("temperature","cycle","sample_num", "replicate","CO2C_mg_corr","CO2C_mg_norm","CO2C_mg_per_day","time_days"))

ggplot(df.soil.control) + geom_point(aes(x = time_days, y = CO2C_mg_per_day)) + facet_wrap(~temperature)

#cumulative
df.soil.control.c <- df.soil.control %>%
  dplyr::group_by(sample_num)%>%
  dplyr::arrange(cycle)%>%
  dplyr::mutate(CO2C_mg_cml= cumsum(CO2C_mg_corr)) %>%
  dplyr::mutate(CO2C_mg_norm_cml= cumsum(CO2C_mg_norm)) %>%
  dplyr::select(c("temperature","cycle","replicate","CO2C_mg_cml","CO2C_mg_norm_cml","time_days"))

ggplot(df.soil.control.c) + geom_point(aes(x = time_days, y = CO2C_mg_cml)) + facet_wrap(~temperature)

##################################################################################################
#####################################################################################
###Merging data frames
df.pyom <- df.soil[,c("temperature", "cycle", "sample_num", "replicate",
                      "pyom_C_mg","pyom_C_mg_norm","pyom_C_mg_per_day", "time_days")]
df.soil <- df.soil %>%
  select(-c("pyom_C_mg", "pyom_C_mg_norm", "pyom_C_mg_per_day"))
#Merging data frames for mineralized C rates
df.merged.cols = c("temperature", "cycle", "sample_num", "replicate", "CO2C_mg","CO2C_mg_norm", 
                   "CO2C_mg_per_day","time_days")

colnames(df.wx) <- df.merged.cols
colnames(df.nonwx) <- df.merged.cols
colnames(df.soil) <- df.merged.cols
colnames(df.soil.control) <- df.merged.cols
colnames(df.pyom) <- df.merged.cols

df.wx['treatment'] = 'water-extractable PyOM-C'
df.nonwx['treatment'] = 'non-water-extractable PyOM-C'
df.soil['treatment'] = 'SOC'
df.soil.control['treatment'] = 'SOC(no PyOM)'
df.pyom['treatment'] = 'PyOM-C'

df.final = rbind(df.wx, df.nonwx, df.soil, df.soil.control, df.pyom)
#write.csv(df.final, "rate_respired.csv")


#Merging data frames for cumulative C figure
df.merged.cols = c("sample_num", "temperature", "cycle", "replicate", "C_mg_cml","C_mg_norm_cml", "time_days")

colnames(df.wx.c) <- df.merged.cols
colnames(df.nonwx.c) <- df.merged.cols
colnames(df.soil.c) <- df.merged.cols
colnames(df.soil.control.c) <- df.merged.cols
colnames(df.pyom.c) <- df.merged.cols

df.wx.c['treatment'] = 'water-extractable PyOM-C'
df.nonwx.c['treatment'] = 'non-water-extractable PyOM-C'
df.soil.c['treatment'] = 'SOC'
df.soil.control.c['treatment'] = 'SOC(no PyOM)'
df.pyom.c['treatment'] = 'PyOM-C'

df.final.c = rbind(df.wx.c, df.nonwx.c, df.soil.c, df.soil.control.c, df.pyom.c)
#write.csv(df.final.c, "cml_respired.csv")

df.final.c.last = df.final.c %>%
  group_by(temperature, treatment) %>%
  arrange(cycle) %>%
  filter(cycle == last(cycle))
#write.csv(df.final.c.last, "cml_respired_lastcyc.csv")



