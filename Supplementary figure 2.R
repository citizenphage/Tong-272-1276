library(tidyverse)
library(readxl)
setwd("C:/Users/et491/OneDrive/Documents/MbyRes/tecan assays/Knockout panel growth curves")

df <- read_excel("tecan sunrise Robin 3.4.25.xlsx", sheet = "RAW")
dictionary <- read.csv("tecan sunrise dictionary 2.csv")

#changing to long format
df_long <- pivot_longer(df, !Time, names_to = "Well", values_to = "Absorbance")

#join the dataframes
df_dictionary <- df_long %>% left_join(dictionary)


#mutating the wells to being factors and more readable for the computer (no I don't really understand how this code works)
df_dictionary$Well= factor(df_dictionary$Well, levels=paste(rep(LETTERS[1:8], each = length(seq(1, 12))), seq(1, 12), sep = ""))
df <- df_dictionary %>% mutate(row=gsub('^([A-H])(.*)', '\\1', Well), column=factor(gsub('^([A-H])(.*)', '\\2', Well), levels=1:12))


#creating a blanks dataframe to normalise your data to the blank 
blks_df <- df %>% filter(Phage == "BLK") %>% group_by(Time) %>% summarise(mean_blk_od = mean(Absorbance))

#normalising your data to the blank 
df <- df %>% left_join(blks_df) %>% mutate(norm_od = Absorbance - mean_blk_od)

df$Time <- as.numeric(as.character(df$Time))

#Creating a growth curve per well
overall_growth_curves <- ggplot(df, aes(x = Time, y = norm_od)) +
  geom_line(color='#0072B2') +
  scale_x_continuous('Time (Hours)') +
  scale_y_continuous('OD600') +
  ylab("Mean OD600") +
  facet_wrap(~Well, ncol = 12)
overall_growth_curves

#Mean Growth Curve with StDev
#You need to create a new dataframe here to group things and mutate it to be the mean 
df_mean <- df %>% group_by(Phage, Time) %>% summarise(Time = Time, Phage = Phage, OD600 = mean(norm_od), STDEV = sd(norm_od))
#removing the blk wells
df_mean <- df_mean %>% filter(!Phage %in% c("BLK", "Protector"))

## calculating the standard error and CI95
df_SE <- df_mean %>% mutate(SE = STDEV/sqrt(6))
df_CI <- df_SE %>% mutate(CI95 = 1.96*SE)
df_CI <- df_CI %>% mutate(CI95_min = OD600 - CI95) %>% mutate(CI95_max = OD600 + CI95)

#changing seconds to hours
df_hours <- df_CI %>% mutate(time_hours = Time/3600)

#reordering mutants
df_hours  <- df_hours %>% mutate(Phage=factor(Phage, levels = c("PAO1", "PAO1 hsdR", "PAO1 oprM", "PAO1 pilA", "PAO1 wbpL", "PAO1 galU", "PAO1 pilA galU", "PAO1 algC", "PAO1 pilA algC")))

#custom labels for mutants

custom_labels <- as_labeller(c(
  "PAO1" = "PAO1",
  "PAO1 hsdR" = "Delta * italic('hsdR')",
  "PAO1 oprM" = "Delta * italic('oprM')",
  "PAO1 pilA" = "Delta * italic('pilA')",
  "PAO1 wbpL" = "Delta * italic('wbpL')",
  "PAO1 galU" = "Delta * italic('galU')",
  "PAO1 pilA galU" = "Delta * italic('pilA') ~ Delta * italic('galU')",
  "PAO1 algC" = "Delta * italic('algC')",
  "PAO1 pilA algC" = "Delta * italic('pilA') ~ Delta * italic('algC')"
), label_parsed)

#cutting down hours
df_16 <- df_hours %>% 
  filter(time_hours <= 16)

all_facetted_growth_curves_16 <- ggplot(data = df_16, aes(x = time_hours, y = OD600)) +
  geom_line(size=1.5, colour = "#56B4E9") +
  geom_ribbon(aes(ymin= CI95_min, ymax = CI95_max), alpha = 0.2, show.legend = FALSE, fill = "#56B4E9", colour = "#56B4E9")+
  scale_y_continuous('Mean OD600 (n=6)') +
  scale_x_continuous('Time (hours)',  breaks = seq(0,16, by=4), expand = expansion(mult = c(0.02, 0.02))) +
  theme(plot.margin = margin(10, 30, 10, 10)) +
  facet_wrap(~Phage, ncol = 3, labeller=custom_labels) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour="#aca9a9", size=0.4),
        panel.grid.minor = element_line(colour="#aca9a9", size=0.05),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        panel.spacing.x = unit(1.2, "lines"))
all_facetted_growth_curves_16

ggsave("C:/Users/et491/OneDrive/Documents/MbyRes/Papers/Figures LPS paper/knockout growth curves.pdf", width = 15.5, height = 11, device = cairo_pdf)



#Using growthcurver to determine whether there is a significant difference in growth curves between deletion mutants and wildtype PAO1

library(growthcurver)
setwd("C:/Users/et491/OneDrive/Documents/MbyRes/tecan assays/Knockout panel growth curves")


#initial input is in a wide format with one column for time and then one column for each well, with OD readings at each time point. A dictionary stating what treatment ("Phage") is in each well
df <- read_excel("tecan sunrise Robin 3.4.25.xlsx", sheet = "RAW")
dictionary <- read.csv("tecan sunrise dictionary 2.csv")

#changing seconds to hours and converting time to a numeric variable
df$Time <- as.numeric(as.character(df$Time))
df <- df %>% mutate(Time = Time/3600)

#changing to long format
df_long <- pivot_longer(df, !Time, names_to = "Well", values_to = "Absorbance")

#join the dictionary to the dataframe
df_dictionary <- df_long %>% left_join(dictionary)

#mutating the wells to being factors and more readable for the computer (no I don't really understand how this code works)
df_dictionary$Well= factor(df_dictionary$Well, levels=paste(rep(LETTERS[1:8], each = length(seq(1, 12))), seq(1, 12), sep = ""))
df <- df_dictionary %>% mutate(row=gsub('^([A-H])(.*)', '\\1', Well), column=factor(gsub('^([A-H])(.*)', '\\2', Well), levels=1:12))

#creating a blanks dataframe to normalise your data to the blank 
blks_df <- df %>% filter(Phage == "BLK") %>% group_by(Time) %>% summarise(mean_blk_od = mean(Absorbance))

#normalising your data to the blank 
df_norm <- df %>% left_join(blks_df) %>% mutate(norm_od = Absorbance - mean_blk_od)
df_norm <- df_norm %>% filter(!Phage %in% c("Protector", "BLK") )

#simplify so the only columns are OD, time, Well
df_simplified <- df_norm %>% 
  select(norm_od, Time, Well)

#converting back to a wide format for growthcurver with columns for each Well and for time
df_wide <- df_simplified %>%
  pivot_wider(names_from = Well, values_from = norm_od)

#creating growthcurver output grouped by treatment, with carrying capacity, intitial size, intrinsic growth rate, area under curve, and doubling time
gc_out_phage <- SummarizeGrowthByPlate(df_wide) %>%
  as.data.frame() %>%
  rename(Well = sample) %>%
  left_join(dictionary) %>%
  group_by(Phage) %>%
  reframe(mean_carrying_capacity = mean(k), sd_carrying_capacity=sd(k), se_carrying_capacity=sd_carrying_capacity/sqrt(6),CI95_carrying_capacity = 1.96*se_carrying_capacity,
          mean_initial_size = mean(n0), sd_initial_size = mean(n0),se_initial_size=sd_initial_size/sqrt(6),CI95_initial_size = 1.96*se_initial_size,
          mean_growth_rate = mean(r), sd_growth_rate = sd(r),se_growth_rate=sd_growth_rate/sqrt(6),CI95_growth_rate = 1.96*se_growth_rate,
          mean_auc = mean(auc_e), sd_auc = sd(auc_e),se_auc=sd_auc/sqrt(6),CI95_auc = 1.96*se_auc,
          mean_DT = mean(log(2)/r), sd_DT = sd(log(2)/r), se_DT=sd_DT/sqrt(6),CI95_DT = 1.96*se_DT,) %>%
  mutate(Phage = factor(Phage, levels = c("PAO1", "PAO1 hsdR", "PAO1 oprM", "PAO1 pilA", "PAO1 wbpL", "PAO1 galU", "PAO1 pilA galU", "PAO1 algC", "PAO1 pilA algC")) )


gc_out_simple <- gc_out_phage %>% select(Phage, mean_carrying_capacity,CI95_carrying_capacity, mean_initial_size, CI95_initial_size, mean_growth_rate, CI95_growth_rate, mean_auc, CI95_auc, mean_DT, CI95_DT)

#creating an output with each well separately for  linear models
gc_out_well <- SummarizeGrowthByPlate(df_wide) %>%
  as.data.frame() %>%
  rename(Well = sample) %>%
  left_join(dictionary) %>%
  reframe(carrying_capacity = k,
          initial_size = n0,
          growth_rate = r,
          auc = auc_e,
          DT = log(2)/r,
          strain = Phage)


#linear models to see if growth parameters differ significantly between strains
lmDT = lm(DT~strain, data = gc_out_well) #Create the linear regression
summary(lmDT) #Review the results

lmk = lm(carrying_capacity~strain, data = gc_out_well) #Create the linear regression
summary(lmk) #Review the results

lmgr = lm(growth_rate~strain, data = gc_out_well) #Create the linear regression
summary(lmgr) #Review the results

lmauc = lm(auc~strain, data = gc_out_well) #Create the linear regression
summary(lmauc) #Review the results