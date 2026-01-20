#load your libraries
library(tidyverse)
library(cowplot)
library(dplyr) 
#set your working directory 
setwd("C:/Users/et491/OneDrive/Documents/MbyRes/tecan assays/tecan assay 1 72 hr")

#setting your dataframes
df = read.csv("reformatted_cross_resistance_30.1.25.csv")
dictionary <- read.csv("tecan_1_dictionary.csv")
#join the dataframes
df <- df %>% left_join(dictionary)
#mutating the wells to being factors and more readable for the computer (no I don't really understand how this code works)
df$well= factor(df$well, levels=paste(rep(LETTERS[1:8], each = length(seq(1, 12))), seq(1, 12), sep = ""))
df <- df %>% mutate(row=gsub('^([A-H])(.*)', '\\1', well), column=factor(gsub('^([A-H])(.*)', '\\2', well), levels=1:12))
#creating a blanks dataframe to normalise your data to the blank 
blks_df <- df %>% filter(treatment == "blk") %>% group_by(time_hours) %>% summarise(mean_blk_od = mean(mean_od))
#normalising your data to the blank 
df <- df %>% left_join(blks_df) %>% mutate(norm_od = mean_od - mean_blk_od)
#woohoo your dataframe should be ready to go to use for these graphs now :) 
#Creating a growth curve per well
overall_growth_curves <- ggplot(df, aes(x = time_hours, y = norm_od)) +
  geom_line(color='#0072B2') +
  scale_x_continuous('Time (Hours)') +
  scale_y_continuous('OD600') +
  ylab("Mean OD600") +
  facet_wrap(~well, ncol = 12)
overall_growth_curves
#Mean Growth Curve with StDev
#You need to create a new dataframe here to group things and mutate it to be the mean 
df_mean <- df %>% group_by(treatment, time_hours) %>% summarise(time_hours = time_hours, treatment = treatment, OD600 = mean(norm_od), STDEV = sd(norm_od))
#removing the blk wells
df_mean <- df_mean %>% filter(!treatment == "blk")

## calculating the standard error and CI95
df_SE <- df_mean %>% mutate(SE = STDEV/sqrt(6))
df_CI <- df_SE %>% mutate(CI95 = 1.96*SE)
df_CI <- df_CI %>% mutate(CI95_min = OD600 - CI95) %>% mutate(CI95_max = OD600 + CI95)


####separating out all data into two datasets that are just the hsdR and just the pilA galU 

df_hsdR_curve<- df_CI %>%  filter (!treatment == "blk") %>% 
  filter(!(treatment %in% c("pilA galU", "pilA galU plus 272", "pilA galU plus 1276", "pilA galU plus 272 and 1276", "1276 control", "272 control"))) %>% 
  mutate(treatment = factor(treatment, 
                            levels = c("hsdR", "hsdR plus 272", "hsdR plus 1276", "hsdR plus 272 and 1276")))

df_pilAgalU_curve<- df_CI %>%  filter (!treatment == "blk") %>% 
  filter(!(treatment %in% c("hsdR", "hsdR plus 272", "hsdR plus 1276", "hsdR plus 272 and 1276", "1276 control", "272 control"))) %>% 
  mutate(treatment = factor(treatment, 
                            levels = c("pilA galU", "pilA galU plus 272", "pilA galU plus 1276", "pilA galU plus 272 and 1276")))



#four separate plots all the same colour, no legend, 95% confidence interval ribbons

custom_labels <- as_labeller(c(
  "hsdR" = "'Parent strain'",
  "hsdR plus 272" = "'Parent strain + Tor'",
  "hsdR plus 1276" = "'Parent strain + Vale'",
  "hsdR plus 272 and 1276" = "'Parent strain + Tor + Vale'",
  "pilA galU" = "Delta * italic('pilA') ~ Delta * italic('galU')",
  "pilA galU plus 272" = "Delta * italic('pilA') ~ Delta * italic('galU') ~ '+ Tor'",
  "pilA galU plus 1276" = "Delta * italic('pilA') ~ Delta * italic('galU') ~ '+ Vale'",
  "pilA galU plus 272 and 1276" = "Delta * italic('pilA') ~ Delta * italic('galU') ~ '+ Tor + Vale'"
), label_parsed)


mean_growth_curve_hsdR <- ggplot(data = df_hsdR_curve, aes(x = time_hours, y = OD600)) +
  geom_line(size=1.5, colour = "#56B4E9") +
  geom_ribbon(aes(ymin= CI95_min, ymax = CI95_max), alpha = 0.2, show.legend = FALSE, fill = "#56B4E9", colour = "#56B4E9")+
  scale_y_continuous('Mean OD600 (n=10)') +
  scale_x_continuous('Time (Hours)', breaks = seq(0,72, by=24), expand = expansion(mult = c(0.02, 0.02))) +
  theme(plot.margin = margin(10, 30, 10, 10)) +
  facet_wrap(~treatment, ncol = 1, labeller = custom_labels) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour="#aca9a9", size=0.4),
        panel.grid.minor = element_line(colour="#aca9a9", size=0.05),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        strip.text = element_text(size=20))
mean_growth_curve_hsdR


mean_growth_curve_pilAgalU <- ggplot(data = df_pilAgalU_curve, aes(x = time_hours, y = OD600)) +
  geom_line(size=1.5, colour = "#56B4E9") +
  geom_ribbon(aes(ymin= CI95_min, ymax = CI95_max), alpha = 0.2, show.legend = FALSE, fill = "#56B4E9", colour = "#56B4E9")+
  scale_y_continuous('Mean OD600 (n=10)') +
  scale_x_continuous('Time (Hours)', breaks = seq(0,72, by=24), expand = expansion(mult = c(0.02, 0.02))) +
  theme(plot.margin = margin(10, 30, 10, 10)) +
  facet_wrap(~treatment, ncol = 1, labeller = custom_labels) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.grid.major = element_line(colour="#aca9a9", size=0.4),
        panel.grid.minor = element_line(colour="#aca9a9", size=0.05),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        strip.text = element_text(size=20))
mean_growth_curve_pilAgalU


library(patchwork)
figure_6 <- mean_growth_curve_hsdR | mean_growth_curve_pilAgalU
figure_6
ggsave("C:/Users/et491/OneDrive/Documents/MbyRes/Papers/Figures LPS paper/figure 6 growth curves.pdf", width = 15, height = 10, device = cairo_pdf)
