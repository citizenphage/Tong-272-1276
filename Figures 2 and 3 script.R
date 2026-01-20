library(ggtext)
library(tidyverse)
library(DescTools)
library(cowplot)


#all strains - oprm, pila, wbpl

setwd("C:/Users/et491/OneDrive/Documents/MbyRes/receptor identification spot assays/96-well plaque assays")
df = read.csv("3r 96 well plaque assay all strains.csv")

#removing any row with blanks
df <- df %>% filter(Phage != "blk")

#remove e.coli and T7
df <- df %>% filter(Strain != "ECBW") %>%
  filter(Phage != "T7")

#now create a df of just the hsdr (assuming WT) strain so we can compare to it.
df_WT <- df %>% 
  filter(Strain == "PAO1_hsdR") %>%
  rename(WT_value = plaque_type) %>% 
  select("Phage", "WT_value")

#and left join that in
df <- df %>% left_join(df_WT, by = "Phage")

#dont need some of the data in df so doing a bit of tidying
df <- df %>% select(-c("well", "isolation_strain"))

#now creating one value for evidence of phage killing
df <- df %>%
  mutate(Evidence_test = ifelse(plaque_type == "N", 0, 1)) %>%
  mutate(WT_test = ifelse(WT_value == "N", 0, 1)) %>%
  select(-c("plaque_type", "WT_value"))

#finally run a comparison to see if they are different#
df <- df %>%
  mutate(Difference = ifelse(Evidence_test == WT_test, "No Difference to Parent Strain", "Infection Inhibited"))

#factor time to alter order of values in plot, also remove 933 as no infection in control, finally remove WT
df <- df %>%
  filter(!Phage %in% c("950","933", "953", "952", "947", "946", "934", "162")) %>%
  mutate(Difference = factor(Difference, levels = c("No Difference to Parent Strain", "Infection Inhibited"))) %>%
  mutate(Strain = factor(Strain, levels = c("PAO1_hsdR_oprM", "PAO1_hsdR_pilA", "PAO1_hsdR_wbpL", "PAO1_hsdR_galU", "PAO1_hsdR_algC", "PAO1_hsdR_pilA_galU", "PAO1_hsdR_pilA_algC"))) %>%
  mutate(Phage = factor(Phage, levels = c("215","252","253","254", "161", "211","212","214","216L","217L","56++","71","165","282","70","72","75","113","173","176","213","272","273","274","275","276","277","278","279","280","281","283","285","286","287","935","936","937","938","939","940","941","942","943","944","945","949","951","954","973"))) %>%
  filter(Strain != "PAO1_hsdr")



y_labels <- setNames(df$fullnames, df$Phage)

ggplot(df, aes(x = Strain, y = Phage, fill = Difference)) +
  geom_tile(colour = "white", lwd = 1, linetype = 1) +
  scale_fill_manual(values = c("lightgrey", "#fb8072")) +
  labs(x = "Bacterial Strain", y = "Phage", fill = "Infection Status") +
  scale_x_discrete(labels= c(
   expression(Delta * italic(oprM)),
   expression(Delta * italic(pilA)),
   expression(Delta * italic(wbpL)),
    expression(Delta * italic(galU)),
    expression(Delta * italic(algC)),
   expression(Delta * italic(pilA) ~ Delta * italic(galU)),
    expression(Delta * italic(pilA) ~ Delta * italic(algC))), expand = c(0,0)) +
  scale_y_discrete(labels=y_labels) +
  theme_cowplot(15) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

ggsave("C:/Users/et491/OneDrive/Documents/MbyRes/Papers/Figures LPS paper/receptor identification heatmap.pdf", width = 13, height = 9.5, device = cairo_pdf)


#for phages isolated on pilA galU
setwd("C:/Users/et491/OneDrive/Documents/MbyRes/Phage enrichment/phage hunt with antibiotic on pila galu/phage on all strains/repeat on 23rd jan 25")
df2 = read.csv("antibiotic phages on all strains 2.csv")


#factor time to alter order of values in plot, also remove 945 as no infection in control, finally remove WT
df2 <- df2 %>%
  mutate(Strain = factor(Strain, levels = c("PAO1_hsdR", "PAO1_hsdR_oprM", "PAO1_hsdR_pilA", "PAO1_hsdR_wbpL", "PAO1_hsdR_galU", "PAO1_hsdR_algC", "PAO1_hsdR_pilA_galU", "PAO1_hsdR_pilA_algC"))) %>%
  filter(Phage != "CPL01283") %>%
  mutate(Phage = factor(Phage, levels = c("CPL01279L", "CPL01280", "CPL01282","(Vale) CPL01276","CPL01279S", "CPL01275", "CPL01278", "CPL01277", "CPL01284","CPL01271"))) %>%
  mutate(Evidence_test = ifelse(plaque_type == "N", "P", "W"))


ggplot(df2, aes(x = Strain, y = Phage, fill = Evidence_test)) +
  geom_tile(colour = "white", lwd = 1, linetype = 1) +
  scale_fill_manual(values = c("lightgrey", "#009e73"), labels = c("no plaques", "plaques")) +
  labs(x = "Bacterial Strain", y = "Phage", fill = "Plaque formation") +
  scale_x_discrete(labels= c(
    "Parent Strain",
    expression(Delta * italic(oprM)),
    expression(Delta * italic(pilA)),
    expression(Delta * italic(wbpL)),
    expression(Delta * italic(galU)),
    expression(Delta * italic(algC)),
    expression(Delta * italic(pilA) ~ Delta * italic(galU)),
    expression(Delta * italic(pilA) ~ Delta * italic(algC))
  ), expand = c(0,0)) +
  scale_y_discrete(labels=y_labels) +
  theme_cowplot(15) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"))

ggsave("C:/Users/et491/OneDrive/Documents/MbyRes/Papers/Figures LPS paper/new phage heatmap.pdf", width = 13, height = 8, device = cairo_pdf)
