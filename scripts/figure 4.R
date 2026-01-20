
library(dplyr)
library(ggplot2)
library(cowplot)
setwd("C:/Users/et491/OneDrive/Documents/MbyRes/LPS ext and sequencing experiment/resistant mutant spot assays/1276 and 272 PFUs on resistant mutants/4 repeats")
df = read.csv("1276 and 272 PFUs with repeats.csv")

df_mean <- df %>% group_by(phage, strain) %>% mutate(mean_EoP = mean(EoP))

df_ordered <- df_mean %>%
  mutate(strain = factor(strain, levels = c("hsdR", "pilAgalU", "R272.1","R272.3","R272.5","R1276.1",  "R1276.2","R1276.4")))


df_adjusted <- df_ordered %>%
  mutate(eop_for_plot = mean_EoP + 1e-4)


P1<- ggplot(df_adjusted, aes(y = as.factor(phage), x = strain, fill = eop_for_plot)) +
  geom_tile(colour = "white", lwd = 1, linetype = 1) +
  geom_text(aes(label = round(mean_EoP, 2)), color = "white", size = 8) + 
  labs(y = "Phage", x = "Host strain", fill = "EoP", title = "Phage resistant mutant susceptibility to Tor and Vale") +
  scale_fill_gradientn(
    colours = c("#C6DBEF","#6BAED6","#2171B5", "#08306B"),
    values = scales::rescale(c(1e-4, 0.01, 0.1, 0.5, 1)),                
    trans = "log",
    breaks = c(1e-4, 0.01, 0.1, 1),
    labels = c("0", "0.01", "0.1", "1")) +
  scale_y_discrete(labels = c("Tor", "Vale"), expand = c(0,0)) +
  scale_x_discrete(labels=c(
    "Parent strain",
    expression(Delta * italic(pilA) * " " * Delta * italic(galU)),
    "Tor_R1",
    "Tor_R2",
    "Tor_R3",
    "Vale_R1",
    "Vale_R2",
    "Vale_R3"), expand = c(0,0)) +
  theme_cowplot(25) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"),
        axis.text.x = element_text(size=25, angle=30, vjust=1,hjust=1))



setwd("C:/Users/et491/OneDrive/Documents/MbyRes/LPS ext and sequencing experiment/double phage resistant mutants/attempt 3")
df2 = read.csv("double resistant mutant PFUs.csv")


df2_mean <- df2 %>% group_by(phage, strain) %>% mutate(mean_EoP = mean(EoP))

df2_ordered <- df2_mean %>%
  mutate(strain = factor(strain, levels = c("hsdR", "pilAgalU", "R272.1","R272.3","R272.5","R1276.1",  "R1276.2","R1276.4")))

df2_adjusted <- df2_ordered %>%
  mutate(eop_for_plot = mean_EoP + 1e-4)

P2 <- ggplot(df2_adjusted, aes(y = as.factor(phage), x = strain, fill = eop_for_plot)) +
  geom_tile(colour = "white", lwd = 1, linetype = 1) +
  geom_text(aes(label = round(mean_EoP, 2)), color = "white", size = 8) + 
  labs(y = "Phage", x = "Host strain", fill = "EoP", title = "Phage resistant mutant susceptibility to Tor and Vale after alternate \nphage exposure") +
  scale_fill_gradientn(
    colours = c("#C6DBEF","#6BAED6","#2171B5", "#08306B"),
    values = scales::rescale(c(1e-4, 0.01, 0.1, 0.5, 1)),                
    trans = "log",
    breaks = c(1e-4, 0.01, 0.1, 1),
    labels = c("0", "0.01", "0.1", "1")) +
  scale_y_discrete(labels = c("Tor", "Vale"), expand = c(0,0)) +
  scale_x_discrete(labels=c(
    "Parent strain",
    expression(Delta * italic(pilA) * " " * Delta * italic(galU)),
    expression(atop("Tor_R1","Vale-exposed")),
    expression(atop("Tor_R2","Vale-exposed")),
    expression(atop("Tor_R3","Vale-exposed")),
    expression(atop("Vale_R1", "Tor-exposed")),
    expression(atop("Vale_R2", "Tor-exposed")),
    expression(atop("Vale_R3", "Tor-exposed"))), expand = c(0,0)) +
  theme_cowplot(25) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"),
        axis.text.x = element_text(size=23, angle=30, vjust=1,hjust=1))


library(patchwork)
figure_4 <- P1 / P2  
figure_4 + plot_annotation(tag_levels = list(c("B", "C")))
ggsave("C:/Users/et491/OneDrive/Documents/MbyRes/Papers/Figures LPS paper/figure 4 heatmaps v2.pdf", width = 17, height = 13, device = cairo_pdf)
