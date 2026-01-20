library(ggplot2)
library(dplyr)


setwd("C:/Users/et491/OneDrive/Documents/MbyRes/Phage oxford 2025")
df.mic = read.csv("antibiotic results for bar graph.csv")
  
ggplot(df.mic, aes(x = strain, y = MIC)) + 
  geom_bar(stat = "identity", fill="#009999") +
  labs(x = "Bacterial strain", y = "MIC of gentamycin (µg/mL)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels=c("Parent strain", "Tor_R1", "Tor_R2", "Tor_R3")) +
  theme_bw(20) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5))

ggsave("C:/Users/et491/OneDrive/Documents/MbyRes/Papers/Figures LPS paper/antibiotic bar chart.pdf", width = 8.5, height = 8.5, device = cairo_pdf)
