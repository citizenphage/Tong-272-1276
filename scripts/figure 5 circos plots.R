#R script to generate circos plots (comparative genomics)
#Improved version with better visualization and structure

#svg(width = 10, height = 10, file = "/Users/lb808/Documents/TempertonLab/EllieTong/circos_plot/PAO1_hsdr/R272_snp_del_V2.svg")

#Load libraries
library(circlize)
library(data.table)

#Set global circos parameters
circos.par(gap.degree = 2, cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 0.01), 
           start.degree = 90, canvas.xlim = c(-1.1, 1.1), canvas.ylim = c(-1.1, 1.1))

bg_color <- "#f7f7f7"     #Background color
del_color <- "#d9d9d9"    #Light gray for deletions
snp_color <- "#000000"    #Black for SNPs
track_border <- "#7f7f7f"

#Load reference genes from BED file, I was thinking on plot the proteins, but at this scale is not visually appropiate. 
hsdR_genes <- read.table("/Users/et491/OneDrive/Documents/MbyRes/Papers/LPS paper/EllieT_circos/hsdR_genes.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(hsdR_genes) <- c("chr", "start", "end")

#Load SNP/DEL data
load_snp_data <- function(file) {
  dt <- fread(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(dt) <- c("Chrom", "Start", "End", "Type")
  return(dt)
}

R272_1 <- load_snp_data("/Users/et491/OneDrive/Documents/MbyRes/Papers/LPS paper/EllieT_circos/PAO1_R272_1_rc_cut.snp.txt")
R272_3 <- load_snp_data("/Users/et491/OneDrive/Documents/MbyRes/Papers/LPS paper/EllieT_circos/PAO1_R272_3_cut.snp.txt")
R272_5 <- load_snp_data("/Users/et491/OneDrive/Documents/MbyRes/Papers/LPS paper/EllieT_circos/PAO1_R272_5_cut.snp.txt")

#Combine all data to get proper xlim
all_data <- rbind(R272_1, R272_3, R272_5)
chromosomes <- unique(all_data$Chrom)
max_pos <- tapply(all_data$End, all_data$Chrom, max)

#First, let's find the maximum position across all chromosomes
max_length <- max(hsdR_genes$end)

#Circos parameters with top orientation (0 on top)
circos.par("start.degree" = 90,  #This puts zero at the top
           "track.height" = 0.15,
           "gap.degree" = 5,
           "cell.padding" = c(0, 0, 0, 0))

#initialize circos with this maximum length for all chromosomes
circos.initialize(factors = chromosomes, 
                 xlim = cbind(rep(0, length(chromosomes)), 
                            rep(max_length, length(chromosomes))))
#Add axis
circos.track(ylim = c(0, 1), bg.col = bg_color, bg.border = NA, 
             panel.fun = function(x, y) {
               chr = CELL_META$sector.index
               xlim = CELL_META$xlim
               #Create breaks every 500kb
               axis_breaks <- seq(0, max_length, by = 5e5)
               #Ensure the last break is exactly at max_length
               if(max(axis_breaks) < max_length) {
                 axis_breaks <- c(axis_breaks, max_length)
               }
               #Create labels (show every 1Mb + the max length)
               axis_labels <- ifelse(axis_breaks %% 1e6 == 0 | axis_breaks == max_length,
                                    paste0(round(axis_breaks/1e6, 1), "Mb"), 
                                    "")
               circos.axis(h = "top", 
                          major.at = axis_breaks,
                          labels = axis_labels, 
                          labels.cex = 1.2,  #Increase here for the axis fonts
                          major.tick.length = 0.2)
               circos.text(mean(xlim), 1.2, chr, cex = 0.8, niceFacing = TRUE)
             }, track.height = 0.08)

#combine deletions and SNPs. Add frame (or background color)
plot_combined_track <- function(data, track_num) {
  circos.track(factors = data$Chrom, ylim = c(0, 1), 
               bg.border = track_border, bg.col = bg_color,
               panel.fun = function(x, y) {
                 current_chr <- CELL_META$sector.index
                 chr_data <- data[data$Chrom == current_chr,]
                 #Plot deletions (gray blocks)
                 dels <- chr_data[chr_data$Type == "DEL",]
                 if(nrow(dels) > 0) {
                   circos.rect(dels$Start, 0, dels$End, 1,
                              col = del_color, border = NA)
                 }
                 #Plot SNPs (black vertical lines)
                 snps <- chr_data[chr_data$Type == "SNP",]
                 if(nrow(snps) > 0) {
                   circos.segments(snps$Start, 0, snps$Start, 1,
                                  col = snp_color, lwd = 1)
                 }
               })
  
  #Track label (inside frame)
  circos.text(-5, 0.5, paste0("R272-", track_num),
             sector.index = chromosomes[1],
             track.index = get.cell.meta.data("track.index"),
             facing = "downward", cex = 0.7)}

#Plot combined tracks
plot_combined_track(R272_1, 1)
plot_combined_track(R272_3, 3)
plot_combined_track(R272_5, 5)

#legend
legend("bottomright", 
       legend = c("Deletions", "SNPs"),
       fill = c(del_color, NA),
       border = c(NA, NA),
       lty = c(NA, 1),
       col = c(NA, snp_color),
       lwd = c(NA, 1),
       bty = "n", cex = 0.8)

circos.clear()

#dev.off()


########R1276#######

#svg(width = 10, height = 10, file = "/Users/lb808/Documents/TempertonLab/EllieTong/circos_plot/PAO1_276/R1276_snp_del.svg")

#Set global circos parameters
circos.par(gap.degree = 2, cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 0.01), 
           start.degree = 90, canvas.xlim = c(-1.1, 1.1), canvas.ylim = c(-1.1, 1.1))

bg_color <- "#f7f7f7"     #Background color
ins_color <- "red"    #Light gray for deletions
snp_color <- "#000000"    #Black for SNPs
track_border <- "#7f7f7f"

#Load reference genes from BED file
R1276_genes <- read.table("/Users/et491/OneDrive/Documents/MbyRes/Papers/LPS paper/EllieT_circos/R276.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(R1276_genes) <- c("chr", "start", "end")

#Load SNP/DEL data
load_snp_data <- function(file) {
  dt <- fread(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(dt) <- c("Chrom", "Start", "End", "Type")
  return(dt)
}

R1276_1 <- load_snp_data("/Users/et491/OneDrive/Documents/MbyRes/Papers/LPS paper/EllieT_circos/PAO1_R1276_1_cut.snp.txt")
R1276_2 <- load_snp_data("/Users/et491/OneDrive/Documents/MbyRes/Papers/LPS paper/EllieT_circos/PAO1_R1276_2_cut.snp.txt")
R1276_4 <- load_snp_data("/Users/et491/OneDrive/Documents/MbyRes/Papers/LPS paper/EllieT_circos/PAO1_R1276_4_cut.snp.txt")

#Combine all data to get proper xlim
all_data <- rbind(R1276_1, R1276_2, R1276_4)
chromosomes <- unique(all_data$Chrom)
max_pos <- tapply(all_data$End, all_data$Chrom, max)

#find the max position across all chromosomes
max_length <- 6278719 #I set manually this to the max length of the genomes

#Circos parameters with top orientation (0 on top)
circos.par("start.degree" = 90,  #This puts zero at the top
           "track.height" = 0.15,
           "gap.degree" = 5,
           "cell.padding" = c(0, 0, 0, 0))

#initialize circos with this maximum length for all chromosomes
circos.initialize(factors = chromosomes, 
                 xlim = cbind(rep(0, length(chromosomes)), 
                            rep(max_length, length(chromosomes))))
#Add axis
circos.track(ylim = c(0, 1), bg.col = bg_color, bg.border = NA, 
             panel.fun = function(x, y) {
               chr = CELL_META$sector.index
               xlim = CELL_META$xlim
               #Create breaks every 500kb
               axis_breaks <- seq(0, max_length, by = 5e5)
               #Ensure the last break is exactly at max_length
               if(max(axis_breaks) < max_length) {
                 axis_breaks <- c(axis_breaks, max_length)
               }
               #Create labels (show every 1Mb + the max length)
               axis_labels <- ifelse(axis_breaks %% 1e6 == 0 | axis_breaks == max_length,
                                    paste0(round(axis_breaks/1e6, 1), "Mb"), 
                                    "")
               circos.axis(h = "top", 
                          major.at = axis_breaks,
                          labels = axis_labels, 
                          labels.cex = 1.2,  #Increase here for the axis fonts
                          major.tick.length = 0.2)
               circos.text(mean(xlim), 1.2, chr, cex = 0.8, niceFacing = TRUE)
             }, track.height = 0.08)

#combine deletions and SNPs.
plot_combined_track <- function(data, track_num) {
  circos.track(factors = data$Chrom, ylim = c(0, 1), 
               bg.border = track_border, bg.col = bg_color,
               panel.fun = function(x, y) {
                 current_chr <- CELL_META$sector.index
                 chr_data <- data[data$Chrom == current_chr,]
                 #Plot deletions (gray blocks)
                 dels <- chr_data[chr_data$Type == "INS",]
                 if(nrow(dels) > 0) {
                   circos.rect(dels$Start, 0, dels$End, 1,
                              col = ins_color, border = NA)
                 }
                 #Plot SNPs (black vertical lines)
                 snps <- chr_data[chr_data$Type == "SNP",]
                 if(nrow(snps) > 0) {
                   circos.segments(snps$Start, 0, snps$Start, 1,
                                  col = snp_color, lwd = 1)
                 }
               })
  
  #Track label (inside frame)
  circos.text(-5, 0.5, paste0("R276-", track_num),
             sector.index = chromosomes[1],
             track.index = get.cell.meta.data("track.index"),
             facing = "downward", cex = 0.7)}

#Plot combined tracks
plot_combined_track(R1276_1, 1)
plot_combined_track(R1276_2, 2)
plot_combined_track(R1276_4, 4)

#legend
legend("bottomright", 
       legend = c("Insrtions", "SNPs"),
       fill = c(ins_color, NA),
       border = c(NA, NA),
       lty = c(NA, 1),
       col = c(NA, snp_color),
       lwd = c(NA, 1),
       bty = "n", cex = 0.8)

circos.clear()

#dev.off()
