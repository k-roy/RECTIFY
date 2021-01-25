library(tidyverse)
library(cowplot)

DIR <- '/Volumes/GoogleDrive/My Drive/Scripts/pAmapping/'

annotated_pA_37 <- read_tsv(paste0(DIR, 'GSM2268784_37degC_RT_annotated_pA_sites.txt'))
annotated_pA_42 <- read_tsv(paste0(DIR, 'GSM2268791_42degC_RT_annotated_pA_sites.txt'))

annotated_pA_37_filtered <- annotated_pA_37 %>% filter(reads > 10)
annotated_pA_42_filtered <- annotated_pA_42 %>% filter(reads > 10)

my_comparisons <- list( c(TRUE, FALSE) )

give.n <- function(x){
  return(c(y = Inf, label = length(x)),)
}

p1 <- ggplot(annotated_pA_37_filtered, aes(x = overlaps_ORF, y = A19) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)
p2 <- ggplot(annotated_pA_37_filtered, aes(x = overlaps_ORF, y = A19 + G19) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)
p3 <- ggplot(annotated_pA_37_filtered, aes(x = overlaps_ORF, y = A19 + G19/2) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)
p4 <- ggplot(annotated_pA_37_filtered, aes(x = overlaps_ORF, y = A19 + A6 ) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)
p5 <- ggplot(annotated_pA_37_filtered, aes(x = overlaps_ORF, y = A19 + A6 + G19/2) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)
p6 <- ggplot(annotated_pA_37_filtered, aes(x = overlaps_ORF, y = A19 + A6 + G6/2 + G19/2) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)
ggsave(filename = paste(DIR, 'annotated_pA_37_filtered_AG_richness_analysis.eps'), width = 9, height = 3) 

p7 <- ggplot(annotated_pA_42_filtered, aes(x = overlaps_ORF, y = A19) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)
p8 <- ggplot(annotated_pA_42_filtered, aes(x = overlaps_ORF, y = A19 + G19) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)
p9 <- ggplot(annotated_pA_42_filtered, aes(x = overlaps_ORF, y = A19 + G19/2) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)
p10 <- ggplot(annotated_pA_42_filtered, aes(x = overlaps_ORF, y = A19 + A6 ) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)
p11 <- ggplot(annotated_pA_42_filtered, aes(x = overlaps_ORF, y = A19 + A6 + G19/2) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)
p12 <- ggplot(annotated_pA_42_filtered, aes(x = overlaps_ORF, y = A19 + A6 + G6/2 + G19/2) ) + geom_boxplot() + stat_summary(fun.data = give.n, geom = "text", vjust = 2)

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)
ggsave(filename = paste(DIR, 'annotated_pA_42_filtered_AG_richness_analysis.eps'), width = 9, height = 3)

plot_grid(p4, p5, p10, p11, nrow = 1)
ggsave(filename = paste(DIR, 'annotated_pA_37_vs_42_filtered_AG_richness_analysis.eps'), width = 9, height = 3)

