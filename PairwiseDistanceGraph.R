
library(rio)
library(tidyverse)
library(readr)
library(graphics)
library(sm)
library(vioplot) 
library(gridExtra)
library(ggplot2)


PairwiseDistances <- import("Regions_Within.csv")

unique(PairwiseDistances$Continent)

PairwiseDistances$Continent <- factor(PairwiseDistances$Continent,levels = c("South America","Oceania","North America","Europe","Asia","Africa"))


VP <- ggplot(PairwiseDistances, aes(x=Continent, y=Dist
                                    , color = Continent), scale = "area", trim = FALSE)
VP <- VP + geom_violin(draw_quantiles = c(0.5), adjust = 2.3) 
VP <- VP + scale_color_manual(values=c("#919191","#919191", "#919191", "#919191", "#919191", "#990e02"))
#VP <- VP + scale_fill_manual(values=c("#919191","#919191", "#919191", "#919191", "#919191", "#990e02"))
VP <- VP + labs(x="Continents", y="Genetic Distance (d)") + theme_light() + coord_flip() 
VP <- VP + theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent", colour = NA),
  plot.title = element_text(hjust = 0.5)
)
VP <- VP + theme(legend.position="none")
VP

ggsave("Continent_PairwiseDist_adjustbw2.3.pdf", VP, width = 8, height = 6)

SouthAmerica <- PairwiseDistances$Dist[PairwiseDistances$Continent=="South America"]
Oceania <- PairwiseDistances$Dist[PairwiseDistances$Continent=="Oceania"]
NorthAmerica <- PairwiseDistances$Dist[PairwiseDistances$Continent=="North America"]
Europe <- PairwiseDistances$Dist[PairwiseDistances$Continent=="Europe"]
Asia <- PairwiseDistances$Dist[PairwiseDistances$Continent=="Asia"]
Africa <- PairwiseDistances$Dist[PairwiseDistances$Continent=="Africa"]


wilcox.test(Africa, SouthAmerica, alternative = "greater") 
wilcox.test(Africa, Oceania, alternative = "greater") 
wilcox.test(Africa, NorthAmerica, alternative = "less") 
wilcox.test(Africa, Europe, alternative = "less")
wilcox.test(Africa, Asia, alternative = "less") 





