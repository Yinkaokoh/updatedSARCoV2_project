
# SARS-CoV2 Project in collaboration with Nidia, Elisha,
#Nicholas et al on SARS-CoV2 in Africa
# Script written by OKOH Olayinka Sunday
# 15th January, 2021

#Dataset for lineages were downloaded from GISAID sent to me by Nidia
# Population data was loaded from worldometer.com

###Load the packages that will be needed for the project
library(maptools)
library(RColorBrewer)
library(maps)
library(mapdata)
library(readxl)
library(ggplot2)
library(qwraps2)
library(dplyr)
library(gridExtra)
library(ggcorrplot)
library(ggpubr)
library(gridExtra)
library(vcd)
library(tidyr)
library(maps)
library(mapdata)
library(scatterpie)
library(ggmap)
library(mapproj)
library(readr) #to read in tsv data
#library(read)
#raw_data <- read_tsv(file.choose()) #This could have been used interactively. But to avoid manual selection each time
#read in the data
raw_data<- read_tsv("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/metadata_2021-01-06_18-26.tsv.zip", col_names = TRUE,
                    col_types = cols(strain ="f", virus = "c", gisaid_epi_isl = "-",
                                     genbank_accession = "-", date = "D", region = "f",
                                     country = "c",  division = "-", location = "-",
                                     region_exposure = "-", country_exposure = "-", division_exposure = "-",
                                     segment = "-", host = "f", age = "i", sex = "-",
                                     Nextstrain_clade = "f", pangolin_lineage = "f", GISAID_clade = "f",
                                     originating_lab = "-", submitting_lab = "-", authors = "-", url = "-", title = "-",
                                     paper_url = "-",date_submitted = "D" 
                    ))
#the above is equivalent to cases and lineages in the previous script

#read in a file containing some data specifically for Africa like population, deaths, recovery etc

africa_ncdc <- read_excel("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/african__data_ncdc.xlsx")
#this is eqivalent to africa in the previous script

#convert the ncdc data to data frame
africa_ncdc <- as.data.frame(africa_ncdc)
africa_ncdc$region <- as.character(africa_ncdc$region)

#Add a column name country and remove the region. This is to avoid confusion, as the world map has continents as region
africa_ncdc$country <- africa_ncdc$region



#make some calculations
africa_ncdc <- africa_ncdc %>% mutate(popn = (Population/1000000), cases_per100H = (SARSCoV2_Cases/Population)*100000, tests_per100H = (T_Tests/Population)*100000, deaths_per100cases = (Total_Deaths/SARSCoV2_Cases)*100, recovery_per100cases = (T_Recovered/SARSCoV2_Cases)*100, positive_per100tests = (SARSCoV2_Cases/T_Tests)*100)

#subset the african map from the world map
africa_map <- map_data("world", region = africa_ncdc$region)

#Add map data to the ncdc data
africa_ncdc <- left_join(africa_map, africa_ncdc, by = "region")

#remove the region column
africa_ncdc <- africa_ncdc %>% select(-(region))


#Convert the sarscov2 from worldometer to data frame
sars_cov2_global <- as.data.frame(raw_data)

#Add number of sequences submitted on GISAID by each country
sequence_country <- sars_cov2_global %>% group_by (country) %>% summarise(sequences = n()) %>% mutate(percentage = round((sequences/sum(sequences))*100))
sars_cov2_global <- left_join(sars_cov2_global, sequence_country, by = "country")

#filter Africa data 
sars_cov2_africa <- sars_cov2_global %>% filter (region == "Africa")

africa_map_details <- left_join (africa_ncdc, sars_cov2_africa, by = "country", all = TRUE)

africa_map_details <- africa_map_details %>% mutate(sequence_per_100cases = (sequences/SARSCoV2_Cases)*100)


#To get number of sequences submitted to GISAID per continent
sequence_continent <- sars_cov2_global %>% group_by (region) %>% summarise(sequences = n()) %>% mutate(percentage = round((sequences/sum(sequences))*100))
next_strain_clade_continent <- sars_cov2_global %>% group_by(region, Nextstrain_clade) %>% summarise(Number = n()) %>% mutate(percentage = round((Number/sum(Number))*100))
pangolin_lineage_continent <- sars_cov2_global %>% group_by(region, pangolin_lineage) %>% summarise(N = n()) %>% mutate(percentage = round((N/sum(N))*100))
GISAID_clade_continent <- sars_cov2_global %>% group_by(region, GISAID_clade) %>% summarise(N = n()) %>% mutate(percentage = round((N/sum(N))*100))
#Add continent population
sequence_continent <- as.data.frame(sequence_continent)

#To get the number of cases in Each country, group by country
next_strain_clade_country <- sars_cov2_global %>% group_by(country, Nextstrain_clade)  %>% summarise(Number = n()) %>% mutate(percentage = round((Number/sum(Number))*100))
sequence_country <- sars_cov2_global %>% group_by (country) %>% summarise(sequences = n()) %>% mutate(percentage = round((sequences/sum(sequences))*100)) 
pangolin_lineage_country <- sars_cov2_global %>% group_by(country, pangolin_lineage) %>% summarise(N = n()) %>% mutate(percentage = round((N/sum(N))*100))
GISAID_clade_country <- sars_cov2_global %>% group_by(country, GISAID_clade) %>% summarise(N = n())  %>% mutate(percentage = round((N/sum(N))*100))

#Add number of sequences to sars-cov2_global data
sars_cov2_global <- left_join(sars_cov2_global, sequence_country, by = "country")
#To get data for each continent
#Africa
africa_sequences <- sars_cov2_global %>% filter(region == "Africa") %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
africa_nstrain_clade <- sars_cov2_global %>% filter(region == "Africa") %>% group_by(Nextstrain_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N))* 100)
africa_pangolin_lineage <- sars_cov2_global %>% filter(region == "Africa") %>% group_by(pangolin_lineage) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
africa_GISAID_clade <- sars_cov2_global %>% filter(region == "Africa") %>% group_by(GISAID_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)

#Asia
asia_sequences <- sars_cov2_global %>% filter(region == "Asia") %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
asia_nstrain_clade <- sars_cov2_global %>% filter(region == "Asia") %>% group_by(Nextstrain_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
asia_pangolin_lineage <- sars_cov2_global %>% filter(region == "Asia") %>% group_by(pangolin_lineage) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
asia_GISAID_clade <- sars_cov2_global %>% filter(region == "Asia") %>% group_by(GISAID_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N))* 100)

#Europe
europe_sequences <- sars_cov2_global %>% filter(region == "Europe") %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
europe_nstrain_clade <- sars_cov2_global %>% filter(region == "Europe") %>% group_by(Nextstrain_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N))* 100)
europe_pangolin_lineage <- sars_cov2_global %>% filter(region == "Europe") %>% group_by(pangolin_lineage) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
europe_GISAID_clade <- sars_cov2_global %>% filter(region == "Europe") %>% group_by(GISAID_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)


#North_America
n_america_sequences <- sars_cov2_global %>% filter(region == "North America") %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
n_america_nstrain_clade <- sars_cov2_global %>% filter(region == "North America") %>% group_by(Nextstrain_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
n_america_pangolin_lineage <- sars_cov2_global %>% filter(region == "North America") %>% group_by(pangolin_lineage) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
n_america_GISAID_clade <- sars_cov2_global %>% filter(region == "North America") %>% group_by(GISAID_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)


#South_America
s_america_sequences <- sars_cov2_global %>% filter(region == "South America") %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N))* 100)
s_america_nstrain_clade <- sars_cov2_global %>% filter(region == "South America") %>% group_by(Nextstrain_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
s_america_pangolin_lineage <- sars_cov2_global %>% filter(region == "South America") %>% group_by(pangolin_lineage) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N))* 100)
s_america_GISAID_clade <- sars_cov2_global %>% filter(region == "South America") %>% group_by(GISAID_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)

#Oceania
oceania_sequences <- sars_cov2_global %>% filter(region == "Oceania") %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N))* 100)
oceania_nstrain_clade <- sars_cov2_global %>% filter(region == "Oceania") %>% group_by(Nextstrain_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
oceania_pangolin_lineage <- sars_cov2_global %>% filter(region == "Oceania") %>% group_by(pangolin_lineage) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)
oceania_GISAID_clade <- sars_cov2_global %>% filter(region == "Oceania") %>% group_by(GISAID_clade) %>% summarise(N = n())  %>% mutate(percentage = (N/sum(N)) * 100)



quantile(pangolin_lineage_continent$N, probs = c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99))
# Quantile  10% 20% 25% 30%   40%   50%   60%   70%   75%   80%   90%   95%  .... 99%
# Value     1.0 2.0 2.0 3.0   6.0   14.0  26.0  47.0  63.0  83.0  205.9 509.2     3213.9
# Sequences 315
lineage_10q <- pangolin_lineage_continent %>% filter(N <= 1.0)
lineage_20q <- pangolin_lineage_continent %>% filter(N <= 2.0) 
lineage_top_10 <- pangolin_lineage_continent %>% filter(N >= 205.9)
lineage_top_1 <- pangolin_lineage_continent %>% filter(N >= 3213.9)
lineage_top_2 <- pangolin_lineage_continent %>% filter(N >= 1238.06)

aspect_ratio <- 2.5
######## Charts at the Continental levels
#Figure 1
#Bar chart (pdf)
sequences_continent_bar <- ggplot(sequence_continent, aes(y = reorder(region, sequences), x = sequences))  +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(percentage, "%")), position = position_dodge(width = 1),vjust = 0.5) +
  labs(title = "SARS-CoV-2 sequences submitted on GISAID", y = "Continents", x = "Number of sequences submitted", caption = "Source: GISAID" ) +
  theme(plot.title = element_text(size = 15))
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Continent_sequences_bar.pdf", height = 120, width = 200, units = 'mm')

#Bar chart (jpg)
sequences_continent_bar <- ggplot(sequence_continent, aes(y = reorder(region, sequences), x = sequences))  +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(percentage, "%")), position = position_dodge(width = 1),vjust = 0.5) +
  labs(title = "SARS-CoV-2 sequences submitted on GISAID", y = "Continents", x = "Number of sequences submitted", caption = "Source: GISAID" ) +
  theme(plot.title = element_text(size = 15))
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Continent_sequences_bar.jpg", height = 120, width = 200, units = 'mm')


####NextStrain Chart
#Figure 2A
fig2a <- ggplot(next_strain_clade_continent, aes(y = Nextstrain_clade, x = Number)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(y = "Next Strain Clade", x = "Number of sequences") +
  facet_wrap(~ region) 
#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_nxtstrain_bar.pdf")


####NestStrain for Africa
#Figure 2B
africa_nstrain_clade <- na.omit(africa_nstrain_clade)
fig2b <- ggplot(africa_nstrain_clade, aes(y = reorder(Nextstrain_clade, N), x = N))  +
  geom_bar(stat = "identity", fill = 'blue') +
  geom_text(aes(label = paste0(round(percentage), "%")), position = position_dodge(width = 1),vjust = 0.5) +
  labs(y = "Next Strain Clades", x = "Number of sequences submitted") 
#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/africa_nxtrain_bar.pdf")
#Figure 2
pdf("updated_output/Clade_diversity.pdf")
grid.arrange(fig2a, fig2b, nrow = 2, top = "NextStrain diversity",
             bottom = "Source: GISAID")
dev.off()

jpeg("updated_output/Clade_diversity.jpg")
grid.arrange(fig2a, fig2b, nrow = 2, top = "NextStrain diversity",
             bottom = "Source: GISAID")
dev.off()
##Top 1% Pango lineage diversity
#Figure 3
lineage_diversity_top_1 <- ggplot(lineage_top_1, aes(y = reorder(pangolin_lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  #geom_text(aes(label = paste0(round(percentage), "%")), position = position_dodge(width = 1),vjust = 0.5) +
  labs(title = "Top 1% lineage diversity in the continents", y = "Pango lineage", x = "Number of sequences", caption = "Source: GISAID") +
  facet_wrap(~ region) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Continent_top1_lineage_bar.pdf", height = 120, width = 200, units = 'mm' )


lineage_diversity_top_1 <- ggplot(lineage_top_1, aes(y = reorder(pangolin_lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  #geom_text(aes(label = paste0(round(percentage), "%")), position = position_dodge(width = 1),vjust = 0.5) +
  labs(title = "Top 1% lineage diversity in the continents", y = "Pango lineage", x = "Number of sequences", caption = "Source: GISAID") +
  facet_wrap(~ region) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Continent_top1_lineage_bar.jpg", height = 120, width = 200, units = 'mm' )

##Top Pango lineage diversity in Africa
#Figure 4
#Get the top lineage in Africa
################################
quantile(africa_pangolin_lineage$N, probs = c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99))
# Quantile  10% 20% 25% 30%   40%   50%   60%   70%   75%   80%   90%   95%   96%     97%       98%     99%
# Value     1.0 1.0 2.0 2.0   3.0   5.0  10.0  18.10  22.0  29.0  67.0 179.35 214.12  326.38    486.12  533.96
# Sequences 315
afr_lineage_10q <- africa_pangolin_lineage %>% filter(N <= 1.0)
afr_lineage_20q <- africa_pangolin_lineage %>% filter(N <= 1.0) 
afr_lineage_top_10 <- africa_pangolin_lineage %>% filter(N >= 67.0)
afr_lineage_top_1 <- africa_pangolin_lineage %>% filter(N >= 533.9)
afr_lineage_top_2 <- africa_pangolin_lineage %>% filter(N >= 486.12)
#top 5 lineages in Africa
afr_lineage_top_5 <- africa_pangolin_lineage %>% filter(N >= 349)

###### Data for proportional area stacked barplot (only time is left to be added)
#Get other lineages not in the top 5%
afr_lineage_below_top_5 <- africa_pangolin_lineage %>% filter(N < 179.35)

#Get the sum of the lineages not in the top 5%
other_lineages_sum <- sum(afr_lineage_below_top_5$N)

#Select lineage and N only for both top 5 and others
top5 <- afr_lineage_top_5 %>% select(pangolin_lineage, N)
below_top5 <- data.frame(pangolin_lineage = "others", N = other_lineages_sum)

##join both top 5 and below top 5
top_and_below_5 <- rbind(top5, below_top5)
#Add a percentage column to the above
top_and_below_5 <- top_and_below_5 %>% mutate(percentage = (N/sum(N)*100))
##Bar chart of Top 10% lineages in Africa
##Figure 
fig4a <- ggplot(afr_lineage_top_10, aes(y = reorder(pangolin_lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(round(percentage), "%")), position = position_dodge(width = 1),vjust = 0.5) +
  #labs(title = "Top 10% lineage diversity in Africa", y = "Pango lineage", x = "Number of sequences", caption = "Source: GISAID") +
  labs(y = "Pango lineage", x = "Number of sequences")

  #ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/africa_top10_lineage_bar.pdf")


##Proportional Stacked chart of Top 5% lineages in Africa
#####Proportional Stacked chart
df <- as.data.frame(raw_data)
afr_data <- na.omit(df) %>% filter (region == "Africa") %>% select(date, pangolin_lineage) %>% filter(pangolin_lineage == c("B.1.5", "B.1", "B.1.1", "B.1.1.206", "B.1.351")) %>% group_by(date, pangolin_lineage) %>% summarise(N = n()) %>% mutate(percentage = N/sum(N))

#plot the graph 
#Stack chart area1
fig4b <- ggplot(afr_data, aes(x = date, y = percentage, fill = pangolin_lineage)) +
  geom_area(size = 0.3, color = 'white') +
  #labs(title = "The trend of the top 5 circulating SARS-CoV-2 lineage in Africa", y = "Percentage", x = "Date") + 
labs(y = "Percentage", x = "Date") 
#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/fig4c_africa_prop_stacked_area1.pdf")


#stacked chart bar1
fig4c <- ggplot(afr_data, aes(x = date, y = percentage, fill = pangolin_lineage)) +
  geom_bar(stat = 'identity', position = 'stack')
labs(title = "The trend of the top 5 circulating SARS-CoV-2 lineage in Africa", y = "Percentage", x = "Date") 
#labs(y = "Percentage", x = "Date") 

#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/africa_prop_stacked_bar1.pdf")
#Figure 4 charts Lineage Diversity in Africa
pdf("updated_output/Pango_lineage.pdf")
grid.arrange(fig4a, fig4b, nrow = 2, top = "Pango lineage diversity in Africa", bottom = textGrob("Source: GISAID"))
dev.off()

jpeg("updated_output/Pango_lineage.jpg")
grid.arrange(fig4a, fig4b, nrow = 2, top = "Pango lineage diversity in Africa", bottom = textGrob("Source: GISAID"))
dev.off()

###############################
##To calculate global average et al
#Group the sars-cov2 data
continental_table <- na.omit(sars_cov2_global) %>% group_by(region) %>% summarise(sequences = sum(sequences))
#add population
continental_table$population <- c(4641054775, 1340598147, 747636026, 592072212, 430759766, 43111704)
#global mean per 100H
continental_table <- continental_table %>% mutate(per100H = (sequences/population)*100000)
global_mean_per100H <-  (sum(continental_table$sequences)/sum(continental_table$population))*100000
continental_table <- as.data.frame(continental_table)
##chart of number of persons infected per 100H

Figure5_absolute_cases <- ggplot(continental_table, aes(x = reorder(region, sequences), y = sequences))  +
  geom_bar(stat = "identity", fill = 'orange') +
  geom_text(aes(label = round(sequences)), position = position_dodge(width = 1),vjust = 0.5) +
  labs(title = "Absolute number of reported COVID-19 cases across the continent", y = "Number of cases", x = "Continents", caption = "Source: worldometer" ) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Absolute_continents.pdf", height = 120, width = 200, units = 'mm')

Figure5_absolute_cases <- ggplot(continental_table, aes(x = reorder(region, sequences), y = sequences))  +
  geom_bar(stat = "identity", fill = 'orange') +
  geom_text(aes(label = round(sequences)), position = position_dodge(width = 1),vjust = 0.5) +
  labs(title = "Absolute number of reported COVID-19 cases across the continent", y = "Number of cases", x = "Continents", caption = "Source: worldometer" ) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Absolute_continents.jpg", height = 120, width = 200, units = 'mm')


Figure6_cases_per100H <- ggplot(continental_table, aes(x = reorder(region, per100H), y = per100H))  +
  geom_bar(stat = "identity", fill = 'orange') +
  geom_text(aes(label = round(per100H)), position = position_dodge(width = 1),vjust = 0.5) +
  labs(title = "Number of reported COVID-19 cases per 100,000 population across the continents", y = "Number of cases", x = "Continents", caption = "Source: worldometer" ) +
  geom_hline(yintercept = global_mean_per100H, col = 'red', size = 2)

#theme(plot.title = element_text(size = 15))
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Cases_per100H_continents.pdf", height = 120, width = 200, units = 'mm')

Figure6_cases_per100H <- ggplot(continental_table, aes(x = reorder(region, per100H), y = per100H))  +
  geom_bar(stat = "identity", fill = 'orange') +
  geom_text(aes(label = round(per100H)), position = position_dodge(width = 1),vjust = 0.5) +
  labs(title = "Number of reported COVID-19 cases per 100,000 population across the continents", y = "Number of cases", x = "Continents", caption = "Source: worldometer" ) +
  geom_hline(yintercept = global_mean_per100H, col = 'red', size = 2)

#theme(plot.title = element_text(size = 15))
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Cases_per100H_continents.jpg", height = 120, width = 200, units = 'mm')




fig7a <- ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = SARSCoV2_Cases), color = 'black') +
  coord_map("bonne", parameters = 45) +
  #labs(title = "Absolute number of COVID-19 reported cases in African countries", caption = "Source: worldometer") +
  #labs(caption = "Source: worldometer") +
  scale_fill_continuous(name = "Cases", low = "white", high = "red", 
                        limits = c(0, 800000), breaks = c(0, 10000, 100000, 250000, 500000, 750000)) +
  theme_void()
#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_cases_Map.pdf")   

fig7b <- ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = cases_per100H), color = 'black') +
  coord_map("bonne", parameters = 45) +
  #labs(title = "COVID-19 reported cases per 100,000 \n population in African countries", caption = "Source: worldometer") +
  labs(caption = "Source: worldometer") +
  scale_fill_continuous(name = "Cases", low = "white", high = "red", 
                        limits = c(0, 1500), breaks = c(0, 500, 1000, 1500)) +
  theme_void()

#Figure 7
pdf("updated_output/Cases_in_africa_map.pdf")
grid.arrange(fig7a, fig7b, nrow = 1, top ="COVID-19 cases reported in African countries")
dev.off()

jpeg("updated_output/Cases_in_africa_map.jpg")
grid.arrange(fig7a, fig7b, nrow = 1, top ="COVID-19 cases reported in African countries")
dev.off()

#Africa Map showing SARS-CoV2 sequences submitted in GISAID from various African countries.
#MAP
#Figure 8
fig8a <- ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = sequences), color = 'black') +
  coord_map("bonne", parameters = 45) +
  #labs(title = "SAR-CoV-2 Sequences submitted to GISAID from African countries", caption = "Source: GISAID") +
  scale_fill_continuous(name = "Sequences", low = "white", high = "green", 
                        limits = c(0, 3000), breaks = c(0, 500, 1000, 2000, 3000)) +
  theme_void()
#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_sequences_on_GISAID_Map.pdf")   


########Sequences submitted per 1000 Cases reported
#Figure 6B
fig8b <- ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = sequence_per_100cases*10), color = 'black') +
  coord_map("bonne", parameters = 45) +
  #labs(title = "SARS-CoV2 sequences submitted per 1000 cases reported", caption = "Source: GISAID and worldometer") +
  scale_fill_continuous(name = "Sequences", low = "white", high = "green", 
                        limits = c(0, 125), breaks = c(0, 25, 50, 75, 100, 125)) +
  theme_void()
#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_seqs_per100_cases_GISAID_Map.pdf")   

#Figure 8
pdf("updated_output/Africa_sequeunces_submitted_map.pdf")
grid.arrange(fig8a, fig8b, nrow = 1, top = "Sequences submitted to GISAID from African countries", bottom = "Source: GISAID and worldometer")
dev.off()

jpeg("updated_output/Africa_sequeunces_submitted_map.jpg")
grid.arrange(fig8a, fig8b, nrow = 1, top = "Sequences submitted to GISAID from African countries", bottom = "Source: GISAID and worldometer")
dev.off()

###Deaths
#Figure 9

##### Total Deaths from SARS-CoV2
fig9a <- ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = Total_Deaths), color = 'black') +
  coord_map("bonne", parameters = 45) +
  #labs(title = "Deaths from SARS-CoV-2 in Africa", caption = "Source: worldometer") +
  scale_fill_continuous(name = "Deaths", low = "white", high = "red", 
                        limits = c(0, 30000), breaks = c(0, 500, 5000, 10000, 20000, 30000)) +
  theme_void()
#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_Deaths_Map.pdf")   

#figure 9b Deaths per 1000 cases
fig9b <- ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = deaths_per100cases*10), color = 'black') +
  coord_map("bonne", parameters = 45) +
  #labs(title = "Deaths per 1000 SARS-CoV-2 cases in Africa", caption = "Source: worldometer") +
  scale_fill_continuous(name = "Deaths", low = "white", high = "red", 
                        limits = c(0, 100), breaks = c(0, 10, 20, 30, 40, 50, 60, 100)) +
  theme_void()
#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_deaths_per100cases_Map.pdf")   
#Figure 9
pdf("updated_output/Africa_reported_deaths_map.pdf")
grid.arrange(fig9a, fig9b, nrow = 1, top = "Reported deaths from SARS-CoV-2 in African countries", bottom = "worldometer")
dev.off()

jpeg("updated_output/Africa_reported_deaths_map.jpg")
grid.arrange(fig9a, fig9b, nrow = 1, top = "Reported deaths from SARS-CoV-2 in African countries", bottom = "worldometer")
dev.off()
######## Total Tests per 100,000 Population

#Figure 10
#Tests
fig10a <- ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = T_Tests), color = 'black') +
  coord_map("bonne", parameters = 45) +
  #labs(title = "SARS-CoV2 tests per 100,000 population in Africa", caption = "Source: worldometer") +
  scale_fill_continuous(name = "Tests", low = "white", high = "blue", 
                        limits = c(0, 6000000), breaks = c(0, 250000, 1000000, 2500000, 5000000)) +
  theme_void()
#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_test_100H_Map.pdf")   


#Figure 10b
#Test per 100H

fig10b <- ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = tests_per100H), color = 'black') +
  coord_map("bonne", parameters = 45) +
  #labs(title = "SARS-CoV2 tests per 100,000 population in Africa", caption = "Source: worldometer") +
  scale_fill_continuous(name = "Tests", low = "white", high = "blue", 
                        limits = c(0, 25000), breaks = c(0, 1000, 5000, 10000, 15000, 20000, 25000)) +
  theme_void()
#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_test_100H_Map.pdf")   

pdf("updated_output/Africa_tests_map.pdf")
grid.arrange(fig10a, fig10b, nrow = 1, top = "SARS-CoV-2 tests conducted in African countries", bottom = "worldometer")
dev.off()

jpeg("updated_output/Africa_tests_map.jpg")
grid.arrange(fig10a, fig10b, nrow = 1, top = "SARS-CoV-2 tests conducted in African countries", bottom = "worldometer")
dev.off()

######## Positive per 1000 Tests
#figure 11
#Positive per 1000

Figure11_Positive_per1000 <- ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = positive_per100tests*10), color = 'black') +
  coord_map("bonne", parameters = 45) +
  labs(title = "Number of positive tests per 1000 SARS-CoV-2 in African countries", caption = "Source: worldometer") +
  #labs(caption = "Source: worldometer") +
  scale_fill_continuous(name = "Positive tests", low = "white", high = "red", 
                        limits = c(0, 250), breaks = c(0, 50, 100, 150, 200, 250)) +
  theme_void()
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_positive_per_1000_Map.pdf")   

Figure11_Positive_per1000 <- ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = positive_per100tests*10), color = 'black') +
  coord_map("bonne", parameters = 45) +
  labs(title = "Number of positive tests per 1000 SARS-CoV-2 in African countries", caption = "Source: worldometer") +
  #labs(caption = "Source: worldometer") +
  scale_fill_continuous(name = "Positive tests", low = "white", high = "red", 
                        limits = c(0, 250), breaks = c(0, 50, 100, 150, 200, 250)) +
  theme_void()
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_positive_per_1000_Map.jpg")   


##To list all the Pango lineages circulating in Africa
write.csv(africa_pangolin_lineage, file = "updated_output/pango_lineages_africa.csv", sep = "")

#################################END of FIGURES in the Manuscript#

########## Recovery per 100 cases
ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = recovery_per100cases), color = 'black') +
  coord_map("bonne", parameters = 45) +
  labs(title = "Recovery per 100 SARS-CoV2 cases", caption = "Source: worldometer and NCDC") +
  scale_fill_continuous(name = "Recovery", low = "white", high = "green", 
                        limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  theme_void()
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_recovery_per100cases_Map.pdf")   



###################################END OF AFRICA SCRIPT###############################
#Supplementary Figures



lineage_2nd_top_1 <- pangolin_lineage_continent %>% filter(N < 3213.9 & N >= 1238.06)
lineage_3rd_top_1 <- pangolin_lineage_continent %>% filter(N < 1238.06 & N >= 849)
lineage_4th_top_1 <- pangolin_lineage_continent %>% filter(N < 849 & N >= 619.6)

#Figure 3B
lineage_2nd_top_1_bar <- ggplot(lineage_2nd_top_1, aes(y = reorder(pangolin_lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Lineage diversity of the 2nd top 1% sequences across the continents", y = "Lineage", x = "Number of sequences", caption = "Source: gisaid") +
  facet_wrap(~ region) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_2nd_top1_lineage_bar.pdf")

#Figure 3C
lineage_3rd_top_1_bar <- ggplot(lineage_3rd_top_1, aes(y = reorder(pangolin_lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Lineage diversity of the 3rd top 1% sequences across the continents", y = "Lineage", x = "Number of sequences", caption = "Source: gisaid") +
  facet_wrap(~ region) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_3rd_top1_lineage_bar.pdf")

#Figure 3D
lineage_4th_top_1_bar <- ggplot(lineage_4th_top_1, aes(y = reorder(pangolin_lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Lineage diversity of the 3rd top 1% sequences across the continents", y = "Lineage", x = "Number of sequences", caption = "Source: gisaid") +
  facet_wrap(~ region) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_4th_top1_lineage_bar.pdf")


#####Proportional Stacked chart
df <- as.data.frame(raw_data)
afr_data <- na.omit(df) %>% filter (region == "Africa") %>% select(date, pangolin_lineage) %>% filter(pangolin_lineage == c("B.1.5", "B.1", "B.1.1", "B.1.1.206", "B.1.351")) %>% group_by(date, pangolin_lineage) %>% summarise(N = n()) %>% mutate(percentage = N/sum(N))

#plot the graph 
#Stack chart area1
ggplot(afr_data, aes(x = date, y = percentage, fill = pangolin_lineage)) +
  geom_area(size = 1) +
  labs(title = "The trend of the top 5 circulating SARS-CoV-2 lineage in Africa", y = "Percentage", x = "Date") 
  ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/africa_prop_stacked_area1.pdf")

#stacked chart area2
ggplot(afr_data, aes(x = date, y = N, fill = pangolin_lineage)) +
  geom_area(size = 1) +
  labs(title = "The trend of the top 5 circulating SARS-CoV-2 lineage in Africa", y = "Number of lineages", x = "Date") +
  ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/africa_prop_stacked_area2.pdf")

#stacked chart bar1
ggplot(afr_data, aes(x = date, y = percentage, fill = pangolin_lineage)) +
  geom_bar(stat = 'identity', position = 'stack')
  labs(title = "The trend of the top 5 circulating SARS-CoV-2 lineage in Africa", y = "Percentage", x = "Date") 
  ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/africa_prop_stacked_bar1.pdf")

#stacked chart bar2
  ggplot(afr_data, aes(x = date, y = N, fill = pangolin_lineage)) +
    geom_bar(position = "stack", stat = 'identity')
  labs(title = "The trend of the top 5 circulating SARS-CoV-2 lineage in Africa", y = "Number of lineages", x = "Date") 
    ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/africa_prop_stacked_bar2.pdf")
########################################################## 
    mean(continent_data$N) # which is 26.14
    continent_data_mean <- continent_data %>% filter (N >= 26.14)
    continent_lineage_above_15 <-continent_data %>% filter (N >= 15)
    continent_lineage_less_15 <-continent_data %>% filter (N < 15)
    
    #Lineages greater than 15 (percentage) across continents
    ggplot(continent_lineage_less_15, aes(y = lineage, x = percentage)) +
      geom_bar(stat = "identity", fill = 'steelblue') +
      labs(title = "SARS-CoV-2 lineages (cases less than 15) distribution across continents", y = "Lineage", x = "Percent (%)", caption = "Source: https://www.worldometers.info/coronavirus") +
      facet_wrap(~ continent)
    ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV2_Project_2021/outputs/Lineages_across_continents_less15_percent_bar.pdf")
    
    
########################################################  
######## filtering data for a specific country
specify_country <- "Mauritius"    
africa_map_details %>% filter(country == specify_country) %>% select(country, SARSCoV2_Cases, cases_per100H,sequences, sequence_per_100cases, Total_Deaths, deaths_per100cases, T_Tests, tests_per100H, positive_per100tests) %>% head() %>% View()
ABC<- africa_map_details %>% select(country, SARSCoV2_Cases, cases_per100H,sequences, sequence_per_100cases, Total_Deaths, deaths_per100cases, T_Tests, tests_per100H, positive_per100tests)
write.csv(ABC, file = "updated_output/Afr_data.csv", sep = "")
write.table(ABC, file = "updated_output/Afr.txt")
#Africa with data on GISAID
afr_with_GISAID_data <- sars_cov2_global %>% filter(region == "Africa") %>% select(country) %>% unique()

##We try to zoom down on March and November to see if there was error
##This is because we had broken lines in the  area.
##However, I discovered it was becasue it is area
#not bar. Area continues to other dates unlike bar
March <- afr_data %>% filter(date < "2020-04-01")
Nov_Dec <- afr_data %>% filter(date > "2020-10-31")

write.table(March, file = "updated_output/March_data.txt", sep =  "")
write.table(Nov_Dec, file = "updated_output/Nov_Dec_data.txt", sep =  "")

write.csv(March, file = "updated_output/March_data.csv", sep =  "")
write.csv(Nov_Dec, file = "updated_output/Nov_Dec_data.csv", sep =  "")



####### Used for writing of methods
## specify the field you want in A
A <- SARSCoV2_Cases
ggplot(africa_map_details, aes(y = reorder(region, SARSCoV2_Cases), x = SARSCoV2_Cases))  +
    geom_bar(stat = "identity", fill = 'orange') +
    geom_text(aes(label = SARSCoV2_Cases), position = position_dodge(width = 1),vjust = 0.5) 
    
  
  
  

### GISAID 
#GISAID_clade_continent <- na.omit(GISAID_clade_continent)
#GISAID_clade_continent_bar <- ggplot(GISAID_clade_continent, aes(y = GISAID_clade, x = N)) +
# geom_bar(stat = "identity", fill = "steelblue") +
#  geom_text(aes(label = paste0(round(percentage), "%")), position = position_dodge(width = 1),vjust = 0.5) +
# labs(title = "GISAID clade diversity in the continents", y = "GISAID Clade", x = "Number of sequences", caption = "Source: gisaid") +
#  facet_wrap(~ region) 
#ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_GISAID_bar.pdf")

### The Pangolin lineages are too many. I decided to get the quantiles:



###### SPATIAL ANALYSIS OF SARS-COV2 IN DIFFERENT AFRICAN COUNTRIES ##########
##Figure 5
####### SARS-CoV2 Cases in Africa

Nigeria_data <- africa_map_details %>% filter(country == "Nigeria") %>% group_by(strain)  %>%select (pangolin_lineage, strain, Nextstrain_clade,  GISAID_clade, length, age, date, date_submitted)
Ghana_data <- africa_map_details %>% filter(country == "Ghana")  %>% group_by(strain) %>% select (pangolin_lineage, strain, Nextstrain_clade,  GISAID_clade, length, age, date, date_submitted)

strains_Nigeria <- unique(Nigeria_data$strain)
strains_Ghana <- unique(Ghana_data$strain)

length(strains_Nigeria) #There 223 sequences submitted from Nigeria in GISAID
length(strains_Ghana) #There 70 sequences submitted from Ghana in GISAID
#Export Ghana and Nigeria data
write.csv(Nigeria_data, file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/new_nigeria_data_gisaid.csv")
write.csv(Ghana_data, file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/new_ghana_data_gisaid.csv")

#Classify the Sequence into Regions/State
write.csv(strains_Nigeria, file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updated_strains_nigeria_gisaid.csv")
write.csv(strains_Ghana, file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updated_strains_ghana_gisaid.csv")

library(stringr)
#use str_count(dataset, "pattern") to indentify the patterns
sum(str_count(strains_Nigeria, "AB")) #enter the state code of each to know the number

##Commit this and note to update later