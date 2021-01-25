
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


######## Charts at the Continental levels
#Pie chart - WORK MORE ON THIS TO ADD THE PERCENTAGE
sequences_continent_pie <- ggplot(sequence_continent, aes(x = "", y = sequences, fill = region))  +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) + 
  #geom_text(aes(label = paste0(round(value*100), "%")), position_stack(vjust = 0.5)) +
  labs(x = NULL, y = NULL, fill = NULL, title = "SARS-CoV2 sequences submitted on GISAID from different continents", caption = "Source: gisaid" ) +
  theme_classic() + theme(axis.line = element_blank(),
                          axis.text = element_blank())
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_sequences_pie.pdf")

#Bar chart
sequences_continent_bar <- ggplot(sequence_continent, aes(y = reorder(region, sequences), x = sequences))  +
  geom_bar(stat = "identity") +
  labs(title = "SARS-CoV2 sequences submitted on GISAID from different continents", y = "Continent", x = "Number of sequences submitted", caption = "Source: gisaid" ) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_sequences_bar.pdf")

####NextStrain Chart
next_strain_clade_continent_bar <- ggplot(next_strain_clade_continent, aes(y = Nextstrain_clade, x = Number)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Next strain clade diversity in the continents", y = "Next Strain Clade", x = "Number of sequences", caption = "Source: gisaid") +
  facet_wrap(~ region) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_nxtstrain_bar.pdf")


### GISAID 
GISAID_clade_continent_bar <- ggplot(GISAID_clade_continent, aes(y = GISAID_clade, x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "GISAID clade diversity in the continents", y = "GISAID Clade", x = "Number of sequences", caption = "Source: gisaid") +
  facet_wrap(~ region) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_GISAID_bar.pdf")

### The Pangolin lineages are too many. I decided to get the quantiles:
quantile(pangolin_lineage_continent$N, probs = c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99))
# Quantile  10% 20% 25% 30%   40%   50%   60%   70%   75%   80%   90%   95%  .... 99%
# Value     1.0 2.0 2.0 3.0   6.0   14.0  26.0  47.0  63.0  83.0  205.9 509.2     3213.9
# Sequences 315
lineage_10q <- pangolin_lineage_continent %>% filter(N <= 1.0)
lineage_20q <- pangolin_lineage_continent %>% filter(N <= 2.0) 
lineage_top_10 <- pangolin_lineage_continent %>% filter(N >= 205.9)
lineage_top_1 <- pangolin_lineage_continent %>% filter(N >= 3213.9)
lineage_top_2 <- pangolin_lineage_continent %>% filter(N >= 1238.06)

lineage_top_1_bar <- ggplot(lineage_top_1, aes(y = reorder(pangolin_lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Top 1% lineage diversity in the continents", y = "Lineage", x = "Number of sequences", caption = "Source: gisaid") +
  facet_wrap(~ region) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_top1_lineage_bar.pdf")

lineage_2nd_top_1 <- pangolin_lineage_continent %>% filter(N < 3213.9 & N >= 1238.06)
lineage_3rd_top_1 <- pangolin_lineage_continent %>% filter(N < 1238.06 & N >= 849)
lineage_4th_top_1 <- pangolin_lineage_continent %>% filter(N < 849 & N >= 619.6)

lineage_2nd_top_1_bar <- ggplot(lineage_2nd_top_1, aes(y = reorder(pangolin_lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Lineage diversity of the 2nd top 1% sequences across the continents", y = "Lineage", x = "Number of sequences", caption = "Source: gisaid") +
  facet_wrap(~ region) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_2nd_top1_lineage_bar.pdf")

#
lineage_3rd_top_1_bar <- ggplot(lineage_3rd_top_1, aes(y = reorder(pangolin_lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Lineage diversity of the 3rd top 1% sequences across the continents", y = "Lineage", x = "Number of sequences", caption = "Source: gisaid") +
  facet_wrap(~ region) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_3rd_top1_lineage_bar.pdf")


lineage_4th_top_1_bar <- ggplot(lineage_4th_top_1, aes(y = reorder(pangolin_lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Lineage diversity of the 3rd top 1% sequences across the continents", y = "Lineage", x = "Number of sequences", caption = "Source: gisaid") +
  facet_wrap(~ region) 
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/continent_4th_top1_lineage_bar.pdf")


###### SPATIAL ANALYSIS OF SARS-COV2 IN DIFFERENT AFRICAN COUNTRIES ##########
#Africa Map showing SARS-CoV2 sequences submitted in GISAID from various African countries.
#MAP
ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = sequences), color = 'black') +
  coord_map("bonne", parameters = 45) +
  labs(title = "SARS-CoV2 sequences submitted in GISAID from Africa", caption = "Source: GISAID") +
  scale_fill_continuous(name = "Sequences", low = "white", high = "green", 
                        limits = c(0, 3000), breaks = c(0, 500, 1000, 2000, 3000)) +
  theme_void()
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_sequences_on_GISAID_Map.pdf")   

########Sequences submitted per 100 Cases reported

ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = sequence_per_100cases), color = 'black') +
  coord_map("bonne", parameters = 45) +
  labs(title = "SARS-CoV2 sequences submitted per 100 cases reported", caption = "Source: GISAID and worldometer") +
  scale_fill_continuous(name = "Sequences", low = "white", high = "green", 
                        limits = c(0, 10), breaks = c(0, 1, 2.5, 5.0, 7.5, 10)) +
  theme_void()
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_seqs_per100_cases_GISAID_Map.pdf")   


####### SARS-CoV2 Cases in Africa
ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = SARSCoV2_Cases), color = 'black') +
  coord_map("bonne", parameters = 45) +
  labs(title = "SARS-CoV2 reported cases in Africa", caption = "Source: worldometer") +
  scale_fill_continuous(name = "Cases", low = "white", high = "red", 
                        limits = c(0, 800000), breaks = c(0, 10000, 200000, 500000)) +
  theme_void()
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_cases_Map.pdf")   


##### Total Deaths from SARS-CoV2
ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = Total_Deaths), color = 'black') +
  coord_map("bonne", parameters = 45) +
  labs(title = "Deaths from SARS-CoV2 in Africa", caption = "Source: worldometer") +
  scale_fill_continuous(name = "Deaths", low = "white", high = "red", 
                        limits = c(0, 30000), breaks = c(0, 500, 5000, 10000, 20000, 30000)) +
  theme_void()
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_Deaths_Map.pdf")   

######## Total Tests per 100,000 Population
ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = tests_per100H), color = 'black') +
  coord_map("bonne", parameters = 45) +
  labs(title = "SARS-CoV2 tests per 100,000 population in Africa", caption = "Source: worldometer") +
  scale_fill_continuous(name = "Tests", low = "white", high = "blue", 
                        limits = c(0, 14000), breaks = c(0, 300, 500, 1000, 5000, 10000, 14000)) +
  theme_void()
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_test_100H_Map.pdf")   

######## Deaths per 100 cases
ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = deaths_per100cases), color = 'black') +
  coord_map("bonne", parameters = 45) +
  labs(title = "Deaths per 100 SARS-CoV2 cases in Africa", caption = "Source: worldometer") +
  scale_fill_continuous(name = "Deaths", low = "white", high = "red", 
                        limits = c(0, 10), breaks = c(0, 1, 2, 3, 4, 5, 6, 10)) +
  theme_void()
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_deaths_per100cases_Map.pdf")   

########## Recovery per 100 cases
ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = recovery_per100cases), color = 'black') +
  coord_map("bonne", parameters = 45) +
  labs(title = "Recovery per 100 SARS-CoV2 cases", caption = "Source: worldometer and NCDC") +
  scale_fill_continuous(name = "Recovery", low = "white", high = "green", 
                        limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  theme_void()
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_recovery_per100cases_Map.pdf")   


######## Positive per 100 Tests
ggplot(africa_map_details) +
  geom_polygon(aes(long, lat, group = group, fill = positive_per100tests), color = 'black') +
  coord_map("bonne", parameters = 45) +
  labs(title = "Number of SARS-CoV2 positive test per 100 persons tested", caption = "Source: worldometer") +
  scale_fill_continuous(name = "Number positive", low = "white", high = "red", 
                        limits = c(0, 20), breaks = c(0, 1, 5, 10, 15, 20)) +
  theme_void()
ggsave(file = "C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/updated_output/Africa_positive_per100tests_Map.pdf")   


###This is an after thought analysis hence we had to use the raw data again
##Plot a mosaic plots of Nextstrain clades
#color
pdf("updated_output/nclade_mosaic.pdf")
mosaicplot(~ region.x + Nextstrain_clade, data = sars_cov2_global_mega,
           shade = TRUE, xlab ="Continent", ylab = "Nextstrain clades", main = "Distribution of Nextstrain clades")
dev.off()

#grey
pdf("updated_output/nclade_mosaic2.pdf")
mosaicplot(~ region.x + Nextstrain_clade, data = sars_cov2_global_mega,
           shade = FALSE, xlab ="Continent", ylab = "Nextstrain clades", main = "Distribution of Nextstrain clades")
dev.off()

## Plot a mosaic plot of GISAID clade
#color
pdf("updated_output/GISAIDclade_mosaic.pdf")
mosaicplot(~ region.x + GISAID_clade, data = sars_cov2_global_mega,
           shade = TRUE, xlab ="Continent", ylab = "GISAID clades", main = "Distribution of GISAID clades")
dev.off()

#grey
pdf("updated_output/GISAIDclade_mosaic2.pdf")
mosaicplot(~ region.x + GISAID_clade, data = sars_cov2_global_mega,
           shade = FALSE, xlab ="Continent", ylab = "GISAID clades", main = "Distribution of GISAID clades")
dev.off()



## Plot a mosaic plot of pangolin lineage
#color
pdf("updated_output/lineage_mosaic.pdf")
mosaicplot(~ region.x + pangolin_lineage, data = sars_cov2_global_mega,
           shade = TRUE, xlab ="Continent", ylab = "Lineage", main = "Distribution of SARS-CoV2 lineages in the continents")
dev.off()

#grey
pdf("updated_output/lineage_mosaic2.pdf")
mosaicplot(~ region.x + pangolin_lineage, data = sars_cov2_global_mega,
           shade = FALSE, xlab ="Continent", ylab = "Lineage", main = "Distribution of SARS-CoV2 lineages in the continents")
dev.off()

#### Plot mosaic plots for GISAID and Nextstrain clades
#color
pdf("updated_output/GISandNclade_mosaic.pdf")
mosaicplot(~ region.x + Nextstrain_clade + GISAID_clade, data = sars_cov2_global_mega,
           shade = TRUE, xlab ="Continent", ylab = "Clades", main = "Distribution of GISAID and Nextstrain clades")
dev.off()

#grey
pdf("updated_output/GISandNclade_mosaic2.pdf")
mosaicplot(~ region.x + Nextstrain_clade + GISAID_clade, data = sars_cov2_global_mega,
           shade = FALSE, xlab ="Continent", ylab = "Clades", main = "Distribution of GISAID and Nextstrain clades")
dev.off()






###################################END OF AFRICA SCRIPT###############################
