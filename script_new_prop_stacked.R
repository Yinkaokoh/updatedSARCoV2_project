library(dplyr)
library(tidyr)
library(ggplot2)
## To plot the stacked proportional bar showing missing dates

Top_lineages <- metadata_2021_01_06_18_26 %>% select(pangolin_lineage, date_submitted)  %>% filter(pangolin_lineage == c("B.1", "B.1.1", "B.1.5", "B.1.351", "B.1.1206")) %>% group_by(date_submitted, pangolin_lineage) %>% summarise (N = n()) %>% mutate(Percentage = N/sum(N))


Top_lineages2 <- Top_lineages

# Specify the dates columns as date format
library(lubridate) # to work with date
library(aweek)

Top_lineages2$date_submitted <- ymd(Top_lineages2$date_submitted)
Top_lineages2$DATE <- ymd(Top_lineages2$date_submitted)


# aggregate days to weeks
Top_lineages2$DATE = week2date(trunc(date2week(Top_lineages2$DATE)))

### ADDED LINES START HERE
Top_lineages2 =  Top_lineages2 %>% select(pangolin_lineage, N, DATE) %>% group_by(DATE, pangolin_lineage) %>% summarise (N = n()) %>%
  mutate(Percentage = N/sum(N))

# Need to add 0's to anchor missing weeks/weeks with 0
Top_lineages2_new = Top_lineages2 %>% group_by(pangolin_lineage) %>% 
  complete(DATE = seq(min(DATE - 7), max(DATE + 7), "week"), fill = list(Value = 0)) 
Top_lineages2_new$N[is.na(Top_lineages2_new$N)] <- 0

Top_lineages2_new$Percentage[is.na(Top_lineages2_new$Percentage)] <- 0
Top_lineages2 = Top_lineages2_new


ggplot(Top_lineages2, 
       aes(x = DATE, 
           y = Percentage, 
           fill = pangolin_lineage)
)+
  geom_area(size = 0.3, 
            color = 'white'
  ) +
  labs(y = "Percentage", x = "Date")
ggsave("new_Top5_incidence3.jpg", height =  120, width = 200, units = "mm")
ggsave("new_Top5_incidence3_bar.pdf", height =  120, width = 200, units = "mm")

### use absolute number, not percentage

ggplot(Top_lineages2, 
       aes(x = DATE, 
           y = N, 
           fill = pangolin_lineage)
)+
  geom_area(size = 0.3, 
            color = 'white'
  ) +
  labs(y = "Percentage", x = "Date")
ggsave("new_Top5_incidence4.jpg", height =  120, width = 200, units = "mm")
ggsave("new_Top5_incidence4_bar.pdf", height =  120, width = 200, units = "mm")

