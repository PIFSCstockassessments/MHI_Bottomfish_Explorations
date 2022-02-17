library(tidyverse)

drs <- read.csv("./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/Data for Meg/D7DealerData2000-2021.csv")
head(drs)

## species codes and common names for deep 7 species
com_names <- c("Hapuupuu","Kalekale","Opakapaka","Ehu","Onaga","Ehu","Lehi","Gindai")
frs_id <- data.frame(COMMON_NAME = com_names,  SPECIES = c(15,17,19,21,22,36,58,97))

## species code and common name for deep 7 plus tunas, billfish, sharks, etc. From Yau 2018, Table 12.
sp_code <- read.csv("./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/frs_species_code.csv")

deep7_drs <- drs %>% 
 mutate(common_name = str_replace_all(local_name, 
                                      c("'Ula'ula Koa'e, Onaga, Ula" = "Onaga",
                                      "'Opakapaka" = "Opakapaka",
                                      "'Ula'ula, Ehu" = "Ehu", 
                                      "HÄ\u0081pu'u, HÄ\u0081pu'upu'u, Shapon" = "Hapuupuu",
                                      "Kalekale, Kalikali" = "Kalekale",
                                      "'Ukikiki, Gindai, Tai" = "Gindai"))) %>% 
  filter(pieces_bought == 1) %>% 
  separate(report_date, into = c("Year", "Month", "Day"), sep = "-", convert = TRUE) %>% 
  mutate(FYear = ifelse(Month < 7, Year, Year + 1)) %>% 
  select(c(Year, 
           FYear,
           Month,
           common_name, 
           pieces_bought, 
           pounds_bought, 
           total_value,
           data_source,
           report_type, 
           dealer_number)) %>% 
  filter(common_name == "Ehu" & pounds_bought < 12 
         | common_name == "Gindai" & pounds_bought < 5 
         | common_name == "Hapuupuu" & pounds_bought < 70 
         | common_name == "Kalekale" & pounds_bought < 4 
         | common_name == "Lehi" & pounds_bought < 33 
         | common_name == "Onaga" & pounds_bought < 35 
         | common_name == "Opakapaka" & pounds_bought < 19) %>% 
  filter(pounds_bought > 0) %>% 
  rename("COMMON_NAME" = "common_name")
          

## Save tidy dataset
write.csv(deep7_drs, file = "./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/Data for Meg/DRS_deep7.csv")



