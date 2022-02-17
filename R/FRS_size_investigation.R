### FRS length comp investigation  
library(tidyverse)
library(sf)


#FRS data broken up into 3 data sets. Need to combine and filter before using in the FRS_Opaka_Investigation.RMarkdown or the Deep7_FRS_Investigation.Rmd
early <- read.csv("./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/Data for Meg/picdr_112849_fy48_15.csv")
late <- read.csv("./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/Data for Meg/picdr_112849_fy16_18.csv")
d2019 <- read.csv("./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/Data for Meg/picdr_112970_fy19.csv")
areas_frs <- read.csv("./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/Data for Meg/BF_Area_Grid_Reggie.csv",header=T)
areas_frs <- areas_frs %>% 
  filter(Valid. == "") %>% 
  select(-Valid.)

HDAR_areas <- st_read("./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/GIS/Shapefiles/DAR_Reporting_grids_all.shp")
HDAR_areas <- HDAR_areas %>% 
  filter()

frs <- rbind(early, late, d2019)
head(frs)

## species codes and common names for deep 7 species
com_names <- c("Hapuupuu","Kalekale","Opakapaka","Ehu","Onaga","Ehu","Lehi","Gindai")
frs_id <- data.frame(COMMON_NAME = com_names,  SPECIES = c(15,17,19,21,22,36,58,97))


## species code and common name for deep 7 plus tunas, billfish, sharks, etc. From Yau 2018, Table 12.
sp_code <- read.csv("./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/frs_species_code.csv")


deep7_frs <- frs %>% 
  filter(SPECIES %in% frs_id$SPECIES) %>% 
  mutate(SPECIES = replace(SPECIES, 36, 21)) %>% 
  filter(AREA %in% areas_frs$area) %>% 
  filter(SUBAREA != "A|B") %>% 
  filter(AREA != 16123 & !is.na(SUBAREA)) %>% 
  filter(CAUGHT == 1) %>% 
  filter(GEAR == 3) %>% 
  select(SPECIES, HOURS, CAUGHT, LBS, NUM_SOLD, LBS_SOLD, F_SEASON, F_LBS, CYEAR, FYEAR, AREA) %>% 
  left_join(frs_id, by = "SPECIES") %>% 
  filter(COMMON_NAME == "Ehu" & LBS < 12 
         | COMMON_NAME == "Gindai" & LBS < 5 
         | COMMON_NAME == "Hapuupuu" & LBS < 563 
         | COMMON_NAME == "Kalekale" & LBS < 4 
         | COMMON_NAME == "Lehi" & LBS < 33 
         | COMMON_NAME == "Onaga" & LBS < 35 
         | COMMON_NAME == "Opakapaka" & LBS < 19)

#deep7 <- unique(camera_lengths$COMMON_NAME)[-c(8:9)]

deep7_frs$Region <- cut(as.numeric(deep7_frs$AREA), 
                        breaks = c(99,300,400,500,20000), 
                        labels = c("BIG ISLAND","MAUI NUI","OAHU",'KAUAI-NIIHAU'))

write.csv(deep7_frs, file = "./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/Deep_7_FRS.csv")

opaka_frs <- deep7_frs %>% 
  filter(str_detect(COMMON_NAME, "Opakapaka")) %>% 
  mutate(Decade = floor(FYEAR / 10) * 10,
         Decade = factor(Decade, levels = seq(1940, 2010, by = 10)))
  
summary(opaka_frs)


head(frs)
colnames(frs)

frs %>% 
  left_join(sp_code, by = "SPECIES") %>% 
  #mutate(TIMELINK = round_date(as.POSIXct(frs$TIMELINK[1], tz = "HST"), "10 minutes")) %>% 
  group_by(VESSEL, FISHED, common_name) %>% summarise(n())
  select(LICENSE, TIMELINK, VESSEL, FISHED, common_name, LBS, CAUGHT, DEPTH_BEG, DEPTH_END)


library(lubridate)
round_date(as.POSIXct(frs$TIMELINK[1], tz = "HST"), "10 minutes")

summary(frs)
