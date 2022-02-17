library(tidyverse)
library(nmfspalette)

dir <- "C:/Users/Megumi.Oshima/Documents/FRMD-SAP-MOSHIMA-SS3_Opakapaka_Assessment"
#BFISH datasets
camera_lengths <- read.csv(paste0(dir, "/Data/BFISH 2016-2019 Camera Lengths.csv"))
fishing_lengths <- read.csv(paste0(dir, "/Data/BFISH 2016-2019 Research Fishing Lengths.csv"))

#Life-history data
lh <- read.csv("./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/Life-history database.csv")

#FRS data broken up into 3 datasets
early <- read.csv(paste0(dir, "/Data/Data for Meg/picdr_112849_fy48_15.csv"))
late <- read.csv(paste0(dir, "/Data/Data for Meg/picdr_112849_fy16_18.csv"))
d2019 <- read.csv(paste0(dir, "/Data/Data for Meg/picdr_112970_fy19.csv"))
areas_frs <- read.csv(paste0(dir, "/Data/Data for Meg/BF_Area_Grid_Reggie.csv"), header=T)
areas_frs <- areas_frs %>% 
  filter(Valid. == "") %>% 
  select(-Valid.)
frs <- rbind(early, late, d2019)
head(frs)

#DRS data 2000-2021
deep7_drs <- read.csv(paste0(dir, "/Data/Data for Meg/DRS_deep7.csv"))

com_names <- c("Hapuupuu","Kalekale","Opakapaka","Ehu","Onaga","Ehu","Lehi","Gindai")
frs_id <- data.frame(COMMON_NAME = com_names,  SPECIES = c(15,17,19,21,22,36,58,97))
deep7 <- unique(camera_lengths$COMMON_NAME)[-c(8:9)]

deep7_frs <- frs %>% 
  filter(SPECIES %in% frs_id$SPECIES) %>% 
  mutate(SPECIES = replace(SPECIES, 36, 21)) %>% 
  filter(AREA %in% areas_frs$area) %>% 
  filter(SUBAREA != "A|B") %>% 
  filter(AREA != 16123 & !is.na(SUBAREA)) %>% 
  select(SPECIES, HOURS, CAUGHT, LBS, NUM_SOLD, LBS_SOLD, F_SEASON, F_LBS, CYEAR, FYEAR) %>% 
  left_join(frs_id, by = "SPECIES")

cam_sum <- camera_lengths %>% 
  filter(COMMON_NAME %in% deep7) %>% 
  select(COMMON_NAME, OFFICIAL_DEPTH_M, LENGTH_CM, BFISH) %>% 
  mutate(Year = as.integer(str_extract_all(BFISH, "[0-9]+", simplify = TRUE))) %>% 
  group_by(COMMON_NAME, Year) %>% #for first figure remove Year from group_by
  summarise(N = n()) %>% 
  mutate(Data_Source = "BFISH_camera_length")

rfish_sum <- fishing_lengths %>% 
  filter(COMMON_NAME %in% deep7) %>% 
  select(COMMON_NAME, LENGTH_CM, BFISH) %>% 
  mutate(Year = as.integer(str_extract_all(BFISH, "[0-9]+", simplify = TRUE))) %>% 
  group_by(COMMON_NAME, Year) %>% 
  summarise(N = n()) %>% 
  mutate(Data_Source = "BFISH_fishing_length")

frs_sum <- deep7_frs %>% 
  rename(Year = FYEAR) %>% 
  group_by(COMMON_NAME, Year) %>% 
  summarise(N = n()) %>% 
  mutate(Data_Source = "FRS_lbs_caught",
         N = N/100) 

drs_sum <- deep7_drs %>% 
  #rename("COMMON_NAME" = "common_name") %>% 
  group_by(COMMON_NAME, Year) %>% 
  summarise(N = n()) %>% 
  mutate(Data_Source = "DRS_lbs_caught",
         N = N/100)

lh_sum <- lh %>% 
  mutate(COMMON_NAME = str_replace_all(TAXONNAME, 
                                       c("Pristopomoides filamentosus" = "Opakapaka",
                                         "Pristipomoides filamentosus" = "Opakapaka",
                                         "Hyporthodus quernus" = "Hapuupuu",
                                         "Etelis coruscans" = "Onaga",
                                         "Etelis carbunculus" = "Ehu",
                                         "Pristipomoides sieboldii" = "Kalekale",
                                         "Apareus rutilans" = "Lehi",
                                         "Aphareus rutilans" = "Lehi",
                                         "Pristopomoides zonatus" = "Gindai",
                                         "Pristipomoides zonatus" = "Gindai"))) %>% 
  filter(str_detect(LOCATION, "HI|Hawaii")) %>% 
  select(Year, COMMON_NAME, N, N.1, N.2) %>% 
  rename("Age-Length" = "N",
         "Maturity" = "N.1",
         "Length-Weight" = "N.2") %>% 
  pivot_longer(cols = -c(COMMON_NAME, Year), names_to = "Data_Source", values_to = "N") %>% 
  filter(!is.na(N)) %>% 
  group_by(COMMON_NAME, Data_Source, Year) %>% 
  summarise(N = n()) %>% 
  mutate(N = N*5) %>% 
  select(COMMON_NAME, Year, N, Data_Source)
  
# create plot with life history data
data.df <- bind_rows(cam_sum, rfish_sum, frs_sum, drs_sum, lh_sum) %>% 
  mutate(COMMON_NAME = factor(COMMON_NAME, 
                              levels = c("Opakapaka", 
                                         "Onaga", 
                                         "Ehu", 
                                         "Kalekale", 
                                         "Hapuupuu", 
                                         "Gindai", 
                                         "Lehi")),
         Data_Source = factor(Data_Source, levels = c("BFISH_camera_length", 
                                                      "BFISH_fishing_length",
                                                      "FRS_lbs_caught",
                                                      "DRS_lbs_caught",
                                                      "Age-Length",
                                                      "Length-Weight",
                                                      "Maturity"))) %>% 
  group_by(Data_Source) 

# create plot with no life history data
data.df <- bind_rows(cam_sum, rfish_sum, frs_sum, drs_sum) %>% 
  mutate(COMMON_NAME = factor(COMMON_NAME, 
                              levels = c("Opakapaka", 
                                         "Onaga", 
                                         "Ehu", 
                                         "Kalekale", 
                                         "Hapuupuu", 
                                         "Gindai", 
                                         "Lehi")),
         Data_Source = factor(Data_Source, levels = c("BFISH_camera_length", 
                                                      "BFISH_fishing_length",
                                                      "FRS_lbs_caught",
                                                      "DRS_lbs_caught"))) %>% 
  group_by(Data_Source) 

#get x coordinates for text labels to the side
x.pos <- group_indices(data.df)


##
data_sources <- ggplot(data = data.df, aes(x = Data_Source, y = COMMON_NAME, alpha = N)) + 
  geom_point(aes(size = N+5), color = "#0093D0", show.legend = FALSE) +
  geom_text(aes(label = N, x = x.pos + .23), alpha = 1.0, size = 4, color = "#00467F") +
  scale_alpha_continuous(range = c(0.3, 0.7)) + 
  scale_size_area(max_size = 20) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  labs(x = "Data Source", 
       y = "Common Name",
       caption = "Amount of data available per species (Nrows). Note FRS_lbs_caught and DRS_lbs_caught are x100.") +
  theme_classic() 
ggsave(filename = paste0(dir, "/Data/Data_sources_N.png"), plot = data_sources)

  
data_sources_N <- ggplot(data = data.df, aes(x = Year, y = COMMON_NAME, alpha = N, color = Data_Source)) +
  geom_point(aes(size = N*2.5), show.legend = FALSE) +
  scale_alpha_continuous(range = c(0.3, 0.7)) + 
  scale_size_area(max_size = 10) +
  scale_color_nmfs("regional web") +
  theme_bw() +
  facet_wrap(~Data_Source, strip.position = "top", scales = "free_x") +
  theme(panel.spacing = unit(.5, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.line = element_blank(),            # disable axis lines
        axis.title = element_blank(),           # disable axis titles
        panel.border = element_blank(),         # disable panel border
        panel.grid.major.x = element_blank(),   # disable lines in grid on X-axis
        panel.grid.minor.x = element_blank()) +
  labs(caption = "Note FRS_lbs_caught and DRS_lbs_caught are x100.") 

ggsave(filename = paste0(dir, "/Data/Data_Sources_N_Year.png"), plot = data_sources_N, width = 30, height = 15, units = "cm")

