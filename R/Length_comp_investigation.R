library(tidyverse)
library(magrittr)

camera_lengths <- read.csv("./Data/BFISH 2016-2019 Camera Lengths.csv")
fishing_lengths <- read.csv("./Data/BFISH 2016-2019 Research Fishing Lengths.csv")
PSU <- read.csv("./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/Data/BFISH PSU lookup table.cvs")
head(camera_lengths)
head(fishing_lengths)
head(PSU)

sp_code <- "PRFI"


### What is causing the mode in the smaller sizes of the BFISH camera data?
cam <- camera_lengths %>% 
filter(str_detect(SPECIES_CD, sp_code)) %>% 
filter(str_detect(BFISH, "_S", negate = TRUE)) %>% 
select(-X)
head(cam)

summary(cam)


### Match up PSU with cam and fishing data

PSU %<>% 
  select(-c(Island, 
            Lat_DD_MM, 
            Long_DD_MM, 
            Lat_DD_MM_SS, 
            Long_DD_MM_SS, 
            STRATA_2020, 
            Depth_MEAN_m, 
            Depth_MIN_m, 
            Depth_MAX_m, 
            Depth_STD_m, 
            med_slp, 
            med_acr, 
            Depth_px_n, 
            Slp_px_n, 
            BS_px_nj, 
            BS_px_over_136j,
            BS_pct_over_136j, 
            pctHB, 
            pctHS))


cam <- merge(cam, PSU, by = "PSU")
head(cam)

cam %>% 
ggplot(aes(x = LENGTH_CM)) +
geom_density()

cam %>% 
filter(LENGTH_CM < 20) %>% 
ggplot(aes(x = LENGTH_CM)) +
geom_density()

small_mode_drops <- cam %>% 
filter(LENGTH_CM < 15 & LENGTH_CM > 10) %>% 
group_by(DROP_CD) %>% 
summarise(n = n()) %>% 
filter(n >= 10) %>% 
pull(DROP_CD)

small_mode_drops

cam %>% 
filter(str_detect(DROP_CD, paste0(small_mode_drops, collapse = "|"))) %>% 
filter(LENGTH_CM < 20) %>% 
group_by(DROP_CD) %>% 
summarise(n())

## 4 camera drops had 10 or more fish < 20cm 
## The samples came from Maui Nui (1 site) and Oahu (3 sites)
## Those small fish came from depths of 80 - 110 m (first 1st quantile of depths)
cam %>% 
filter(str_detect(DROP_CD, paste0(small_mode_drops, collapse = "|")))  %>% 
ggplot(aes(x = LENGTH_CM, y = OFFICIAL_DEPTH_M)) + 
geom_point(aes(colour = DROP_CD), size = 3.5) +
scale_y_reverse() +
theme_classic()

### What is causing the mode in the larger sizes of the BFISH camera data?
cam %>% 
filter(LENGTH_CM > 50 & LENGTH_CM < 60) %>% 
ggplot(aes(x = LENGTH_CM)) +
geom_density()

cam %>% 
filter(LENGTH_CM > 50 & LENGTH_CM < 60) %>% 
group_by(DROP_CD) %>% 
summarise(n = n()) %>% 
filter(n > 1) %>% 
ggplot(aes(x = DROP_CD, y = n)) + geom_point()

large_mode_drops <- cam %>% 
filter(LENGTH_CM > 50 & LENGTH_CM < 60) %>% 
group_by(DROP_CD) %>% 
summarise(n = n()) %>% 
filter(n > 5) %>% 
pull(DROP_CD)

cam %>% 
filter(str_detect(DROP_CD, paste0(large_mode_drops, collapse = "|")))  %>% 
ggplot(aes(x = LENGTH_CM, y = OFFICIAL_DEPTH_M)) + 
geom_point(aes(colour = DROP_CD, shape = Island), size = 3.5) +
scale_y_reverse() +
theme_classic()


fish <- fishing_lengths %>% 
  filter(str_detect(SPECIES_CD, sp_code)) %>% 
  filter(str_detect(BFISH, "_S", negate = TRUE)) %>% 
  select(-X)

small_mode_fish <- fish %>% 
filter(LENGTH_CM > 20 & LENGTH_CM < 25)  %>% 
group_by(SAMPLE_ID) %>% 
summarise(n = n()) %>% 
filter(n > 1) %>% 
pull(SAMPLE_ID)

summary(fish)

fish %>% 
filter(str_detect(SAMPLE_ID, paste0(small_mode_fish, collapse = "|")))  %>% 
filter(LENGTH_CM > 20 & LENGTH_CM < 25) %>% 
ggplot(aes(x = LENGTH_CM, y = SAMPLE_MEAN_DEPTH_M)) + 
geom_point(aes(colour = SAMPLE_ID, shape = Island), size = 3.5) +
scale_y_reverse() +
theme_classic()

cam.merge <- cam %>% 
  select(PSU, DROP_CD, SPECIES_CD, OFFICIAL_DEPTH_M, LENGTH_CM, Island) %>% 
  rename(MEAN_DEPTH = "OFFICIAL_DEPTH_M",
         SAMPLE_ID = "DROP_CD") %>% 
  mutate(Year = str_sub(SAMPLE_ID, 1, 4),
         Month = str_sub(SAMPLE_ID, 5,6),
         Day = str_sub(SAMPLE_ID, 7,8),
         Gear = "Camera") 

fish.merge <- fish %>% 
  select(PSU, SAMPLE_ID, SPECIES_CD, SAMPLE_MEAN_DEPTH_M, LENGTH_CM, Island) %>% 
  rename(MEAN_DEPTH = "SAMPLE_MEAN_DEPTH_M") %>% 
  mutate(Year = str_sub(SAMPLE_ID, 1, 4),
    Month = str_sub(SAMPLE_ID, 5, 6),
    Day = str_sub(SAMPLE_ID, 7, 8),
    Gear = "Research Fishing") 

lengths.combo <- bind_rows(cam.merge, fish.merge)

lengths.combo %>% 
  mutate(SAMPLE_ID = factor(SAMPLE_ID)) %>% 
  distinct(SAMPLE_ID, .keep_all = TRUE) %>% 
  group_by(Gear, Year) %>% 
  summarise(N = n()) %>% 
  pivot_wider(names_from = Gear, values_from = N) %>% 
  gt() %>% 
  tab_header(title = "Number of Sampling Events")
  

lengths.combo %>% 
  ggplot(aes(x = LENGTH_CM)) +
  geom_density(aes(colour = Gear)) +
  facet_wrap(~Year) +
  theme_classic()


lengths.combo %>% 
  mutate(SAMPLE_ID = factor(SAMPLE_ID)) %>% 
  distinct(SAMPLE_ID, .keep_all = TRUE) %>% 
  group_by(Gear, Island) %>% 
  summarise(N = n()) %>% 
  pivot_wider(names_from = Gear, values_from = N) %>% 
  gt() %>% 
  tab_header(title = "Number of Sampling Events by Island")

lengths.combo %>% 
  ggplot(aes(LENGTH_CM)) +
  geom_density(aes(colour = Gear)) + 
  facet_wrap(~Island) +
  theme_classic() +
  labs(main = "Length Distributions by Island")


lengths.combo %>% 
  ggplot(aes(x = Year, y = LENGTH_CM, color = Gear)) +
  geom_violin() +
  theme_classic() +
  scale_fill_nmfs(palette = "regional web")

lengths.combo %>% 
  ggplot(aes(LENGTH_CM)) +
  geom_density(aes(colour = Year)) +
  facet_wrap(~Gear) +
  theme_classic()


### Deep 6 Investigations 

deep_6 <- unique(camera_lengths$COMMON_NAME)[-c(1,8:9)]
camera_lengths %>% 
  filter(COMMON_NAME %in% deep_6) %>% 
  distinct(SPECIES_CD, .keep_all = TRUE) %>% 
  select(SPECIES_CD, SCIENTIFIC_NAME, COMMON_NAME)

cam_6 <- camera_lengths %>% 
  filter(str_detect(BFISH, "_S", negate = TRUE)) %>% 
  filter(COMMON_NAME %in% deep_6) %>% 
  separate(BFISH, into = c("BFISH", "Year", "Seas"), sep = "_") %>% 
  rename(SAMPLE_ID = "DROP_CD",
         DEPTH = OFFICIAL_DEPTH_M) %>% 
  select(PSU, 
         SAMPLE_ID, 
         Year, 
         SPECIES_CD, 
         COMMON_NAME, 
         DEPTH, 
         LENGTH_CM, 
         Island) %>% 
  mutate(Gear = "Camera")


cam_6 %>% 
  ggplot(aes(x = LENGTH_CM, color = COMMON_NAME)) +
  geom_density() +
  theme_classic()

cam_6 %>% 
  group_by(COMMON_NAME, Year) %>% 
  summarise(n = n(), 
            mean_length = mean(LENGTH_CM),
            CV = sd(LENGTH_CM)/mean(LENGTH_CM)) %>% 
  gt(groupname_col = "COMMON_NAME",
     rowname_col = "Year") %>% 
  fmt_number(columns = c("mean_length", "CV"), decimals = 2)

cam_6 %>% 
  group_by(COMMON_NAME) %>% 
  summarise(n = n(), 
            mean_length = mean(LENGTH_CM),
            CV = sd(LENGTH_CM)/mean(LENGTH_CM)) %>% 
  arrange(by_group = FALSE, desc(n)) %>% 
  gt() %>% 
  fmt_number(columns = c("mean_length", "CV"), decimals = 2)

cam_6 %>% 
  group_by(COMMON_NAME, Year) %>% 
  summarise(N = n()) %>% 
  pivot_wider(names_from = COMMON_NAME, values_from = N, values_fill = 0) %>% 
  gt() %>% 
  tab_header(title = "Number of Individuals Caught per Year") %>% 
  grand_summary_rows(columns = c("Ehu", "Gindai", "Hapuupuu", "Kalekale", "Lehi", "Onaga"), fns = list(Total = "sum"))


fish_6 <- fishing_lengths %>% 
  filter(str_detect(BFISH, "_S", negate = TRUE)) %>% 
  filter(COMMON_NAME %in% deep_6) %>% 
  separate(BFISH, into = c("BFISH", "Year", "Seas"), sep = "_") %>% 
  rename(DEPTH = SAMPLE_MEAN_DEPTH_M) %>% 
  select(PSU, 
         SAMPLE_ID,
         Year, 
         SPECIES_CD, 
         COMMON_NAME, 
         DEPTH, 
         LENGTH_CM, 
         Island) %>% 
  mutate(Gear = "Research Fishing")

head(cam_6)
head(fish_6)

deep_6_df <- cam_6 %>% 
  bind_rows(fish_6)


deep_6_df %>% 
  ggplot(aes(x = LENGTH_CM)) +
  geom_density(aes(colour = Gear)) +
  geom_text(data = n.df, aes(x = 90, y = .1, label = paste("N = ", N))) +
  facet_wrap(~COMMON_NAME) +
  theme_classic()


n.df <- data.frame(COMMON_NAME = c("Kalekale", "Ehu", "Lehi", "Onaga", "Hapuupuu", "Gindai"), N = c(278, 50, 42, 14, 9, 6))


same_psu <- unique(camera_lengths[which(camera_lengths$PSU %in% fishing_lengths$PSU),2])


cam_depths <- camera_lengths %>% 
  filter(PSU %in% same_psu) %>% 
  select(PSU, OFFICIAL_DEPTH_M)

fish_depths <- fishing_lengths %>% 
  filter(PSU %in% same_psu) %>% 
  select(PSU, SAMPLE_MEAN_DEPTH_M)

left_join(cam_depths, fish_depths, by = "PSU") %>% 
  distinct() %>% 
  mutate(PSU = as.factor(PSU)) %>% 
  ggplot(aes(x = PSU, y = OFFICIAL_DEPTH_M)) +
  geom_point(color = "red") +
  geom_point(aes(x = PSU, y = SAMPLE_MEAN_DEPTH_M)) +
  theme_classic() +
  ggsave(filename = "./FRMD-SAP-MOshima-SS3_Opakapaka_Assessment/PSU_depths.png")

camera_lengths %>% group_by(BFISH) %>% distinct(DROP_CD) %>% summarise(n())

cam %>% 
  filter(str_detect(BFISH, "_F")) %>% 
  separate(BFISH, into = c("BFISH", "Year", "Seas"), sep = "_") %>% 
  filter(Year == 2016) %>% 
  group_by(DROP_CD) %>% 
  ggplot(aes(x = LENGTH_CM, fill = DROP_CD)) +
  #geom_density() 
  geom_histogram(binwidth = 1) 
  facet_wrap(~Year)

  head(cam)
  
  
cam %>% 
  separate(BFISH, into = c("BFISH", "Year", "Seas"), sep = "_") %>%
  filter(STRATA == "HB_L_S" & Year == 2017) %>% 
  group_by(STRATA) %>% 
  summarise(sd(LENGTH_CM))

