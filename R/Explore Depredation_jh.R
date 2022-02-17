library(tidyverse)
library(lubridate)
library(hrbrthemes)

# load data, lookup tables, and join fields
d1 <- read_csv("Deep 7 Uku 2000-2009.csv")
d2 <- read_csv("Deep 7 Uku 2010-2021.csv")

data <- rbind(d1,d2)
colnames(data)

data <- data %>%
          select(license_number, trip_start_date, trip_end_date, vessel_id, departure_port_code, landing_port_code, charter_trip, comments, day_fished, wind_speed, buoy_code, area_code, area, subarea, method_code, method_name, method_abbreviation, method_type_code, method_type, hours_fished, number_of_gear_units, food_species_code, number_caught, pounds_caught, number_released, number_lost_to_predator, predator)

# predation data collection started - October 1st 2002 and limit to MHI
data <- filter(data, day_fished >= as.Date("2002-10-01"))
                     
data <- rename(data, SPECIES = food_species_code)

data <- data %>%
  mutate(Year = year(day_fished),
         Month = month(day_fished))

sp_lookup <- read_csv("SPECIES.csv")
colnames(sp_lookup)
sp_lookup <- sp_lookup %>% select(SPECIES,NAME, Sp_Cat)
data <- left_join(data, sp_lookup, by = "SPECIES")

area <- read_csv("AREA new.csv")
colnames(area)
data <- left_join(data, area, by = "area")

# filter MHI
data <- filter(data, REGION == "MHI")

# convert na to 0 for lost to predator 
data <- data %>% mutate(number_lost_to_predator = ifelse(is.na(number_lost_to_predator), 0, number_lost_to_predator))


####################################################################################
# summarize fishing events deep 7 & Uku
# for estimation of pounds lost, use avg weight for day fished
# need to replace NA with avg fish weight when no fish caught, only lost
# but for now, remove 0 caught

Sum_by_day <- data %>% filter(number_caught > 0) %>%
  group_by(Year, Month, Island, license_number,day_fished, Sp_Cat) %>%
  summarize(NumCaught = sum(number_caught, na.rm=TRUE),
            LbsCaught = sum(pounds_caught, na.rm=TRUE),
            LostPred = sum(number_lost_to_predator, na.rm=TRUE),
            Dep = LostPred / (NumCaught + LostPred),
            NumHooked = NumCaught + LostPred,
            Weight = LbsCaught / NumCaught,
            Lbs_Lost = Weight * LostPred
  )

Sum_by_day <- Sum_by_day %>% 
  mutate(ZeroDep = ifelse(LostPred == 0,1,0))

# distribution of Depradation rates
Sum_by_day %>% filter(Sp_Cat == "Uku") %>% 
  ggplot( aes(x=Dep)) +
  geom_histogram( ) +
  ggtitle("Depredation Rates individual days fished - Uku") +
  theme_bw() +
  facet_wrap(~Year)

Sum_by_day %>% filter(Sp_Cat == "Deep 7") %>% 
  ggplot( aes(x=Dep)) +
  geom_histogram( ) +
  ggtitle("Depredation Rates individual days fished - Deep 7") +
  theme_bw() +
  facet_wrap(~Year)


# proportion of zeroes
Sum_by_day <- Sum_by_day %>% 
  mutate(ZeroDep = ifelse(LostPred == 0,1,0))


# summarize by license & year
Sum_by_license <- Sum_by_day %>% 
  group_by(Year, license_number, Sp_Cat) %>% 
  summarize(Days = n(),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            CPUE_lbs = LbsCaught / Days,
            CPUE_num = NumCaught / Days,
            CPUE_hooked = NumHooked / Days,
            CPUE_Lbs_hooked = (LbsCaught + Lbs_Lost) / Days
  )

# How many fishers never reported depredation
Sum_by_license <- Sum_by_license %>% 
  mutate(Perc_ZeroDep = ZeroDep / Days,
         NoDepred = ifelse(Perc_ZeroDep == 1,1,0))

Lic_no_depredation <- filter(Sum_by_license, Perc_ZeroDep == 1)

Lic_no_depredation %>% filter(Sp_Cat == "Uku") %>% 
  ggplot( aes(x=LbsCaught)) +
  geom_histogram( ) +
  ggtitle("Distribution of reported Lbs Caught Uku by fishers that never reported any depredation") +
  theme_bw() +
  facet_wrap(~Year)

Lic_no_depredation %>% filter(Sp_Cat == "Deep 7") %>% 
  ggplot( aes(x=LbsCaught)) +
  geom_histogram( ) +
  ggtitle("Distribution of reported Lbs Caught Deep 7 by fishers that never reported any depredation") +
  theme_bw() +
  facet_wrap(~Year)

# summarize by year
Sum_by_year <- Sum_by_license %>% 
  group_by(Year, Sp_Cat) %>% 
  summarize(Licenses = n(),
            Days = sum(Days),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep_mean, na.rm=TRUE),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            NoDepred = sum(NoDepred),
            CPUE_lbs = mean(CPUE_lbs, na.rm=TRUE),
            CPUE_num = mean(CPUE_num, na.rm=TRUE),
            CPUE_hooked = mean(CPUE_hooked, na.rm=TRUE),
            CPUE_Lbs_hooked = mean(CPUE_Lbs_hooked, na.rm=TRUE)
  )


# Investigate proportion of zero depredation
Sum_by_year <- Sum_by_year %>% 
  mutate(Perc_NoDepred = NoDepred / Licenses,
         ZeroDep_perc = ZeroDep / Days)

# plot
Sum_by_year  %>% 
  ggplot(aes(x=Year, y=ZeroDep_perc, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="Proportion of reports with no depredation") +
  theme_bw()

Sum_by_year  %>% 
  ggplot(aes(x=Year, y=Perc_NoDepred, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="Proportion of licenses with no depredation") +
  theme_bw()

# question - should we believe that fishers are accurately reporting depredation?
# or are high percentage of fishers not reporting depredation a sign of inaccurate reporting
# look into spatial differences in reporting

# plot mean depredation
Sum_by_year  %>% 
  ggplot(aes(x=Year, y=Dep_mean, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="Mean depredation") +
  theme_bw()

##########################################################
# sumarize by year & island

# summarize by license, island, & year
Sum_by_license_isl <- Sum_by_day %>% 
  group_by(Year, Island, license_number, Sp_Cat) %>% 
  summarize(Days = n(),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            CPUE_lbs = LbsCaught / Days,
            CPUE_num = NumCaught / Days,
            CPUE_hooked = NumHooked / Days,
            CPUE_Lbs_hooked = (LbsCaught + Lbs_Lost) / Days
  )

Sum_by_license_isl <- Sum_by_license_isl %>% 
  mutate(Perc_ZeroDep = ZeroDep / Days,
         NoDepred = ifelse(Perc_ZeroDep == 1,1,0))

Sum_by_isl_year <- Sum_by_license_isl %>% 
  group_by(Year, Island, Sp_Cat) %>% 
  summarize(Licenses = n(),
            Days = sum(Days),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep_mean, na.rm=TRUE),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            NoDepred = sum(NoDepred),
            CPUE_lbs = mean(CPUE_lbs, na.rm=TRUE),
            CPUE_num = mean(CPUE_num, na.rm=TRUE),
            CPUE_hooked = mean(CPUE_hooked, na.rm=TRUE),
            CPUE_Lbs_hooked = mean(CPUE_Lbs_hooked, na.rm=TRUE)
  )

Sum_by_isl_year <- Sum_by_isl_year %>% 
  mutate(Perc_NoDepred = NoDepred / Licenses,
         ZeroDep_perc = ZeroDep / Days)

# plot
Sum_by_isl_year  %>% 
  ggplot(aes(x=Year, y=ZeroDep_perc, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="Proportion of reports with no depredation") +
  facet_wrap(~ Island) +
  theme_bw()

Sum_by_isl_year  %>% 
  ggplot(aes(x=Year, y=Perc_NoDepred, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="Proportion of licenses with no depredation") +
  facet_wrap(~ Island) +
  theme_bw()

# look at mean depredation by year

# ignoring spatial differences, depredation mostly less than 5%
Sum_by_year  %>% 
  ggplot(aes(x=Year, y=Dep_mean, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="Mean depredation") +
  theme_bw()

# offshore/seamount trips are unique and sporadic, probably should be ommitted from assessment or atleast investigated more closely
# some really high reports of depredation

Sum_by_isl_year  %>% filter(Island != "Offshore") %>%
  ggplot(aes(x=Year, y=Dep_mean, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="Mean depredation") +
  facet_wrap(~ Island) +
  theme_bw()

# Uku depredation at Penguin bank above 15% fro 2021
####################################################################################
#  Explore seasonal patterns - focusing on Penguin Bank

Sum_by_license_isl_month <- Sum_by_day %>% 
  group_by(Year, Month, Island, license_number, Sp_Cat) %>% 
  summarize(Days = n(),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            CPUE_lbs = LbsCaught / Days,
            CPUE_num = NumCaught / Days,
            CPUE_hooked = NumHooked / Days,
            CPUE_Lbs_hooked = (LbsCaught + Lbs_Lost) / Days
  )

Sum_by_license_isl_month <- Sum_by_license_isl_month %>% 
  mutate(Perc_ZeroDep = ZeroDep / Days,
         NoDepred = ifelse(Perc_ZeroDep == 1,1,0))

Sum_by_isl_month_year <- Sum_by_license_isl_month %>% 
  group_by(Year, Month, Island, Sp_Cat) %>% 
  summarize(Licenses = n(),
            Days = sum(Days),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep_mean, na.rm=TRUE),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            NoDepred = sum(NoDepred),
            CPUE_lbs = mean(CPUE_lbs, na.rm=TRUE),
            CPUE_num = mean(CPUE_num, na.rm=TRUE),
            CPUE_hooked = mean(CPUE_hooked, na.rm=TRUE),
            CPUE_Lbs_hooked = mean(CPUE_Lbs_hooked, na.rm=TRUE)
  )

Sum_by_isl_month_year <- Sum_by_isl_month_year %>% 
  mutate(Perc_NoDepred = NoDepred / Licenses,
         ZeroDep_perc = ZeroDep / Days)

# investigate seasonality at Penguin Banks

PB_uku <- filter(Sum_by_isl_month_year, Sp_Cat == "Uku" & Island == "Penguin Bank" & Year > 2002)
PB_uku <- PB_uku %>% add_column(Time = seq(1,224,1))

ggplot(PB_uku, aes(x=Time, y=Dep_mean)) +
  geom_line() +
  geom_point() +
  labs(title="Mean Depredation - Uku - Penguin Bank") +
  scale_x_continuous(breaks = seq(1,224,12), labels=c(2003:2021))

# fill in missing months for Deep 7
# missing months due to closures? but need to check to confirm
yr_mnth <- tibble(Month   = c(1:12,1:7),
                  Year   = c(2003:2021))

yr_mnth <- yr_mnth %>% expand(Year, Month)

PB_deep7 <- filter(Sum_by_isl_month_year, Sp_Cat == "Deep 7" & Island == "Penguin Bank" & Year > 2002)
PB_deep7 <- left_join(yr_mnth, PB_deep7, by = c("Year" = "Year", "Month" = "Month"))
PB_deep7 <- PB_deep7 %>% add_column(Time = seq(1,228,1))

ggplot(PB_deep7, aes(x=Time, y=Dep_mean)) +
  geom_line() +
  geom_point() +
  labs(title="Mean Depredation - Deep 7 - Penguin Bank)") +
  scale_x_continuous(breaks = seq(1,224,12), labels=c(2003:2021))

PB_uku  %>%
  ggplot(aes(x=Year, y=Dep_mean)) +
  geom_line() +
  geom_point() +
  labs(title= "Mean depredation - Deep 7 - Penguin Bank") +
  facet_wrap(~ Month) +
  theme_bw()

PB_deep7  %>%
  ggplot(aes(x=Year, y=Dep_mean)) +
  geom_line() +
  geom_point() +
  labs(title= "Mean depredation - Deep 7 - Penguin Bank") +
  facet_wrap(~ Month) +
  theme_bw()

# Plot depredation monthly mean for timeseries

Sum_by_isl_monthly_avg <- Sum_by_isl_month_year %>% 
  group_by(Month, Island, Sp_Cat) %>% 
  summarize(Licenses = mean(Licenses),
            Days = mean(Days),
            NumCaught = mean(NumCaught),
            LbsCaught = mean(LbsCaught),
            LostPred = mean(LostPred),
            Lbs_Lost = mean(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep_mean, na.rm=TRUE),
            NumHooked = mean(NumHooked),
            ZeroDep = mean(ZeroDep),
            NoDepred = mean(NoDepred),
            CPUE_lbs = mean(CPUE_lbs, na.rm=TRUE),
            CPUE_num = mean(CPUE_num, na.rm=TRUE),
            CPUE_hooked = mean(CPUE_hooked, na.rm=TRUE),
            CPUE_Lbs_hooked = mean(CPUE_Lbs_hooked, na.rm=TRUE)
  )

Sum_by_isl_monthly_avg  %>% filter(Island != "Offshore") %>%
  ggplot(aes(x=Month, y=Dep_mean, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="Long-term Monthly Mean depredation") +
  facet_wrap(~ Island) +
  scale_x_continuous(breaks=c(1:12),
                   labels=c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")) +
  theme_bw()

####################################################################################
# summarize fishing events deep 7 by species

deep7 <- filter(data, Sp_Cat == "Deep 7")

deep7_Sum_by_day <- deep7 %>%
  group_by(Year, Month, Island, license_number,day_fished, NAME) %>%
  summarize(NumCaught = sum(number_caught, na.rm=TRUE),
            LbsCaught = sum(pounds_caught, na.rm=TRUE),
            LostPred = sum(number_lost_to_predator, na.rm=TRUE),
            Dep = LostPred / (NumCaught + LostPred),
            NumHooked = NumCaught + LostPred,
            Weight = LbsCaught / NumCaught,
            Lbs_Lost = Weight * LostPred
  )

# proportion of zeroes
deep7_Sum_by_day <- deep7_Sum_by_day %>% 
  mutate(ZeroDep = ifelse(LostPred == 0,1,0))


# summarize by license & year
deep7_Sum_by_license <- deep7_Sum_by_day %>% 
  group_by(Year, license_number, NAME) %>% 
  summarize(Days = n(),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            CPUE_lbs = LbsCaught / Days,
            CPUE_num = NumCaught / Days,
            CPUE_hooked = NumHooked / Days,
            CPUE_Lbs_hooked = (LbsCaught + Lbs_Lost) / Days
  )

# How many fishers never reported depredation
deep7_Sum_by_license <- deep7_Sum_by_license %>% 
  mutate(Perc_ZeroDep = ZeroDep / Days,
         NoDepred = ifelse(Perc_ZeroDep == 1,1,0))

deep7_Lic_no_depredation <- filter(deep7_Sum_by_license, Perc_ZeroDep == 1)

deep7_Lic_no_depredation %>% filter(NAME == "Opakapaka") %>% 
  ggplot( aes(x=LbsCaught)) +
  geom_histogram( ) +
  ggtitle("Distribution of reported Lbs Caught Uku by fishers that never reported any depredation") +
  theme_bw() +
  facet_wrap(~Year)

deep7_Lic_no_depredation %>% filter(NAME == "Onaga") %>% 
  ggplot( aes(x=LbsCaught)) +
  geom_histogram( ) +
  ggtitle("Distribution of reported Lbs Caught Deep 7 by fishers that never reported any depredation") +
  theme_bw() +
  facet_wrap(~Year)

# summarize by year
deep7_Sum_by_year <- deep7_Sum_by_license %>% 
  group_by(Year, NAME) %>% 
  summarize(Licenses = n(),
            Days = sum(Days),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep_mean, na.rm=TRUE),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            NoDepred = sum(NoDepred),
            CPUE_lbs = mean(CPUE_lbs, na.rm=TRUE),
            CPUE_num = mean(CPUE_num, na.rm=TRUE),
            CPUE_hooked = mean(CPUE_hooked, na.rm=TRUE),
            CPUE_Lbs_hooked = mean(CPUE_Lbs_hooked, na.rm=TRUE)
  )

# Investigate proportion of zero depredation
deep7_Sum_by_year <- deep7_Sum_by_year %>% 
  mutate(Perc_NoDepred = NoDepred / Licenses,
         ZeroDep_perc = ZeroDep / Days)

# plot
deep7_Sum_by_year  %>% 
  ggplot(aes(x=Year, y=ZeroDep_perc)) +
  geom_line() +
  geom_point() +
  labs(title="Proportion of reports with no depredation") +
  theme_bw() + 
  facet_wrap(~ NAME)

deep7_Sum_by_year  %>% 
  ggplot(aes(x=Year, y=Perc_NoDepred)) +
  geom_line() +
  geom_point() +
  labs(title="Proportion of licenses with no depredation") +
  theme_bw() + 
  facet_wrap(~ NAME) 

deep7_Sum_by_year  %>% 
  ggplot(aes(x=Year, y=Dep_mean)) +
  geom_line() +
  geom_point() +
  labs(title="Mean depredation") +
  facet_wrap(~ NAME) +
  theme_bw() 

# explore seasonality by species

#####################################
# sumarize by deep 7 species & island

# summarize by license & year
deep7_Sum_by_isl_license <- deep7_Sum_by_day %>% 
  group_by(Year, Island, license_number, NAME) %>% 
  summarize(Days = n(),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            CPUE_lbs = LbsCaught / Days,
            CPUE_num = NumCaught / Days,
            CPUE_hooked = NumHooked / Days,
            CPUE_Lbs_hooked = (LbsCaught + Lbs_Lost) / Days
  )

# How many fishers never reported depredation
deep7_Sum_by_isl_license <- deep7_Sum_by_isl_license %>% 
  mutate(Perc_ZeroDep = ZeroDep / Days,
         NoDepred = ifelse(Perc_ZeroDep == 1,1,0))

deep7_Isl_Lic_no_depredation <- filter(deep7_Sum_by_isl_license, Perc_ZeroDep == 1)

deep7_Isl_Lic_no_depredation %>% filter(NAME == "Opakapaka" & Island == "Penguin Bank") %>% 
  ggplot( aes(x=LbsCaught)) +
  geom_histogram( ) +
  ggtitle("Distribution of reported Lbs Caught Uku by fishers that never reported any depredation") +
  theme_bw() +
  facet_wrap(~Year)

deep7_Isl_Lic_no_depredation %>% filter(NAME == "Onaga" & Island == "Penguin Bank") %>% 
  ggplot( aes(x=LbsCaught)) +
  geom_histogram( ) +
  ggtitle("Distribution of reported Lbs Caught Deep 7 by fishers that never reported any depredation") +
  theme_bw() +
  facet_wrap(~Year)

# summarize by year
deep7_Sum_by_isl_year <- deep7_Sum_by_isl_license %>% 
  group_by(Year, Island, NAME) %>% 
  summarize(Licenses = n(),
            Days = sum(Days),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep_mean, na.rm=TRUE),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            NoDepred = sum(NoDepred),
            CPUE_lbs = mean(CPUE_lbs, na.rm=TRUE),
            CPUE_num = mean(CPUE_num, na.rm=TRUE),
            CPUE_hooked = mean(CPUE_hooked, na.rm=TRUE),
            CPUE_Lbs_hooked = mean(CPUE_Lbs_hooked, na.rm=TRUE)
  )

# Investigate proportion of zero depredation
deep7_Sum_by_isl_year <- deep7_Sum_by_isl_year %>% 
  mutate(Perc_NoDepred = NoDepred / Licenses,
         ZeroDep_perc = ZeroDep / Days)

# plot
deep7_Sum_by_isl_year  %>% 
  ggplot(aes(x=Year, y=ZeroDep_perc)) +
  geom_line() +
  geom_point() +
  labs(title="Proportion of reports with no depredation") +
  theme_bw() + 
  facet_grid(Island ~ NAME)

# without offshore areas
deep7_Sum_by_isl_year  %>% filter(Island != "Offshore") %>%
  ggplot(aes(x=Year, y=ZeroDep_perc)) +
  geom_line() +
  geom_point() +
  labs(title="Proportion of reports with no depredation") +
  theme_bw() + 
  facet_grid(Island ~ NAME)

deep7_Sum_by_isl_year  %>% 
  ggplot(aes(x=Year, y=Perc_NoDepred)) +
  geom_line() +
  geom_point() +
  labs(title="Proportion of licenses with no depredation") +
  theme_bw() + 
  facet_grid(Island ~ NAME) 

# explore spatial differences in depredation by species
# note that fishers are highly unlikely to be able to 100% predict what species was lost to a predator unless they retrieve a portion of the fish
# most data is based on assumption of target catch and catch composition of the day
# expereinced fishers can tell by feeling the bite as different species bite differently

deep7_Sum_by_isl_year  %>% filter(Island != "Offshore") %>%
  ggplot(aes(x=Year, y=Dep_mean)) +
  geom_line() +
  geom_point() +
  labs(title="Mean depredation") +
  facet_grid(Island ~ NAME) +
  theme_bw() 

########################
# what does this mean at the "assessment" level?
# look at change in CPUE (incorporating lost fish to predation)

# first look at cpue based on numbers

Sum_by_year_long_num <- Sum_by_year %>%
  select(Year, Sp_Cat, CPUE_num,CPUE_hooked) %>%
  gather(CPUE_meth,CPUE,CPUE_num:CPUE_hooked) 

ggplot(Sum_by_year_long_num, aes(x=Year, y=CPUE, group=CPUE_meth)) +
  geom_line(aes(color=CPUE_meth)) +
  geom_point(aes(color=CPUE_meth)) +
  facet_wrap(~ Sp_Cat)

# look at cpue based on lbs

Sum_by_year_long_lbs <- Sum_by_year %>%
  select(Year, Sp_Cat, CPUE_lbs,CPUE_Lbs_hooked) %>%
  gather(CPUE_meth,CPUE,CPUE_lbs:CPUE_Lbs_hooked) 
  
ggplot(Sum_by_year_long_lbs, aes(x=Year, y=CPUE, group=CPUE_meth)) +
  geom_line(aes(color=CPUE_meth)) +
  geom_point(aes(color=CPUE_meth)) +
  facet_wrap(~ Sp_Cat)
  
# plot total landings and pounds lost
# note, the big question here is whether reported depredation is believable 
# or whether there is reason to beleive it is not reported accirately and 
# therefore requires some kind of extrapolatin based on reported losses

Sum_by_year_long_lbs_lost <- Sum_by_year %>%
  select(Year, Sp_Cat, LbsCaught,Lbs_Lost) %>%
  mutate(LbsCaught_Lost = LbsCaught + Lbs_Lost) %>%
  select(-Lbs_Lost) %>%
  gather(Landings_class,Lbs,LbsCaught:LbsCaught_Lost) 

ggplot(Sum_by_year_long_lbs_lost, aes(x=Year, y=Lbs, group=Landings_class)) +
  geom_line(aes(color=Landings_class)) +
  geom_point(aes(color=Landings_class)) +
  facet_wrap(~ Sp_Cat)

###############################################################
# investigate method differences in depredation

Sum_by_day_method <- data %>% filter(number_caught > 0) %>%
  group_by(Year, Month, Island, license_number,day_fished, method_name, Sp_Cat) %>%
  summarize(NumCaught = sum(number_caught, na.rm=TRUE),
            LbsCaught = sum(pounds_caught, na.rm=TRUE),
            LostPred = sum(number_lost_to_predator, na.rm=TRUE),
            Dep = LostPred / (NumCaught + LostPred),
            NumHooked = NumCaught + LostPred,
            Weight = LbsCaught / NumCaught,
            Lbs_Lost = Weight * LostPred
  )

Sum_by_day_method <- Sum_by_day_method %>% 
  mutate(ZeroDep = ifelse(LostPred == 0,1,0))

# proportion of zeroes
Sum_by_day_method <- Sum_by_day_method %>% 
  mutate(ZeroDep = ifelse(LostPred == 0,1,0))


# summarize by license & year
Sum_by_license_method <- Sum_by_day_method %>% 
  group_by(Year, license_number, method_name, Sp_Cat) %>% 
  summarize(Days = n(),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            CPUE_lbs = LbsCaught / Days,
            CPUE_num = NumCaught / Days,
            CPUE_hooked = NumHooked / Days,
            CPUE_Lbs_hooked = (LbsCaught + Lbs_Lost) / Days
  )

# How many fishers never reported depredation
Sum_by_license_method <- Sum_by_license_method %>% 
  mutate(Perc_ZeroDep = ZeroDep / Days,
         NoDepred = ifelse(Perc_ZeroDep == 1,1,0))

Lic_no_depredation <- filter(Sum_by_license_method, Perc_ZeroDep == 1)

# summarize by year
Sum_by_year_method <- Sum_by_license_method %>% 
  group_by(Year, method_name, Sp_Cat) %>% 
  summarize(Licenses = n(),
            Days = sum(Days),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Lbs_Lost = sum(Lbs_Lost, na.rm=TRUE),
            Dep_mean = mean(Dep_mean, na.rm=TRUE),
            NumHooked = sum(NumHooked),
            ZeroDep = sum(ZeroDep),
            NoDepred = sum(NoDepred),
            CPUE_lbs = mean(CPUE_lbs, na.rm=TRUE),
            CPUE_num = mean(CPUE_num, na.rm=TRUE),
            CPUE_hooked = mean(CPUE_hooked, na.rm=TRUE),
            CPUE_Lbs_hooked = mean(CPUE_Lbs_hooked, na.rm=TRUE)
  )

# Investigate proportion of zero depredation
Sum_by_year_method <- Sum_by_year_method %>% 
  mutate(Perc_NoDepred = NoDepred / Licenses,
         ZeroDep_perc = ZeroDep / Days)

# plot depredation by method
Sum_by_year_method  %>% filter(Sp_Cat == "Uku") %>%
  ggplot(aes(x=Year, y=Dep_mean)) +
  geom_line() +
  geom_point() +
  labs(title="Mean depredation - Uku") +
  theme_bw() +
  facet_wrap(~ method_name)

# many of the methods are data errors, refine list
Sum_by_year_method  %>% filter(Sp_Cat == "Uku" & method_name %in% c("Deep-Sea Handline, Bottom Handline","Trolling - Bait","Trolling - Lures", "Inshore Handline, Cowrie Shell (Tako)","Casting, Rod & Reel, Jigging")) %>%
  ggplot(aes(x=Year, y=Dep_mean)) +
  geom_line() +
  geom_point() +
  labs(title="Mean depredation - Uku") +
  theme_bw() +
  facet_wrap(~ method_name)

Sum_by_year_method  %>% filter(Sp_Cat == "Deep 7") %>%
  ggplot(aes(x=Year, y=Dep_mean)) +
  geom_line() +
  geom_point() +
  labs(title="Mean depredation - Deep 7") +
  theme_bw() +
  facet_wrap(~ method_name)

# many of the methods are data errors, refine list
Sum_by_year_method  %>% filter(Sp_Cat == "Deep 7" & method_name %in% c("Deep-Sea Handline, Bottom Handline","Palu 'Ahi")) %>%
  ggplot(aes(x=Year, y=Dep_mean)) +
  geom_line() +
  geom_point() +
  labs(title="Mean depredation - Deep 7") +
  theme_bw() +
  facet_wrap(~ method_name)

# supposedly method gear change where some fishers are using up to 6,7 lines
# explore whether this is having any effect in Uku fishery

uku <- data %>% filter(Sp_Cat == "Uku") %>% mutate(Gear_cat = ifelse(number_of_gear_units >6, "Gr 7", as.character(number_of_gear_units)))

uku_num_na <- filter(uku, is.na(number_caught))

Uku_Sum_by_gear_day <- uku %>%
  group_by(Year, license_number,day_fished, method_name, Gear_cat) %>%
  summarize(NumCaught = sum(number_caught, na.rm=TRUE),
            LbsCaught = sum(pounds_caught, na.rm=TRUE),
            LostPred = sum(number_lost_to_predator, na.rm=TRUE),
            Dep = LostPred / (NumCaught + LostPred),
            NumHooked = NumCaught + LostPred
  )

Uku_Sum_by_gear_lic <- Uku_Sum_by_gear_day %>%
  group_by(Year, license_number, method_name, Gear_cat) %>%
  summarize(Days = n(),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Dep_mean = mean(Dep, na.rm=TRUE),
            NumHooked = sum(NumHooked),
            CPUE_lbs = LbsCaught / Days,
            CPUE_num = NumCaught / Days,
            CPUE_hooked = NumHooked / Days
  )

Uku_Sum_by_gear <- Uku_Sum_by_gear_lic %>%
  group_by(Year, method_name,Gear_cat) %>%
  summarize(Licenses = n(),
            Days = sum(Days),
            NumCaught = sum(NumCaught),
            LbsCaught = sum(LbsCaught),
            LostPred = sum(LostPred),
            Dep_mean = mean(Dep_mean, na.rm=TRUE),
            NumHooked = sum(NumHooked),
            CPUE_lbs = mean(CPUE_lbs, na.rm=TRUE),
            CPUE_num = mean(CPUE_num, na.rm=TRUE),
            CPUE_hooked = mean(CPUE_hooked, na.rm=TRUE)
  )

# plot depredation by method and number of gear units
Uku_Sum_by_gear %>% filter(method_name %in% c("Deep-Sea Handline, Bottom Handline","Trolling - Bait","Trolling - Lures", "Inshore Handline, Cowrie Shell (Tako)","Casting, Rod & Reel, Jigging")) %>%
  ggplot(aes(x=Year, y=Dep_mean, color=Gear_cat)) +
  geom_line() +
  geom_point() +
  labs(title="Mean depredation - Uku") +
  facet_wrap(~ method_name, ncol=2)

########################################################################################################
# Explore CPUE

Sum_by_year  %>% 
  ggplot(aes(x=Year, y=CPUE_lbs, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="CPUE (lbs / day") +
  theme_bw() 

Sum_by_year  %>% 
  ggplot(aes(x=Year, y=CPUE_num, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="CPUE (# / day)") +
  theme_bw() 

Sum_by_isl_year  %>% 
  ggplot(aes(x=Year, y=CPUE_lbs, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="CPUE (lbs / day") +
  facet_wrap(~ Island) +
  theme_bw() 

Sum_by_isl_year  %>% filter(Island != "Offshore") %>%
  ggplot(aes(x=Year, y=CPUE_lbs, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="CPUE (lbs / day") +
  facet_wrap(~ Island) +
  theme_bw() 

Sum_by_isl_year  %>% filter(Island != "Offshore") %>%
  ggplot(aes(x=Year, y=CPUE_num, col = Sp_Cat)) +
  geom_line() +
  geom_point() +
  labs(title="CPUE (# / day)") +
  facet_wrap(~ Island) +
  theme_bw() 

# NOTE this way of calculating CPUE is how WPacFIN does it for the SAFE report
# it only factors trips when the species was reported as caught, this is totally misleading
# one should not expect to catch 20 lbs of hapupu per day! It is a rare catch
# this basicly says if you do catch a hapupupu, you catch ~ 20 lbs
# need to first classify as a deep 7 trip based on some type of criteria (catch composition)

deep7_Sum_by_year  %>% 
  ggplot(aes(x=Year, y=CPUE_lbs)) +
  geom_line() +
  geom_point() +
  labs(title="CPUE (lbs / day)") +
  facet_wrap(~ NAME) +
  theme_bw() 

deep7_Sum_by_year  %>% 
  ggplot(aes(x=Year, y=CPUE_num)) +
  geom_line() +
  geom_point() +
  labs(title="CPUE (# / day)") +
  facet_wrap(~ NAME) +
  theme_bw()

Uku_Sum_by_gear %>% filter(method_name %in% c("Deep-Sea Handline, Bottom Handline","Trolling - Bait", "Trolling - Lures", "Inshore Handline, Cowrie Shell (Tako)")) %>%
  ggplot(aes(x=Year, y=CPUE_lbs, color=Gear_cat)) +
  geom_line() +
  geom_point() +
  labs(title="Mean CPUE - Uku") +
  facet_wrap(~ method_name, ncol=2)

# investigate high CPUE for trolling-bait
c <- Uku_Sum_by_gear_lic %>% filter(method_name %in% c("Trolling - Bait") & CPUE_lbs > 400)

