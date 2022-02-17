library(tidyverse)
library(r4ss)
library(magrittr)

camera_lengths <- read.csv("./Data/BFISH 2016-2019 Camera Lengths.csv")
fishing_lengths <- read.csv("./Data/BFISH 2016-2019 Research Fishing Lengths.csv")

head(camera_lengths)
head(fishing_lengths)

sp_code <- "PRFI"

opaka_lengths <- camera_lengths %>% 
filter(str_detect(SPECIES_CD, sp_code)) %>% 
separate(BFISH, into = c("BFISH", "YEAR", "SEASON"), sep = "_") %>% 
filter(str_detect(SEASON, "F")) %>% 
select(YEAR, LENGTH_CM) %>% 
mutate(GEAR = "Camera")

opaka_lengths <- fishing_lengths %>% 
filter(str_detect(SPECIES_CD, sp_code)) %>% 
separate(BFISH, into = c("BFISH", "YEAR", "SEASON"), sep = "_") %>% 
filter(str_detect(SEASON, "F")) %>% 
select(YEAR, LENGTH_CM) %>% 
filter(!is.na(LENGTH_CM)) %>% 
mutate(GEAR = "Fishing") %>% 
bind_rows(opaka_lengths) %>% 
mutate(#bin = cut(LENGTH_CM, seq(min(LENGTH_CM), max(LENGTH_CM) + 2, 2), right = FALSE),
SS_YEAR = as.numeric(YEAR) + 1)


opaka_lengths %>% 
ggplot(aes(LENGTH_CM, colour = GEAR)) +
geom_freqpoly(aes(y = stat(count / sum(count)))) +
facet_wrap(~YEAR, scales = "free") +
theme_classic()


plot_fun = function(data, x, z, bw) {
     ggplot(data, aes(x = .data[[x]], colour = .data[[z]]) ) +
          geom_freqpoly(binwidth = bw) +
          facet_wrap(~YEAR, scales = "free")+
          theme_bw()
}

bw.vec <- seq(1, 5, by = 1)

opaka_lengths %>% 
plot_fun(., x = "LENGTH_CM", z = "GEAR", bw = bw.vec)

plots = map(bw.vec, ~ggplot(data = opaka_lengths, aes(x = LENGTH_CM, colour = GEAR) ) +
          geom_freqpoly(aes(y = stat(count / sum(count))), binwidth = .x)  +
          facet_wrap(~YEAR, scales = "free")+
          theme_bw()
)

walk2(.x = plots, .y = bw.vec, 
  ~ ggsave(filename = paste0("length_bin", .y, ".png"), plot = .x))

head(opaka_lengths)
range(opaka_lengths$LENGTH_CM)
opaka_lengths %>% 
group_by(GEAR, SS_YEAR) %>% 
#ggplot(aes(x = YEAR, y = LENGTH_CM, fill = GEAR)) +
#geom_boxplot()
summarise(min = min(LENGTH_CM),
max = max(LENGTH_CM),
n = n(), 
mean = mean(LENGTH_CM), 
median = median(LENGTH_CM))

bins <- sort(unique(opaka_lengths$bin))

opaka.dat <- SS_readdat_3.30(file = "./Run3/BFISH_Length_Comps/fish_cam/opakadat.ss")

# Adding in Length Comps
##Population_length_bins
### Change it from generating bins, to using databins
opaka.dat$lbin_method <- 1
opaka.dat$binwidth <- NULL
opaka.dat$minimum_size <- NULL
opaka.dat$maximum_size <- NULL
opaka.dat$use_lencomp <- 1
opaka.dat$N_lencomp <- 4
opaka.dat$lbin_vector_pop <- seq(1, 90, by = 2)

head(opaka_lengths)
### Create length comp table
lc.2017 <- opaka_lengths %>% 
filter(SS_YEAR == 2017) %>% 
mutate(bin = cut(LENGTH_CM, seq(1, 90, by = 2), right = FALSE))
lc.2018 <- opaka_lengths %>% 
filter(SS_YEAR == 2018) %>% 
mutate(bin = cut(LENGTH_CM, seq(1, 90, by = 2), right = FALSE))
lc.2019 <- opaka_lengths %>% 
filter(SS_YEAR == 2019) %>% 
mutate(bin = cut(LENGTH_CM, seq(1, 90, by = 2), right = FALSE))
lc.2020 <- opaka_lengths %>% 
filter(SS_YEAR == 2020) %>% 
mutate(bin = cut(LENGTH_CM, seq(1, 90, by = 2), right = FALSE))

length_comps <- rbind(summary(lc.2017$bin), summary(lc.2018$bin), summary(lc.2019$bin), summary(lc.2020$bin))

str_remove_all(colnames(length_comps), ", ") 

lc.cols <- data.frame(Method = rep(2, 4), Yr = sort(unique(opaka_lengths$SS_YEAR)), Seas = rep(1, 4), FltSvy = rep(4, 4), Gender = rep(0,4), Part = rep(2, 4), Nsamp = c(136,264,126,301))

bind_cols(lc.cols, length_comps)

#if FltSvy for commercial fleet is positive, make it negative so that it is not included in the likelihood component right now
size_freq <- opaka.dat$sizefreq_data_list[[1]] 
size_freq %<>% 
mutate(FltSvy = as.numeric(FltSvy)*-1)
head(size_freq)

#opaka.dat$N_sizefreq_methods <- 2
#opaka.dat$nbins_per_method <- c(30, 90)
#opaka.dat$units_per_method <- c(2, 2) #numbers
#opaka.dat$scale_per_method <- c(1, 3) #kg and cm
#opaka.dat$mincomp_per_method <- c(1e-09, 1e-09)
#opaka.dat$Nobs_per_method <- c(69, 4)

opaka.dat$sizefreq_bins_list <- list(c(0.23, 0.68, 1.13, 1.59, 2.04, 2.49, 2.95, 3.40, 3.86, 4.31, 4.76, 5.22, 5.67, 6.12, 6.58, 7.03, 7.48, 7.94, 8.39, 8.85, 9.30, 9.75, 10.21, 10.66, 11.11, 11.57, 12.02, 12.47, 12.93, 13.38), c(seq(1, 90, by = 1)))

length.comps <- data.frame(Method = rep(2, 4), Yr = unique(opaka_lengths$SS_YEAR), Seas = rep(1, 4), FltSvy = rep(4, 4), Gender = rep(0,4), Part = rep(2, 4), Nsamp = c(136,264,126,301))


SS_writedat_3.30(opaka.dat, outfile = "./Run3/BFISH_Length_Comps/fish_cam/opakadat.ss", overwrite = TRUE)

### Weight and length relationship
a <- 6.378e-5
b <- 2.7621
w <- unlist(opaka.dat$sizefreq_bins_list)
wlbs <- w*2.205
l <- opaka.dat$lbin_vector_pop

format((a*l^b)/2.205, scientific = FALSE, digits = 1)
