t1 = Sys.time()

library(sas7bdat)
library(zoo)
library(KFAS)

setwd("/Users/Toby/Documents/GitHub/MHI_Bottomfish_2023/R/HMRFS")
i1 = read.sas7bdat("./MRIP_Ifiles_HI/i1.sas7bdat", debug = TRUE)
observed_catch = read.sas7bdat("./MRIP_Ifiles_HI/i3.sas7bdat", debug = TRUE)
unavailable_catch = read.sas7bdat("./MRIP_Ifiles_HI/i2.sas7bdat", debug = TRUE)

species_df = data.frame("key" = c(8835020413, 8835360304, 8835360302, 8835360704, 8835360706, 8835360901, 8835360707),
                        "common_name" = c("HAWAIIAN GROUPER", "LONGTAILED RED SNAPPER", "RUBY SNAPPER", "PINK SNAPPER", "VON SIEBOLDS SNAPPER", "IRONJAW SNAPPER", "BINGHAMS SNAPPER"),
                        "hawaiian_name" = c("Hapu'upu'u", "Onaga", "Ehu", "Opakapaka", "Kalekale", "Lehi", "Gindai"),
                        "alpha" = c(2.884, 2.673, 3.026, 2.928, 2.932, 2.458, 2.859), # cm
                        "beta" = c(3.065*10^-5, 6.005*10^-5, 1.551*10^-5, 2.311*10^-5, 2.243*10^-5, 1.298*10^-4, 3.526*10^-5))

# for mode: 1 = private boat, 2 = shore
# for area: 1 = ocean (> 3 mi), 2 = ocean (<= 3 mi), 3 = inland
# for disposition: index 1 = sold, 2 = not sold
i1$mode = ifelse(i1$MODE_F == 8, 1, ifelse(i1$MODE_FX %in% 1:5, 2, NA))
i1$area_exp = ifelse(i1$AREA_X == 2, 1, ifelse(i1$AREA_X == 1, 2, ifelse(i1$AREA_X == 5, 3, NA)))
observed_catch = observed_catch[observed_catch$DISP3 %in% c(3, 4, 5) & observed_catch$MODE_F != 7,] # fish that were not released and not from a charter trip
observed_catch$sold = ifelse(observed_catch$DISP3 == 5, T, F)
unavailable_catch = unavailable_catch[unavailable_catch$DISPO %in% c(3, 4, 5) & unavailable_catch$MODE_F != 7,] # fish that were not released and not from a charter trip
unavailable_catch$sold = ifelse(unavailable_catch$DISPO == 5, T, F)

observed_catch$EST_WGT = sapply(1:nrow(observed_catch), function(f) {
  r = observed_catch[f,]
  sp = r$SP_CODE
  if((sp %in% species_df$key) && !is.nan(r$LNGTH) && is.nan(r$WGT)) { # bottomfish with length, but no weight
    i = which(species_df$key == r$SP_CODE)
    return(species_df[i,]$beta * (r$LNGTH / 10) ^ species_df[i,]$alpha)
  } else {
    return(NaN)
  }
})

years = sort(unique(c(observed_catch$YEAR, unavailable_catch$year)), decreasing = F)
trips = unique(i1[!is.na(i1$mode) & !is.na(i1$ID_CODE),]$ID_CODE)

n_years = length(years)
n_waves = 6
n_species = nrow(species_df)
n_modes = 2
n_areas = 3
n_trips = length(trips)
n_dispositions = 2

num_trips = array(0, dim = c(n_years, n_waves, n_modes, n_areas))

for(y in 1:n_years) {
  for(w in 1:n_waves) {
    for(m in 1:n_modes) {
      for(a in 1:n_areas) {
        num_trips[y, w, m, a] = length(unique((i1[i1$YEAR == years[y] & i1$WAVE == w & i1$mode == m & i1$area_exp == a  & !is.na(i1$ID_CODE),])$ID_CODE))
      }
    }
  }
}

anglers_by_trip = array(0, dim = c(n_years, n_waves, n_modes, n_areas, n_trips))
observed_caught_by_trip = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas, n_trips))
observed_total_caught = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas, n_dispositions))
num_weighed = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))
num_estimated = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))
total_weight = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))
unavailable_total_caught = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas, n_dispositions))

for(t in 1:n_trips) {
  trip = trips[t]
  year = as.numeric(substr(trip, 6, 9))
  y = which(years == year)
  w = ceiling(as.numeric(substr(trip, 10, 11)) / 2)
  m = i1[i1$ID_CODE == trip,]$mode
  a = i1[i1$ID_CODE == trip,]$area_exp
  
  anglers_by_trip[y, w, m, a, t] = max(c(i1[i1$ID_CODE == trip,]$CNTRBTRS, observed_catch[observed_catch$ID_CODE == trip,]$CNTRBTRS, unavailable_catch[unavailable_catch$ID_CODE == trip,]$CNTRBTRS))
  
  # observed catch
  # FSHINSP is the number of specimens of a species observed
  #   Each row represents one specimen measured, but all have the same FSHINSP
  # catch can be from a group of fishers, and effort should be in angler trips
  #  CNTRBTRS is the number of fishers the catch represents
  
  observed = observed_catch[observed_catch$ID_CODE == trip,]
  
  if(nrow(observed) > 0) {
    species = unique(observed$SP_CODE)
    species = species[species %in% species_df$key]
    
    if(length(species) > 0) {
      for(i in 1:length(species)) {
        s = which(species_df$key == species[i])
        
        observed_s = observed[observed$SP_CODE == species[i],]
        
        observed_caught_by_trip[y, w, s, m, a, t] = max(observed_s$FSHINSP)
        
        if(all(observed_s$sold)) {
          observed_total_caught[y, w, s, m, a, 1] = observed_total_caught[y, w, s, m, a, 1] + max(observed_s$FSHINSP)
        } else if(any(!observed$sold)) {
          observed_total_caught[y, w, s, m, a, 2] = observed_total_caught[y, w, s, m, a, 2] + max(observed_s$FSHINSP)
        } 
#        else {
          # there are 7 instances where there are multiple dispositions within the same interview/species
          # in each case there is a row for each fish caught, but this could potentially not be the case
          # e.g. ID_CODE = 1700920031105001, SP_CODE = 8835360704
#          count_sold = 0
#          count_unsold = 0
          
#          for(j in 1:nrow(observed_s)) {
#            r = observed_s[j,]
#            if(r$sold) {
#              count_sold = count_sold + 1
#            } else {
#              count_unsold = count_unsold + 1
#            }
            
#            observed_total_caught[y, w, s, m, a, 1] = observed_total_caught[y, w, s, m, a, 1] + max(observed_s$FSHINSP) * count_sold / (count_sold + count_unsold)
#            observed_total_caught[y, w, s, m, a, 2] = observed_total_caught[y, w, s, m, a, 2] + max(observed_s$FSHINSP) * count_unsold / (count_sold + count_unsold)
#          }
#        }
        
        for(j in 1:nrow(observed_s)) {
          r = observed_s[j,]
          
          if(!is.nan(r$WGT)) {
            num_weighed[y, w, s, m, a] = num_weighed[y, w, s, m, a] + 1
            total_weight[y, w, s, m, a] = total_weight[y, w, s, m, a] + r$WGT
          } else if(!is.nan(r$LNGTH)) {
            num_estimated[y, w, s, m, a] = num_estimated[y, w, s, m, a] + 1
            total_weight[y, w, s, m, a] = total_weight[y, w, s, m, a] + r$EST_WGT
          }
        }
      }
    }
  }
  
  # unavailable catch
  # can only represent individual fishers
  # NUM2: sequential numbering of species within interview
  # NUM_TYP2: total number of species within interview
  # NUM_FISH: number of of fish of each species
  
  unavailable = unavailable_catch[unavailable_catch$ID_CODE == trip,]
  
  if(nrow(unavailable) > 0) {
    species = unique(unavailable$SP_CODE)
    species = species[species %in% species_df$key]
    
    if(length(species) > 0) {
      for(i in 1:length(species)) {
        s = which(species_df$key == species[i])
        
        unavailable_s = unavailable[unavailable$SP_CODE == species[i],]
        
        #if(all(unavailable_s$sold)) {
        #  unavailable_total_caught[y, w, s, m, a, 1] = unavailable_total_caught[y, w, s, m, a, 1] + max(unavailable_s$NUM_FISH)
        #} else if(any(unavailable_s$sold)) {
        #  unavailable_total_caught[y, w, s, m, a, 2] = unavailable_total_caught[y, w, s, m, a, 2] + max(unavailable_s$NUM_FISH)
        #}
        # There is an instance of two rows for the same specie within a single trip (ID_CODE = 1701620190402001, SP_CODE = 8835360704). The below code counts both rows, each with 7 fish.
        for(j in 1:nrow(unavailable_s)) {
          r = unavailable_s[j,]
          
          if(r$sold) {
            unavailable_total_caught[y, w, s, m, a, 1] = unavailable_total_caught[y, w, s, m, a, 1] + r$NUM_FISH
          } else {
            unavailable_total_caught[y, w, s, m, a, 2] = unavailable_total_caught[y, w, s, m, a, 2] + r$NUM_FISH
          }
        }
      }
    }
  }
}

mean_weight_wave = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))
mean_weight_year = array(0, dim = c(n_years, n_species, n_modes, n_areas))
mean_weight_all = array(0, dim = c(n_species))
catch_rate = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))

for(s in 1:n_species) {
  mean_weight_wave[, , s, , ] = total_weight[, , s, , ] / (num_weighed[, , s, , ] + num_estimated[, , s, , ])
  mean_weight_year[, s, , ] = apply(total_weight[, , s, , ], c(1, 3, 4), sum) / apply(num_weighed[, , s, , ] + num_estimated[, , s, , ], c(1, 3, 4), sum)
  mean_weight_all[s] = sum(total_weight[, , s, , ]) / (sum(num_weighed[, , s, , ] + num_estimated[, , s, , ]))
  catch_rate[, , s, , ] = apply(observed_total_caught[, , s, , , ] + unavailable_total_caught[, , s, , , ], c(1, 2, 3, 4), sum) / apply(anglers_by_trip, c(1, 2, 3, 4), sum)
}
catch_rate[18, 3, , , ] = catch_rate[17, 3, , , ] # No data available for April - June 2020, so use previous year
catch_rate[is.nan(catch_rate)] = 0

catch_rate_var = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))
unavailable_catch_rate_var = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))
observed_catch_rate_var = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))

for(y in 1:n_years) {
  for(w in 1:n_waves) {
    for(s in 1:n_species) {
      for(m in 1:n_modes) {
        for(a in 1:n_areas) {
          unavailable_catch_rate = sum(unavailable_total_caught[y, w, s, m, a, ]) / num_trips[y, w, m, a]
          unavailable_catch_rate_var[y, w, s, m, a] = 1 / num_trips[y, w, m, a] * sum((unavailable_total_caught[y, w, s, m, a, ] - unavailable_catch_rate) ^ 2 / (num_trips[y, w, m, a] - 1))
          
          observed_catch_rate = (sum(observed_total_caught[y, w, s, m, a, ]) / num_trips[y, w, m, a]) / mean(anglers_by_trip[y, w, m, a, ])
          observed_catch_rate_var[y, w, s, m, a] = 1 / (num_trips[y, w, m, a] * mean(anglers_by_trip[y, w, m, a, ]) ^ 2) * (var(observed_caught_by_trip[y, w, s, m, a, ]) + (sum(observed_total_caught[y, w, s, m, a, ]) / num_trips[y, w, m, a]) ^ 2 * var(anglers_by_trip[y, w, m, a, ]) - 2 * sum(observed_total_caught[y, w, s, m, a, ]) / num_trips[y, w, m, a] * cov(observed_caught_by_trip[y, w, s, m, a, ], anglers_by_trip[y, w, m, a, ]))
          #y_bar = mean(observed_caught_by_trip[y, w, s, m, a, ])#sum(observed_total_caught[y, w, s, m, a, ]) / num_trips[y, w, m, a]
          #x_bar = mean(anglers_by_trip[y, w, m, a, ])
          #var_y = var(observed_caught_by_trip[y, w, s, m, a, ])#var(observed_total_caught[y, w, s, m, a, ]) / (num_trips[y, w, m, a] ^ 2)
          #var_x = var(anglers_by_trip[y, w, m, a, ])
          #cov_xy = cov(observed_caught_by_trip[y, w, s, m, a, ], anglers_by_trip[y, w, m, a, ])
          
          #observed_catch_rate_var = (y_bar / x_bar) ^ 2 * (var_y / (y_bar ^ 2) + var_x / (x_bar ^ 2) - 2 * cov_xy / (x_bar * y_bar))
          
          catch_rate_var[y, w, s, m, a] = unavailable_catch_rate_var[y, w, s, m, a] + observed_catch_rate_var[y, w, s, m, a]
        }
      }
    }
  }
}

# effort

effort_df = read.csv("./MRIP_Ifiles_HI/mrip_effort_series.csv")
effort_df$Wave = ordered(effort_df$Wave, levels = c("JANUARY/FEBRUARY", "MARCH/APRIL", "MAY/JUNE", "JULY/AUGUST", "SEPTEMBER/OCTOBER", "NOVEMBER/DECEMBER"))
effort_df$Wave_num = as.numeric(effort_df$Wave)
effort_df$mode = ifelse(effort_df$Fishing.Mode == "PRIVATE/RENTAL BOAT", 1, ifelse(effort_df$Fishing.Mode == "SHORE", 2, NA))
effort_df$area = ifelse(effort_df$Fishing.Area == "OCEAN (> 3 MI)", 1, ifelse(effort_df$Fishing.Area == "OCEAN (<= 3 MI)", 2, ifelse(effort_df$Fishing.Area == "INLAND", 3, NA)))

effort = array(0, dim = c(n_years, n_waves, n_modes, n_areas))
effort_var = array(0, dim = c(n_years, n_waves, n_modes, n_areas))

for(y in 1:n_years) {
  for(w in 1:n_waves) {
    for(m in 1:n_modes) {
      for(a in 1:n_areas) {
        rows = effort_df[effort_df$Year == years[y] & effort_df$Wave_num == w & effort_df$mode == m & effort_df$area == a,]
        
        if(nrow(rows) > 0) {
          effort[y, w, m, a] = sum(rows$Angler.Trips)
          
          for(i in 1:nrow(rows)) {
            r = rows[i,]
            
            effort_var[y, w, m, a] = effort_var[y, w, m, a] + (r$Angler.Trips  * r$PSE / 100) ^ 2
          }
        }
      }
    }
  }
}

# total catch

total_catch_num = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))
total_catch_num_var = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))

for(s in 1:n_species) {
  total_catch_num[, , s, , ] = catch_rate[, , s, , ] * effort
  total_catch_num_var[, , s, , ] = catch_rate_var[, , s, , ] * effort ^ 2 + effort_var * catch_rate[, , s, , ] ^ 2 - catch_rate_var[, , s, , ] * effort_var
}

total_catch_num_var_annual = apply(total_catch_num_var, c(1, 3, 4, 5), sum)
total_catch_pse_annual = sqrt(total_catch_num_var_annual) / apply(total_catch_num, c(1, 3, 4, 5), sum) * 100

total_catch_weight = sweep(apply(total_catch_num, c(1, 3, 4, 5), sum), 2, mean_weight_all, FUN = "*")

proportion_kept = (observed_total_caught[, , , , , 2] + unavailable_total_caught[, , , , , 2]) / (apply(observed_total_caught + unavailable_total_caught, c(1, 2, 3, 4, 5), sum))
proportion_kept[is.nan(proportion_kept)] = 0

total_catch_num_kept = total_catch_num * proportion_kept
total_catch_weight_kept = sweep(apply(total_catch_num_kept, c(1, 3, 4, 5), sum), 2, mean_weight_all, FUN = "*")

#checking results

# checking catch disposition against Appendix Table 1c from Ma submitted manuscript
total_catch_proportion_kept = 1-apply(total_catch_num_kept, c(1, 3), sum, na.rm = T)/apply(total_catch_num, c(1, 3), sum, na.rm = T) # [year, species]

official_total_catch_proportion_kept = array(c(NA, 0, 0, NA, 11, 0, 100, 0, 100, 16, 15, 40, 74, 0, NA, NA, NA, NA, NA,
                                               NA, 71, 51, 89, 54, 38, 44, 53, 62, 59, 53, 69, 39, 66, NA, NA, NA, NA, NA,
                                               NA, 61, 31, 68, 68, 45, 46, 79, 77, 37, 37, 62, 10, 64, NA, NA, NA, NA, NA,
                                               NA, 69, 27, 23, 55, 56, 61, 78, 65, 57, 55, 83, 60, 69, NA, NA, NA, NA, NA,
                                               NA, 0, 0, 0, 5, 0, 0, 82, 71, 37, 48, 45, 75, 56, NA, NA, NA, NA, NA,
                                               NA, 100, 100, NA, 100, 0, 0, 0, 100, 0, 100, NA, 0, 40, NA, NA, NA, NA, NA,
                                               NA, 0, 0, NA, 60, 47, 0, 0, 0, 78, 0, 100, 82, 0, NA, NA, NA, NA, NA) / 100, dim = c(n_years, n_species))

for(s in 1:n_species) {
  plot(x = years, y = official_total_catch_proportion_kept[, s], type = "l", xlab = "Year", ylab = paste0(species_df[s,]$hawaiian_name, " Catch Number % Kept"))
  lines(x = years, y = total_catch_proportion_kept[, s], col = "red")
  legend(x = "topright", col = c("black", "red"), lty = 1, legend = c("Ma et al.", "Me"), bty = "n")
}

# checking catch number against https://www.fisheries.noaa.gov/data-tools/recreational-fisheries-statistics-queries

official_catch = read.csv("./MRIP_Ifiles_HI/HMRFS catch.csv")
official_catch$Total.Harvest..A.B1. = as.numeric(gsub(",", "", official_catch$Total.Harvest..A.B1.))
official_catch$wave = match(official_catch$Wave, c("JANUARY/FEBRUARY", "MARCH/APRIL", "MAY/JUNE", "JULY/AUGUST", "SEPTEMBER/OCTOBER", "NOVEMBER/DECEMBER"))

official_catch_num = array(0, dim = c(n_years, n_species))
official_catch_num_wave = array(0, dim = c(n_years, n_waves, n_species))
official_catch_num_var = array(0, dim = c(n_years, n_species))
official_catch_num_var_wave = array(0, dim = c(n_years, n_waves, n_species))

for(y in 1:n_years) {
  for(s in 1:n_species) {
    rows = official_catch[official_catch$Year == years[y] & official_catch$Common.Name == species_df$common_name[s],]
    
    if(nrow(rows) > 0) {
      official_catch_num[y, s] = sum(rows$Total.Harvest..A.B1.)
      
      for(i in 1:nrow(rows)) {
        r = rows[i,]
        
        official_catch_num_var[y, s] = official_catch_num_var[y, s] + (r$Total.Harvest..A.B1.  * r$PSE / 100) ^ 2
      }
      
      for(w in 1:n_waves) {
        rows = rows[rows$wave == w,]
        
        if(nrow(rows) > 0) {
          official_catch_num_wave[y, w, s] = sum(rows$Total.Harvest..A.B1.)
          
          for(i in 1:nrow(rows)) {
            official_catch_num_var_wave[y, w, s] = official_catch_num_var_wave[y, w, s] + (r$Total.Harvest..A.B1.  * r$PSE / 100) ^ 2
          }
        }
      }
    }
  }
}

for(s in 1:n_species) {
  plot(x = years, y = official_catch_num[, s], type = "l", xlab = "Year", ylab = paste0(species_df[s,]$hawaiian_name, " Catch"))
  lines(x = years, y = apply(total_catch_num[, , s, , ], c(1), sum), col = "red")
  legend(x = "topright", col = c("black", "red"), lty = 1, legend = c("HMRFS", "Me"), bty = "n")
}

# checking catch number variance against https://www.fisheries.noaa.gov/data-tools/recreational-fisheries-statistics-queries

sqrt(apply(total_catch_num_var_annual, c(1, 2), sum, na.rm = T))/sqrt(official_catch_num_var)

# Annual CPUE smoothing

# in the second-to-last dimension, index 1 = ocean (> 3 mi), index 2 = ocean (<= 3 mi), and index 3 = aggregated across fishing areas
# in the last dimension, index 1 = MA-3, index 2 = Kalman
year_smoothed_catch_rate = array(0, dim = c(n_years, n_species, 3, 2))
year_smoothed_total_catch = array(0, dim = c(n_years, n_species, 3, 2))

n_smooth_years = 1 # number of years on either side to smooth across
q = 0.01 # random walk parameter for Kalman filter
h = 0.01 # observation process parameter for Kalman filter

total_caught = observed_total_caught + unavailable_total_caught

for(s in 1:n_species) {
  year_smoothed_catch_rate[, s, 3, 1] = rollapplyr(apply(total_caught[, , s, 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum), 1 + 2 * n_smooth_years, FUN = mean, partial = T, align = "center")
  year_smoothed_total_catch[, s, 3, 1] = year_smoothed_catch_rate[, s, 3, 1] * apply(effort[, , 1, 1:2], c(1), sum)
  
  for(a in 1:(n_areas - 1)) {
    year_smoothed_catch_rate[, s, a, 1] = rollapplyr(apply(total_caught[, , s, 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum), 1 + 2 * n_smooth_years, FUN = mean, partial = T, align = "center")
    year_smoothed_total_catch[, s, a, 1] = year_smoothed_catch_rate[, s, a, 1] * apply(effort[, , 1, a], c(1), sum)
  }
  
  cpues = apply(total_caught[, , s, 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum)
  year_smoothed_catch_rate[, s, 3, 2] = KFS(SSModel(cpues ~ SSMtrend(1, Q = q), H = h))$alphahat
  year_smoothed_total_catch[, s, 3, 2] = year_smoothed_catch_rate[, s, 3, 2] * sapply(1:n_years, function(f) sum(effort[f, , 1, 1:2]))
  
  for(a in 1:(n_areas - 1)) {
    cpues = apply(total_caught[, , s, 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum)
    year_smoothed_catch_rate[, s, a, 2] = KFS(SSModel(cpues ~ SSMtrend(1, Q = q), H = h))$alphahat
    year_smoothed_total_catch[, s, a, 2] = year_smoothed_catch_rate[, s, a, 2] * sapply(1:n_years, function(f) sum(effort[f, , 1, a]))
  }
}

unaccounted_total_catch = apply(total_catch_num[, , , 2,], c(1, 3), sum, na.rm = T)
ma_smoothed_catch_plot = year_smoothed_total_catch[, , 3, 1] + unaccounted_total_catch
ma_smoothed_by_area_catch_plot = apply(year_smoothed_total_catch[, , 1:2, 1], c(1, 2), sum) + unaccounted_total_catch
kalman_smoothed_catch_plot = year_smoothed_total_catch[, , 3, 2] + unaccounted_total_catch
kalman_smoothed_by_area_catch_plot = apply(year_smoothed_total_catch[, , 1:2, 2], c(1, 2), sum) + unaccounted_total_catch

s = 6
y_max = max(c(ma_smoothed_catch_plot[, s], ma_smoothed_by_area_catch_plot[, s], kalman_smoothed_catch_plot[, s], kalman_smoothed_by_area_catch_plot[, s], official_catch_num[, s]))
plot(years, ma_smoothed_catch_plot[, s], ylim = c(0, y_max), col = "red", type = "l", xlab = "Year", ylab = paste0(species_df[s,]$hawaiian_name, " Catch"))
lines(years, ma_smoothed_by_area_catch_plot[, s], col = "red", lty = "dashed")
lines(years, kalman_smoothed_catch_plot[, s], col = "green")
lines(years, kalman_smoothed_by_area_catch_plot[, s], col = "green", lty = "dashed")
lines(years, official_catch_num[, s], col = "black")
legend(x = "topright", col = c("red", "red", "green", "green", "black"), lty = c("solid", "dashed", "solid", "dashed", "solid"), legend = c("MA-3", "MA-3 by area", "Kalman", "Kalman by area", "Original"), bty = "n")

# Wave-level CPUE smoothing

n_year_waves = n_years * n_waves

wave_smoothed_catch_rate = array(0, dim = c(n_year_waves, n_species, 3, 2))
wave_smoothed_total_catch = array(0, dim = c(n_year_waves, n_species, 3, 2))

n_smooth_waves = 6 # number of waves on either side to smooth across
q = 0.01 # random walk parameter for Kalman filter
h = 0.01 # observation process parameter for Kalman filter

total_caught = observed_total_caught + unavailable_total_caught

total_caught_yw = array(0, dim = c(n_year_waves, n_species, n_modes, n_areas, n_dispositions))
anglers_by_trip_yw = array(0, dim = c(n_year_waves, n_modes, n_areas, n_trips))
effort_yw = array(0, dim = c(n_year_waves, n_modes, n_areas))

for(y in 1:n_years) {
  for(w in 1:n_waves) {
    index = (y - 1) * 6 + w
    
    total_caught_yw[index, , , , ] = total_caught[y, w, , , , ]
    anglers_by_trip_yw[index, , , ] = anglers_by_trip[y, w, , , ]
    effort_yw[index, , ] = effort[y, w, , ]
  }
}

for(s in 1:n_species) {
  wave_smoothed_catch_rate[, s, 3, 1] = rollapplyr(apply(total_caught_yw[, s, 1, 1:2, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, 1:2, ], c(1), sum), 1 + 2 * n_smooth_waves, FUN = mean, partial = T, align = "center")
  wave_smoothed_total_catch[, s, 3, 1] = wave_smoothed_catch_rate[, s, 3, 1] * apply(effort_yw[, 1, 1:2], c(1), sum)
  
  for(a in 1:(n_areas - 1)) {
    wave_smoothed_catch_rate[, s, a, 1] = rollapplyr(apply(total_caught_yw[, s, 1, a, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, a, ], c(1), sum), 1 + 2 * n_smooth_waves, FUN = mean, partial = T, align = "center")
    wave_smoothed_total_catch[, s, a, 1] = wave_smoothed_catch_rate[, s, a, 1] * effort_yw[, 1, a]
  }
  
  cpues = apply(total_caught_yw[, s, 1, 1:2, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, 1:2, ], c(1), sum)
  wave_smoothed_catch_rate[, s, 3, 2] = KFS(SSModel(cpues ~ SSMtrend(1, Q = q), H = h))$alphahat
  wave_smoothed_total_catch[, s, 3, 2] = wave_smoothed_catch_rate[, s, 3, 2] * sapply(1:n_year_waves, function(f) sum(effort_yw[f, 1, 1:2]))
  
  for(a in 1:(n_areas - 1)) {
    cpues = apply(total_caught_yw[, s, 1, a, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, a, ], c(1), sum)
    wave_smoothed_catch_rate[, s, a, 2] = KFS(SSModel(cpues ~ SSMtrend(1, Q = q), H = h))$alphahat
    wave_smoothed_total_catch[, s, a, 2] = wave_smoothed_catch_rate[, s, a, 2] * sapply(1:n_year_waves, function(f) sum(effort_yw[f, 1, a]))
  }
}

year_smoothed_total_catch = array(0, dim = c(n_years, n_species, 3, 2))

for(s in 1:n_species) {
  for(area in 1:3) {
    for(smoother in 1:2) {
      year_smoothed_total_catch[, s, area, smoother] = sapply(1:n_years, function(f) sum(wave_smoothed_total_catch[((f - 1) * 6 + 1):((f - 1) * 6 + 6), s, area, smoother], na.rm = T))
    }
  }
}

unaccounted_total_catch = apply(total_catch_num[, , , 2,], c(1, 3), sum, na.rm = T)
ma_smoothed_catch_plot = year_smoothed_total_catch[, , 3, 1] + unaccounted_total_catch
ma_smoothed_by_area_catch_plot = apply(year_smoothed_total_catch[, , 1:2, 1], c(1, 2), sum) + unaccounted_total_catch
kalman_smoothed_catch_plot = year_smoothed_total_catch[, , 3, 2] + unaccounted_total_catch
kalman_smoothed_by_area_catch_plot = apply(year_smoothed_total_catch[, , 1:2, 2], c(1, 2), sum) + unaccounted_total_catch

s = 4
y_max = max(c(ma_smoothed_catch_plot[, s], ma_smoothed_by_area_catch_plot[, s], kalman_smoothed_catch_plot[, s], kalman_smoothed_by_area_catch_plot[, s], official_catch_num[, s]))
plot(years, ma_smoothed_catch_plot[, s], ylim = c(0, y_max), col = "red", type = "l", xlab = "Year", ylab = paste0(species_df[s,]$hawaiian_name, " Catch"))
lines(years, ma_smoothed_by_area_catch_plot[, s], col = "red", lty = "dashed")
lines(years, kalman_smoothed_catch_plot[, s], col = "green")
lines(years, kalman_smoothed_by_area_catch_plot[, s], col = "green", lty = "dashed")
lines(years, official_catch_num[, s], col = "black")
legend(x = "topright", col = c("red", "red", "green", "green", "black"), lty = c("solid", "dashed", "solid", "dashed", "solid"), legend = c("MA-3", "MA-3 by area", "Kalman", "Kalman by area", "Original"), bty = "n")

t2 = Sys.time()
t2 - t1