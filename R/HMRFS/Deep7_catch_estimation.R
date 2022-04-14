library(sas7bdat)

setwd("/Users/Toby/Documents/GitHub/MHI_Bottomfish_2023/R/HMRFS")
i1 = read.sas7bdat("./MRIP_Ifiles_HI/i1.sas7bdat", debug=TRUE)
observed_catch = read.sas7bdat("./MRIP_Ifiles_HI/i3.sas7bdat", debug=TRUE)
unavailable_catch = read.sas7bdat("./MRIP_Ifiles_HI/i2.sas7bdat", debug=TRUE)

species_df = data.frame("key" = c(8835020413, 8835360304, 8835360302, 8835360704, 8835360706, 8835360901, 8835360707),
                        "common_name" = c("Hawaiian Grouper", "Long-tail Red Snapper", "Ruby Snapper", "Pink Snapper", "Von Siebold's Snapper", "Ironjaw Snapper", "Brigham's Snapper"),
                        "hawaiian_name" = c("Hapu'upu'u", "Onaga", "Ehu", "Opakapaka", "Kalekale", "Lehi", "Gindai"),
                        "alpha" = c(2.884, 2.673, 3.026, 2.928, 2.932, 2.458, 2.859), # cm
                        "beta" = c(3.065*10^-5, 6.005*10^-5, 1.551*10^-5, 2.311*10^-5, 2.243*10^-5, 1.298*10^-4, 3.526*10^-5))

# for mode: 1 = private boat, 2 = shore
# for disposition: index 1 = sold, 2 = not sold
i1$mode = ifelse(i1$MODE_F == 8, 1, ifelse(i1$MODE_FX %in% 1:5, 2, NA))
observed_catch = observed_catch[observed_catch$DISP3 %in% c(3, 4, 5) & observed_catch$MODE_F != 7,] # fish that were not released and not from a charter trip
observed_catch$sold = ifelse(observed_catch$DISP3 == 5, T, F)
unavailable_catch = unavailable_catch[unavailable_catch$DISPO %in% c(3, 4, 5) & unavailable_catch$MODE_F != 7,] # fish that were not released and not from a charter trip
unavailable_catch$sold = ifelse(unavailable_catch$DISPO == 5, T, F)

observed_catch$EST_WGT = sapply(1:nrow(observed_catch), function(r) {
  sp = observed_catch[r,]$SP_CODE
  if((sp %in% species_df$key) && !is.nan(observed_catch[r,]$LNGTH) && is.nan(observed_catch[r,]$WGT)) { # bottomfish with length, but no weight
    i = which(species_df$key == observed_catch[r,]$SP_CODE)
    return(species_df[i,]$beta * (observed_catch[r,]$LNGTH / 10) ^ species_df[i,]$alpha)
  } else {
    return(NaN)
  }
})

years = sort(unique(c(observed_catch$YEAR, unavailable_catch$year)), decreasing = F)
trips = unique(i1[!is.na(i1$mode) & !is.na(i1$ID_CODE),]$ID_CODE)

n_waves = 6
n_modes = 2
n_dispositions = 2

num_weighed = array(0, dim = c(length(years), n_waves, nrow(species_df), n_modes))
num_estimated = array(0, dim = c(length(years), n_waves, nrow(species_df), n_modes))
num_caught_by_trip = array(0, dim = c(length(years), n_waves, nrow(species_df), n_modes, length(trips)))
anglers_by_trip = array(0, dim = c(length(years), n_waves, n_modes, length(trips)))
total_weight = array(0, dim = c(length(years), n_waves, nrow(species_df), n_modes))
observed_total_caught = array(0, dim = c(length(years), n_waves, nrow(species_df), n_modes, n_dispositions))
unavailable_total_caught = array(0, dim = c(length(years), n_waves, nrow(species_df), n_modes, n_dispositions))
mean_weight_wave = array(0, dim = c(length(years), n_waves, nrow(species_df), n_modes))
mean_weight_year = array(0, dim = c(length(years), nrow(species_df), n_modes))
num_trips = array(0, dim = c(length(years), n_waves, n_modes))
catch_rate = array(0, dim = c(length(years), n_waves, nrow(species_df), n_modes))
catch_rate_var = array(0, dim = c(length(years), n_waves, nrow(species_df), n_modes))

for(y in 1:length(years)) {
  for(w in 1:n_waves) {
    for(m in 1:n_modes) {
      num_trips[y, w, m] = length(unique((i1[i1$YEAR == years[y] & i1$WAVE == w & i1$mode == m & !is.na(i1$ID_CODE),])$ID_CODE))
    }
  }
}

for(i in 1:length(trips)) {
  trip = trips[i]
  year = as.numeric(substr(trip, 6, 9))
  y = which(years == year)
  w = ceiling(as.numeric(substr(trip, 10, 11)) / 2)
  m = i1[i1$ID_CODE == trip,]$mode
  
  anglers_by_trip[y, w, m, i] = i1[i1$ID_CODE == trip,]$CNTRBTRS
  
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
      for(j in 1:length(species)) {
        s = which(species_df$key == species[j])
        
        observed_s = observed[observed$SP_CODE == species[j],]
        
        num_caught_by_trip[y, w, s, m, i] = num_caught_by_trip[y, w, s, m, i] + max(observed_s$FSHINSP)
        
        if(all(observed_s$sold)) {
          observed_total_caught[y, w, s, m, 1] = observed_total_caught[y, w, s, m, 1] + max(observed_s$FSHINSP)
        } else if(all(!observed$sold)) {
          observed_total_caught[y, w, s, m, 2] = observed_total_caught[y, w, s, m, 2] + max(observed_s$FSHINSP)
        } else {
          # there are 7 instances where there are multiple dispositions within the same interview/species
          # in each case there is a row for each fish caught, but this could potentially not be the case
          # e.g. ID_CODE = 1700920031105001, SP_CODE = 8835360704
          count_sold = 0
          count_unsold = 0
          
          for(k in 1:nrow(observed_s)) {
            r = observed_s[k,]
            if(r$sold) {
              count_sold = count_sold + 1
            } else {
              count_unsold = count_unsold + 1
            }
            
            observed_total_caught[y, w, s, m, 1] = observed_total_caught[y, w, s, m, 1] + max(observed_s$FSHINSP) * count_sold / (count_sold + count_unsold)
            observed_total_caught[y, w, s, m, 2] = observed_total_caught[y, w, s, m, 2] + max(observed_s$FSHINSP) * count_unsold / (count_sold + count_unsold)
          }
        }
        
        for(k in 1:nrow(observed_s)) {
          r = observed_s[k,]
          
          if(!is.nan(r$WGT)) {
            num_weighed[y, w, s, m] = num_weighed[y, w, s, m] + 1
            total_weight[y, w, s, m] = total_weight[y, w, s, m] + r$WGT
          } else if(!is.nan(r$LNGTH)) {
            num_estimated[y, w, s, m] = num_estimated[y, w, s, m] + 1
            total_weight[y, w, s, m] = total_weight[y, w, s, m] + r$EST_WGT
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
      for(j in 1:length(species)) {
        s = which(species_df$key == species[j])
        
        r = unavailable[unavailable$SP_CODE == species[j],]
        
        # There is an instance of two rows for the same specie within a single trip (ID_CODE = 1701620190402001, SP_CODE = 8835360704). The below code counts both rows, each with 7 fish.
        for(k in 1:nrow(r)) {
          q = r[k,]
          num_caught_by_trip[y, w, s, m, i] = num_caught_by_trip[y, w, s, m, i] + q$NUM_FISH
          
          if(q$sold) {
            unavailable_total_caught[y, w, s, m, 1] = unavailable_total_caught[y, w, s, m, 1] + q$NUM_FISH
          } else {
            unavailable_total_caught[y, w, s, m, 2] = unavailable_total_caught[y, w, s, m, 2] + q$NUM_FISH
          }
        }
      }
    }
  }
}

for(i in 1:nrow(species_df)) {
  mean_weight_wave[, , i, ] = total_weight[, , i, ] / (num_weighed[, , i, ] + num_estimated[, , i, ])
  mean_weight_year[, i, ] = apply(total_weight[, , i, ], c(1, 3), sum) / apply(num_weighed[, , i, ] + num_estimated[, , i, ], c(1, 3), sum)
  catch_rate[, , i, ] = apply(observed_total_caught[, , i, , ] + unavailable_total_caught[, , i, , ], c(1, 2, 3), sum) / apply(anglers_by_trip, c(1, 2, 3), sum)
}
catch_rate[is.nan(catch_rate)] = 0

for(i in 1:length(years)) {
  for(j in 1:n_waves) {
    for(k in 1:nrow(species_df)) {
      for(l in 1:n_modes) {
        x = num_caught_by_trip[i, j, k, l, ]
        y = anglers_by_trip[i, j, l, ]
        
        catch_rate_var[i, j, k, l] = (mean(x) / mean(y)) ^ 2 * (var(x) / mean(x) ^ 2 - 2 * cov(x, y) / (mean(x) * mean(y)) + var(y) / (mean(y)) ^ 2)
      }
    }
  }
}

# effort

effort_df = read.csv("./MRIP_Ifiles_HI/mrip_effort_series.csv")
effort_df$Wave = ordered(effort_df$Wave, levels = c("JANUARY/FEBRUARY", "MARCH/APRIL", "MAY/JUNE", "JULY/AUGUST", "SEPTEMBER/OCTOBER", "NOVEMBER/DECEMBER"))
effort_df$Wave_num = as.numeric(effort_df$Wave)
effort_df$mode = ifelse(effort_df$Fishing.Mode == "PRIVATE/RENTAL BOAT", 1, ifelse(effort_df$Fishing.Mode == "SHORE", 2, NA))

effort = array(0, dim = c(length(years), n_waves, n_modes))
effort_var = array(0, dim = c(length(years), n_waves, n_modes))
for(i in 1:length(years)) {
  for(j in 1:n_waves) {
    for(k in 1:n_modes) {
      entry = effort_df[effort_df$Year == years[i] & effort_df$Wave_num == j & effort_df$mode == k,]
      effort[i, j, k] = entry$Angler.Trips
      effort_var[i, j, k] = (entry$Angler.Trips  * entry$PSE / 100) ^ 2
    }
  }
}

# total catch

total_catch_num = array(0, dim = c(length(years), n_waves, nrow(species_df), n_modes))
total_catch_num_var = array(0, dim = c(length(years), n_waves, nrow(species_df), n_modes))
for(i in 1:nrow(species_df)) {
  total_catch_num[, , i, ] = catch_rate[, , i, ] * effort
  total_catch_num_var[, , i, ] = catch_rate_var[, , i, ] * effort ^ 2 + effort_var * catch_rate[, , i, ] ^ 2 - catch_rate_var[, , i, ] * effort_var
}

total_catch_num_var_annual = apply(total_catch_num_var, c(1, 3, 4), sum)
total_catch_pse_annual = sqrt(total_catch_num_var_annual) / apply(total_catch_num, c(1, 3, 4), sum) * 100

total_catch_weight = total_catch_num * mean_weight_wave

proportion_sold = (observed_total_caught[, , , , 1] + unavailable_total_caught[, , , , 1]) / (apply(observed_total_caught + unavailable_total_caught, c(1, 2, 3, 4), sum))
total_sale_num = total_catch_num * proportion_sold
total_sale_weight = total_catch_weight * proportion_sold

# checking results

true_annual_catch = t(matrix(c(333, 406, 1473, 0, 2306, 1072, 335, 1342, 334, 1425, 1753, 1955, 718, 91, 0 ,1901, 3299, 1484, 0, 13878, 13845, 18263, 29913, 7482, 7984, 4359, 19992, 7607, 18519, 14327, 17909, 9917, 5515, 1302, 26473, 1513, 6275, 5902, 0, 12552, 9543, 13196, 7701, 12150, 13575, 49693, 7021, 22016, 30148, 32939, 9178, 8829, 3869, 10921, 30043, 12902, 10856, 25380, 83909, 24247, 27653, 16845, 42544, 30194, 115003, 25380, 46702, 45120, 41332, 22924, 17082, 9739, 56140, 62815, 16486, 12997, 1157, 406, 2783, 5331, 4912, 6098, 819, 28983, 8517, 15739, 7598, 7596, 5650, 9443, 4347, 21872, 12596, 20232, 12789, 0, 10021, 1761, 0, 568, 1855, 697, 417, 271, 440, 335, 0, 1416, 721, 192, 0, 1245, 313, 895, 0, 813, 3674, 0, 850, 8037, 484, 2110, 167, 9254, 2147, 166, 1708, 238, 920, 2151, 20394, 4999, 2133), nrow = length(years), ncol = nrow(species_df), byrow = F))

for(i in 1:nrow(species_df)) {
  plot(x = years, y = true_annual_catch[i, ], type = "l")
  lines(x = years, y = apply(total_catch_num[, , i, 1], c(1), sum), col = "red")
}