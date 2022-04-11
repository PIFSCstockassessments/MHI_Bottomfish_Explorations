# are there ever separate dispositions within species for an interview?

library(sas7bdat)

observed_catch = read.sas7bdat("./MRIP_Ifiles_HI/i3.sas7bdat", debug=TRUE)
unavailable_catch = read.sas7bdat("./MRIP_Ifiles_HI/i2.sas7bdat", debug=TRUE)

species_df = data.frame("key" = c(8835020413, 8835360304, 8835360302, 8835360704, 8835360706, 8835360901, 8835360707),
                        "common_name" = c("Hawaiian Grouper", "Long-tail Red Snapper", "Ruby Snapper", "Pink Snapper", "Von Siebold's Snapper", "Ironjaw Snapper", "Brigham's Snapper"),
                        "hawaiian_name" = c("Hapu'upu'u", "Onaga", "Ehu", "Opakapaka", "Kalekale", "Lehi", "Gindai"),
                        "alpha" = c(2.884, 2.673, 3.026, 2.928, 2.932, 2.458, 2.859), # cm
                        "beta" = c(3.065*10^-5, 6.005*10^-5, 1.551*10^-5, 2.311*10^-5, 2.243*10^-5, 1.298*10^-4, 3.526*10^-5))

observed_catch$sold = ifelse(observed_catch$DISP3 == 5, T, F)
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
trips = unique(c(observed_catch$ID_CODE, unavailable_catch$ID_CODE))

# for disposition: index 1 = sold, 2 = not sold
num_weighed = array(0, dim = c(length(years), 6, nrow(species_df)))
num_estimated = array(0, dim = c(length(years), 6, nrow(species_df)))
num_caught_by_trip = array(0, dim = c(length(years), 6, nrow(species_df), length(trips)))
total_weight = array(0, dim = c(length(years), 6, nrow(species_df)))
observed_total_caught = array(0, dim = c(length(years), 6, nrow(species_df), 2))
unavailable_total_caught = array(0, dim = c(length(years), 6, nrow(species_df), 2))
mean_weight_wave = array(0, dim = c(length(years), 6, nrow(species_df)))
mean_weight_year = array(0, dim = c(length(years), nrow(species_df)))
num_trips = array(0, dim = c(length(years), 6))
catch_rate = array(0, dim = c(length(years), 6, nrow(species_df)))
catch_rate_var = array(0, dim = c(length(years), 6, nrow(species_df)))

for(i in 1:length(trips)) {
  trip = trips[i]
  year = as.numeric(substr(trip, 6, 9))
  y = which(years == year)
  w = ceiling(as.numeric(substr(trip, 10, 11)) / 2)
  
  num_trips[y, w] = num_trips[y, w] + 1
  
  # observed catch
  # FSHINSP is the number of specimens of a species observed
  #   Each row represents one specimen measured, but all have the same FSHINSP
  
  observed = observed_catch[observed_catch$ID_CODE == trip,]
  
  if(length(observed) > 0) {
    species = unique(observed$SP_CODE)
    species = species[species %in% species_df$key]
    
    if(length(species) > 0) {
      for(j in 1:length(species)) {
        s = which(species_df$key == species[j])
        
        observed_s = observed[observed$SP_CODE == species[j],]
        
        num_caught_by_trip[y, w, s, i] = num_caught_by_trip[y, w, s, i] + max(observed_s$FSHINSP)
        
        if(any(observed_s$sold)) {
          observed_total_caught[y, w, s, 1] = observed_total_caught[y, w, s, 1] + max(observed_s$FSHINSP)
        } else {
          observed_total_caught[y, w, s, 2] = observed_total_caught[y, w, s, 2] + max(observed_s$FSHINSP)
        }
        
        for(k in 1:nrow(observed_s)) {
          r = observed_s[k,]
          
          if(!is.nan(r$WGT)) {
            num_weighed[y, w, s] = num_weighed[y, w, s] + 1
            total_weight[y, w, s] = total_weight[y, w, s] + r$WGT
          } else if(!is.nan(r$LNGTH)) {
            num_estimated[y, w, s] = num_estimated[y, w, s] + 1
            total_weight[y, w, s] = total_weight[y, w, s] + r$EST_WGT
          }
        }
      }
    }
    
    # unavailable catch
    # NUM2: sequential numbering of species within interview
    # NUM_TYP2: total number of species within interview
    # NUM_FISH: number of of fish of each species
    
    unavailable = unavailable_catch[unavailable_catch$ID_CODE == trip,]
    
    species = unique(unavailable$SP_CODE)
    species = species[species %in% species_df$key]
      
    if(length(species) > 0) {
      for(j in 1:length(species)) {
        s = which(species_df$key == species[j])
          
        r = unavailable[unavailable$SP_CODE == species[j],]
          
        #for(k in 1:nrow(r)) {
        #  q = r[k,]
        #  num_caught_by_trip[y, w, s, i] = num_caught_by_trip[y, w, s, i] + q$NUM_FISH
          
        #  if(q$sold) {
        #    unavailable_total_caught[y, w, s, 1] = unavailable_total_caught[y, w, s, 1] + q$NUM_FISH
        #  } else {
        #    unavailable_total_caught[y, w, s, 2] = unavailable_total_caught[y, w, s, 2] + q$NUM_FISH
        #  }
        #}
        # There are a few instances of multiple rows for the same specie within a single trip. The below code only counts one row (2531 unavailable vs 2538 with the above)
        num_caught_by_trip[y, w, s, i] = num_caught_by_trip[y, w, s, i] + max(r$NUM_FISH)
          
        if(any(r$sold)) {
          unavailable_total_caught[y, w, s, 1] = unavailable_total_caught[y, w, s, 1] + max(r$NUM_FISH)
        } else {
          unavailable_total_caught[y, w, s, 2] = unavailable_total_caught[y, w, s, 2] + max(r$NUM_FISH)
        }
      }
    }
  }
}

for(i in 1:nrow(species_df)) {
  mean_weight_wave[, , i] = total_weight[, , i] / (num_weighed[, , i] + num_estimated[, , i])
  mean_weight_year[, i] = apply(total_weight[, , i], c(1), sum) / apply(num_weighed[, , i] + num_estimated[, , i], c(1), sum)
  catch_rate[, , i] = apply(observed_total_caught[, , i, ] + unavailable_total_caught[, , i, ], c(1, 2), sum) / num_trips
}

for(i in 1:length(years)) {
  for(j in 1:6) {
    for(k in 1:nrow(species_df)) {
      catch_rate_var[i, j, k] = var(num_caught_by_trip[i, j, k, ] / length(num_caught_by_trip[i, j, k, ]))
    }
  }
}


# effort

effort_df = read.csv("./MRIP_Ifiles_HI/mrip_effort_series.csv")
effort_df$Wave = ordered(effort_df$Wave, levels = c("JANUARY/FEBRUARY", "MARCH/APRIL", "MAY/JUNE", "JULY/AUGUST", "SEPTEMBER/OCTOBER", "NOVEMBER/DECEMBER"))
effort_df$Wave_num = as.numeric(effort_df$Wave)
effort_df = effort_df[effort_df$Fishing.Mode == "PRIVATE/RENTAL BOAT",]

effort = array(0, dim = c(length(years), 6))
effort_var = array(0, dim = c(length(years), 6))
for(i in 1:length(years)) {
  for(j in 1:6) {
    entry = effort_df[effort_df$Year == years[i] & effort_df$Wave_num == j,]
    effort[i, j] = entry$Angler.Trips
    effort_var[i, j] = (entry$Angler.Trips * entry$PSE / 100) ^ 2
  }
}

total_catch_num = array(0, dim = c(length(years), 6, nrow(species_df)))
total_catch_var = array(0, dim = c(length(years), 6, nrow(species_df)))
for(i in 1:nrow(species_df)) {
  total_catch_num[, , i] = catch_rate[, , i] * effort
  total_catch_var[, , i] = catch_rate_var[, , i] * effort ^ 2 + effort_var * catch_rate[, , i] ^ 2 - catch_rate_var[, , i] * effort_var
}




# catch rate and effort are both calculated by wave --> total catch number; multiply by mean weight in the wave for total catch weight
# catch estimates from 2003-2010 in the catch estimation files were adjusted by a factor of 82% (= 1/1.22) for an error in the population household count for Maui County
# The proportions of fish sold or non-sold (Pw) and the catch number estimates (Cw) from each wave (w) in a year were combined to estimate the proportion of fish sold or non-sold
# "We obtained catch in weight by multiplying the number of uku caught in a year by the overall mean weight of uku in the HMRFS dataset between 2003 and 2018 (mean of 3.26 kg from 151 individuals)"
# "We obtained the proportion of uku catch sold by year from HMRFS fisher interviews and multiplied the total HMRFS uku catch by one minus the year-specific proportion"

# Figure 1: Catch number (total catch and catch not sold) by species and year, with SE on total catch number
# Figure 2: Catch weight (total catch and catch not sold) by species and year, with SE on total catch weight
# Table 1: Proportion of estimated catch numbers without catch weight in the default HMRFS catch weight estimates, by species and year
# Table 2: Mean weight by species and year, plus the number of measured and calculated weight values

bottomfish_lw = observed_catch[(observed_catch$SP_CODE %in% species_df$key) & (!is.nan(observed_catch$LNGTH) | !is.nan(observed_catch$WGT)),]
specie_lw = bottomfish_lw[bottomfish_lw$SP_CODE == 8835360704,]
summarise(group_by(specie_lw, YEAR), total = n(), weighed = sum(ifelse(!is.nan(WGT), 1, 0)), length_only = sum(ifelse(!is.nan(EST_WGT), 1, 0)))

# Table 3: Mean weight, number of weight records (measured or estimated?), and SD by year and catch disposition (sold, unsold, combined) for opakapaka
# Table 4: Mean weight by wave (across all years), SD, and number of weight records (measured or estimated) for opakapaka, onaga, ehu, and kalekale