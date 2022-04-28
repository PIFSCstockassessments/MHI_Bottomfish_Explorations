num_weighed = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))
num_estimated = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))
n = num_weighed + num_estimated

total_weight = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))

mean_weight_wave = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))
mean_weight_year = array(0, dim = c(n_years, n_species, n_modes, n_areas))
# for mode: 1 = private boat, 2 = shore
# for area: 1 = ocean (> 3 mi), 2 = ocean (<= 3 mi), 3 = inland

species_df = data.frame("key" = c(8835020413, 8835360304, 8835360302, 8835360704, 8835360706, 8835360901, 8835360707),
                        "common_name" = c("HAWAIIAN GROUPER", "LONGTAILED RED SNAPPER", "RUBY SNAPPER", "PINK SNAPPER", "VON SIEBOLDS SNAPPER", "IRONJAW SNAPPER", "BINGHAMS SNAPPER"),
                        "hawaiian_name" = c("Hapu'upu'u", "Onaga", "Ehu", "Opakapaka", "Kalekale", "Lehi", "Gindai"),
                        "alpha" = c(2.884, 2.673, 3.026, 2.928, 2.932, 2.458, 2.859), # cm
                        "beta" = c(3.065*10^-5, 6.005*10^-5, 1.551*10^-5, 2.311*10^-5, 2.243*10^-5, 1.298*10^-4, 3.526*10^-5))


mean_weight_wave[, , s, , ] = total_weight[, , s, , ] / (num_weighed[, , s, , ] + num_estimated[, , s, , ])
mean_weight_year[, s, , ] = apply(total_weight[, , s, , ], c(1, 3, 4), sum) / apply(num_weighed[, , s, , ] + num_estimated[, , s, , ], c(1, 3, 4), sum)




observed_total_caught = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas, n_dispositions))
unavailable_total_caught = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas, n_dispositions))

a = observed_total_caught + unavailable_total_caught
a = apply(a, c(1, 3, 6), sum)
a[, , 1] / apply(a, c(1, 2), sum)


# Weight data points by species/year

cols = c("black", "deeppink", "deepskyblue", "blueviolet", "lawngreen")
transparent_cols = array("", dim = c(length(cols)))
for(i in 1:length(cols)) {
  col_rgb = col2rgb(cols[i])
  transparent_cols[i] = rgb(col_rgb[1] / 255, col_rgb[2] / 255, col_rgb[3] / 255, alpha = 0.7)
}

x_min = min(years)
x_max = max(years)

par(xpd = NA, mar = c(4, 4, 3, 3))

for(s in 1:nrow(species_df)) {
  observed_s = observed_catch[observed_catch$SP_CODE == species_df$key[s],]
  w_s = data.frame(WGT = c(observed_s$WGT, observed_s$EST_WGT), YEAR = c(observed_s$YEAR, observed_s$YEAR))
  w_s = w_s[!is.na(w_s$WGT),]
  
  y_max = ceiling(max(w_s$WGT) * 1.2)
  y_min = 0
  
  plot(NA, axes = F, xaxs = "i", yaxs = "i", xlim = c(x_min, x_max), ylim = c(y_min, y_max), xlab = "", ylab = "")
  axis(side = 1, at = seq(from = 2005, to = 2020, by = 5))
  lines(x = c(x_min, x_max), y = c(y_min, y_min))
  axis(side = 2, at = pretty(y_min:y_max, n = 4), las = 2)
  lines(x = c(x_min, x_min), y = c(y_min, y_max))
  mtext("Year", side = 1, line = 3)
  mtext("Weight (kg)", side = 2, line = 3)
  mtext(species_df$hawaiian_name[s], side = 3, line = 1)
  
  w_sy_mean = array(NA, dim = c(length(years)))
  w_sy_sd = array(NA, dim = c(length(years)))
  
  for(y in 1:length(years)) {
    w_sy = w_s[w_s$YEAR == years[y],]
    
    if(nrow(w_sy) > 0) {
      points(x = rep(years[y], nrow(w_sy)), y = w_sy$WGT, col = rgb(col2rgb(cols[3])[1] / 255, col2rgb(cols[3])[2] / 255, col2rgb(cols[3])[3] / 255, alpha = min(1, 20 / nrow(w_s))), bg = NA, pch = 19)
      
      w_sy_mean[y] = mean(w_sy$WGT, na.rm = T)
      
      if(nrow(w_sy) > 1) {
        w_sy_sd[y] = sqrt(var(w_sy$WGT, na.rm = T))
        lines(x = c(years[y], years[y]), y = c(w_sy_mean[y] - w_sy_sd[y], w_sy_mean[y] + w_sy_sd[y]), col = rgb(0, 0, 0, alpha = 0.5))
      }
    }
  }
  
  lines(x = years[!is.na(w_sy_mean)], y = w_sy_mean[!is.na(w_sy_mean)], col = cols[1])
  
  w_s_mean = mean(w_s$WGT, na.rm = T)
  lines(x = c(x_min, x_max), y = c(w_s_mean, w_s_mean), lty = "dashed", col = transparent_cols[1])
  
  text(x = x_max, y = w_s_mean, paste0(" n = ", nrow(w_s)), adj = c(0, 0.5))
  
  legend(x = "topright", legend = c("Data Point", "Annual Mean", "Overall Mean"), col = c(transparent_cols[3], "black", transparent_cols[1]), lty = c(NA, "solid", "dashed"), pch = c(19, NA, NA), bty = "n", cex = 0.75)
}

# CPUE by species/year


catch_rate = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas))
observed_total_caught = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas, n_dispositions))
unavailable_total_caught = array(0, dim = c(n_years, n_waves, n_species, n_modes, n_areas, n_dispositions))
anglers_by_trip = array(0, dim = c(n_years, n_waves, n_modes, n_areas, n_trips))

catch_rate[, , s, , ] = apply(observed_total_caught[, , s, , , ] + unavailable_total_caught[, , s, , , ], c(1, 2, 3, 4), sum) / apply(anglers_by_trip, c(1, 2, 3, 4), sum)