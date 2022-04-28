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
  mtext("Weight ± SE (kg)", side = 2, line = 3)
  mtext(species_df$hawaiian_name[s], side = 3, line = 1)
  
  w_sy_mean = array(NA, dim = c(length(years)))
  w_sy_se = array(NA, dim = c(length(years)))
  
  for(y in 1:length(years)) {
    w_sy = w_s[w_s$YEAR == years[y],]
    
    if(nrow(w_sy) > 0) {
      points(x = rep(years[y], nrow(w_sy)), y = w_sy$WGT, col = rgb(col2rgb(cols[3])[1] / 255, col2rgb(cols[3])[2] / 255, col2rgb(cols[3])[3] / 255, alpha = min(1, 20 / nrow(w_s))), bg = NA, pch = 19)
      
      w_sy_mean[y] = mean(w_sy$WGT, na.rm = T)
      
      if(nrow(w_sy) > 1) {
        w_sy_se[y] = sqrt(var(w_sy$WGT[!is.na(w_sy$WGT)]) / length(w_sy$WGT[!is.na(w_sy$WGT)]))
        lines(x = c(years[y], years[y]), y = c(w_sy_mean[y] - w_sy_se[y], w_sy_mean[y] + w_sy_se[y]), col = rgb(0, 0, 0, alpha = 0.5))
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

cols = c("black", "deeppink", "deepskyblue", "blueviolet", "lawngreen")
transparent_cols = array("", dim = c(length(cols)))
for(i in 1:length(cols)) {
  col_rgb = col2rgb(cols[i])
  transparent_cols[i] = rgb(col_rgb[1] / 255, col_rgb[2] / 255, col_rgb[3] / 255, alpha = 0.7)
}

x_min = min(years)
x_max = max(years)

par(xpd = NA, mar = c(4, 4, 3, 1))

for(s in 1:nrow(species_df)) {
  catch_rate_sy = sapply(1:length(years), function(f) sum(observed_total_caught[f, , s, , , ] + unavailable_total_caught[f, , s, , , ]) / sum(anglers_by_trip[f, , , ,]))
  catch_rate_se_sy = sapply(1:length(years), function(f) sqrt(sum(catch_rate_var[f, , s, , ], na.rm = T) / sum(num_trips[f, , , ])))
  
  y_ticks = pretty(c(0, max(catch_rate_sy)), n = 4)
  y_max = max(y_ticks)
  y_min = 0
  
  plot(NA, axes = F, xaxs = "i", yaxs = "i", xlim = c(x_min, x_max), ylim = c(y_min, y_max), xlab = "", ylab = "")
  axis(side = 1, at = seq(from = 2005, to = 2020, by = 5))
  lines(x = c(x_min, x_max), y = c(y_min, y_min))
  axis(side = 2, at = y_ticks[y_ticks >= 0], las = 2)
  lines(x = c(x_min, x_min), y = c(y_min, y_max))
  mtext("Year", side = 1, line = 3)
  mtext("CPUE (# / Angler Trip)", side = 2, line = 3)
  mtext(species_df$hawaiian_name[s], side = 3, line = 1)
  
  w_sy_mean = array(NA, dim = c(length(years)))
  w_sy_sd = array(NA, dim = c(length(years)))
  
  points(x = years[!is.na(catch_rate_sy)], y = catch_rate_sy[!is.na(catch_rate_sy)], col = cols[3])
  
  for(y in 1:length(years)) {
    lines(x = c(years[y], years[y]), y = c(catch_rate_sy[y] - catch_rate_se_sy[y], catch_rate_sy[y] + catch_rate_se_sy[y]), col = transparent_cols[3])
  }
  
  lines(x = years[!is.na(catch_rate_sy)], y = catch_rate_sy[!is.na(catch_rate_sy)])
}





