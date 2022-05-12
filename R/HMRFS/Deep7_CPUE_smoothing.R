# Before running, you must first load "deep7.RData" into the current environment

library(shiny)
library(zoo)
library(KFAS)

cols = c("black", "deeppink", "deepskyblue", "blueviolet", "lawngreen")

ui = fluidPage(
  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #000000;}")), # making the hr() lines black
  ),
  titlePanel("Deep 7 HMRFS CPUE Smoothing"),
  sidebarLayout(
    sidebarPanel(
      selectInput("species", "Species", choices = list("Hapu'upu'u" = 1,
                                                       "Onaga" = 2,
                                                       "Ehu" = 3,
                                                       "Opakapaka" = 4,
                                                       "Kalekale" = 5,
                                                       "Lehi" = 6,
                                                       "Gindai" = 7), selected = 3),
      checkboxGroupInput("lines", "Series to Plot", choices = list("Moving Average" = 1,
                                                                   "Moving Average by Area" = 2,
                                                                   "Kalman" = 3,
                                                                   "Kalman by Area" = 4,
                                                                   "Official Values" = 5), selected = 1:5),
      hr(),
      h4("Annual Smoothing Inputs"),
      numericInput("ma_y", "Moving Average # of Years +/-", value = 1, step = 1),
      numericInput("q_y", "Kalman Q", value = 0.01), # Random walk for annual smoothing
      numericInput("h_y", "Kalman H", value = 0.01), # Observation process for annual smoothing
      hr(),
      h4("Wave Smoothing Inputs"),
      numericInput("ma_w", "Moving Average # of Waves +/-", value = 6, step = 1),
      numericInput("q_w", "Kalman Q", value = 0.01), # Random walk for wave smoothing
      numericInput("h_w", "Kalman H", value = 0.01), # Observation process for wave smoothing

    ),
    mainPanel(
      plotOutput("plot_y"),
      plotOutput("plot_w")
    )
  )
)

server = function(input, output) {
  rvalues = reactiveValues()
  rvalues$year_smoothed_total_catch_ma = array(0, dim = c(n_years, 3))
  rvalues$year_smoothed_total_catch_kalman = array(0, dim = c(n_years, 3))
  rvalues$wave_smoothed_total_catch_ma = array(0, dim = c(n_years, 3))
  rvalues$wave_smoothed_total_catch_kalman = array(0, dim = c(n_years, 3))
  
  observeEvent(input$species, {
    s = as.numeric(input$species)
    
    # annual smoothing
    
    total_caught = observed_total_caught[, , s, , , ] + unavailable_total_caught[, , s, , , ]
    
    rvalues$year_smoothed_total_catch_ma[, 3] = rollapplyr(apply(total_caught[, , 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum), 1 + 2 * input$ma_y, FUN = mean, partial = T, align = "center") * apply(effort[, , 1, 1:2], c(1), sum)
    
    for(a in 1:(n_areas - 1)) {
      rvalues$year_smoothed_total_catch_ma[, a] = rollapplyr(apply(total_caught[, , 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum), 1 + 2 * input$ma_y, FUN = mean, partial = T, align = "center") * apply(effort[, , 1, a], c(1), sum)
    }
    
    cpues = apply(total_caught[, , 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum)
    rvalues$year_smoothed_total_catch_kalman[, 3] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_y), H = input$h_y))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, 1:2]))
    
    for(a in 1:(n_areas - 1)) {
      cpues = apply(total_caught[, , 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum)
      rvalues$year_smoothed_total_catch_kalman[, a] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_y), H = input$h_y))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, a]))
    }
    
    # wave smoothing
    
    total_caught_yws = total_caught_yw[, s, , , ]
    wave_smoothed_total_catch_ma = array(0, dim = c(n_year_waves, 3))
    wave_smoothed_total_catch_kalman = array(0, dim = c(n_year_waves, 3))
    
    wave_smoothed_total_catch_ma[, 3] = rollapplyr(apply(total_caught_yws[, 1, 1:2, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, 1:2, ], c(1), sum), 1 + 2 * input$ma_w, FUN = mean, partial = T, align = "center", na.rm = T) * apply(effort_yw[, 1, 1:2], c(1), sum)
    
    for(a in 1:(n_areas - 1)) {
      wave_smoothed_total_catch_ma[, a] = rollapplyr(apply(total_caught_yws[, 1, a, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, a, ], c(1), sum), 1 + 2 * input$ma_w, FUN = mean, partial = T, align = "center", na.rm = T) * effort_yw[, 1, a]
    }
    
    cpues = apply(total_caught_yws[, 1, 1:2, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, 1:2, ], c(1), sum)
    wave_smoothed_total_catch_kalman[, 3] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_w), H = input$h_w))$alphahat * sapply(1:n_year_waves, function(f) sum(effort_yw[f, 1, 1:2]))
    
    for(a in 1:(n_areas - 1)) {
      cpues = apply(total_caught_yws[, 1, a, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, a, ], c(1), sum)
      wave_smoothed_total_catch_kalman[, a] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_w), H = input$h_w))$alphahat * sapply(1:n_year_waves, function(f) sum(effort_yw[f, 1, a]))
    }
    
    for(area in 1:3) {
      rvalues$wave_smoothed_total_catch_ma[, area] = sapply(1:n_years, function(f) sum(wave_smoothed_total_catch_ma[((f - 1) * 6 + 1):((f - 1) * 6 + 6), area], na.rm = T))
      rvalues$wave_smoothed_total_catch_kalman[, area] = sapply(1:n_years, function(f) sum(wave_smoothed_total_catch_kalman[((f - 1) * 6 + 1):((f - 1) * 6 + 6), area], na.rm = T))
    }
  })
  
  observeEvent(input$ma_y, {
    s = as.numeric(input$species)
    
    # annual smoothing
    
    total_caught = observed_total_caught[, , s, , , ] + unavailable_total_caught[, , s, , , ]

    rvalues$year_smoothed_total_catch_ma[, 3] = rollapplyr(apply(total_caught[, , 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum), 1 + 2 * input$ma_y, FUN = mean, partial = T, align = "center") * apply(effort[, , 1, 1:2], c(1), sum)
    
    for(a in 1:(n_areas - 1)) {
      rvalues$year_smoothed_total_catch_ma[, a] = rollapplyr(apply(total_caught[, , 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum), 1 + 2 * input$ma_y, FUN = mean, partial = T, align = "center") * apply(effort[, , 1, a], c(1), sum)
    }
  })
  
  observeEvent(input$ma_w, {
    s = as.numeric(input$species)
    
    # wave smoothing
    
    total_caught_yws = total_caught_yw[, s, , , ]
    wave_smoothed_total_catch_ma = array(0, dim = c(n_year_waves, 3))

    wave_smoothed_total_catch_ma[, 3] = rollapplyr(apply(total_caught_yws[, 1, 1:2, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, 1:2, ], c(1), sum), 1 + 2 * input$ma_w, FUN = mean, partial = T, align = "center", na.rm = T) * apply(effort_yw[, 1, 1:2], c(1), sum)
    
    for(a in 1:(n_areas - 1)) {
      wave_smoothed_total_catch_ma[, a] = rollapplyr(apply(total_caught_yws[, 1, a, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, a, ], c(1), sum), 1 + 2 * input$ma_w, FUN = mean, partial = T, align = "center", na.rm = T) * effort_yw[, 1, a]
    }
    
    for(area in 1:3) {
      rvalues$wave_smoothed_total_catch_ma[, area] = sapply(1:n_years, function(f) sum(wave_smoothed_total_catch_ma[((f - 1) * 6 + 1):((f - 1) * 6 + 6), area], na.rm = T))
    }
  })
  
  observeEvent(input$q_y, {
    s = as.numeric(input$species)
    
    total_caught = observed_total_caught[, , s, , , ] + unavailable_total_caught[, , s, , , ]
    
    cpues = apply(total_caught[, , 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum)
    rvalues$year_smoothed_total_catch_kalman[, 3] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_y), H = input$h_y))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, 1:2]))
    
    for(a in 1:(n_areas - 1)) {
      cpues = apply(total_caught[, , 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum)
      rvalues$year_smoothed_total_catch_kalman[, a] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_y), H = input$h_y))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, a]))
    }
  })
    
  observeEvent(input$q_w, {
    s = as.numeric(input$species)
    
    total_caught_yws = total_caught_yw[, s, , , ]
    wave_smoothed_total_catch_kalman = array(0, dim = c(n_year_waves, 3))
    
    cpues = apply(total_caught_yws[, 1, 1:2, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, 1:2, ], c(1), sum)
    wave_smoothed_total_catch_kalman[, 3] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_w), H = input$h_w))$alphahat * sapply(1:n_year_waves, function(f) sum(effort_yw[f, 1, 1:2]))
    
    for(a in 1:(n_areas - 1)) {
      cpues = apply(total_caught_yws[, 1, a, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, a, ], c(1), sum)
      wave_smoothed_total_catch_kalman[, a] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_w), H = input$h_w))$alphahat * sapply(1:n_year_waves, function(f) sum(effort_yw[f, 1, a]))
    }
    
    for(area in 1:3) {
      rvalues$wave_smoothed_total_catch_kalman[, area] = sapply(1:n_years, function(f) sum(wave_smoothed_total_catch_kalman[((f - 1) * 6 + 1):((f - 1) * 6 + 6), area], na.rm = T))
    }
  })
  
  observeEvent(input$h_y, {
    s = as.numeric(input$species)
    
    total_caught = observed_total_caught[, , s, , , ] + unavailable_total_caught[, , s, , , ]
    
    cpues = apply(total_caught[, , 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum)
    rvalues$year_smoothed_total_catch_kalman[, 3] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_y), H = input$h_y))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, 1:2]))
    
    for(a in 1:(n_areas - 1)) {
      cpues = apply(total_caught[, , 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum)
      rvalues$year_smoothed_total_catch_kalman[, a] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_y), H = input$h_y))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, a]))
    }
  })
  
  observeEvent(input$h_w, {
    s = as.numeric(input$species)
    
    total_caught_yws = total_caught_yw[, s, , , ]
    wave_smoothed_total_catch_kalman = array(0, dim = c(n_year_waves, 3))
    
    cpues = apply(total_caught_yws[, 1, 1:2, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, 1:2, ], c(1), sum)
    wave_smoothed_total_catch_kalman[, 3] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_w), H = input$h_w))$alphahat * sapply(1:n_year_waves, function(f) sum(effort_yw[f, 1, 1:2]))
    
    for(a in 1:(n_areas - 1)) {
      cpues = apply(total_caught_yws[, 1, a, ], c(1), sum) / apply(anglers_by_trip_yw[, 1, a, ], c(1), sum)
      wave_smoothed_total_catch_kalman[, a] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q_w), H = input$h_w))$alphahat * sapply(1:n_year_waves, function(f) sum(effort_yw[f, 1, a]))
    }
    
    for(area in 1:3) {
      rvalues$wave_smoothed_total_catch_kalman[, area] = sapply(1:n_years, function(f) sum(wave_smoothed_total_catch_kalman[((f - 1) * 6 + 1):((f - 1) * 6 + 6), area], na.rm = T))
    }
  })
  
  output$plot_y = renderPlot({
    s = as.numeric(input$species)
    
    unaccounted_total_catch = apply(total_catch_num[, , s, 2,], c(1), sum, na.rm = T)
    ma_smoothed_catch_plot = rvalues$year_smoothed_total_catch_ma[, 3] + unaccounted_total_catch
    ma_smoothed_by_area_catch_plot = apply(rvalues$year_smoothed_total_catch_ma[, 1:2], c(1), sum) + unaccounted_total_catch
    kalman_smoothed_catch_plot = rvalues$year_smoothed_total_catch_kalman[, 3] + unaccounted_total_catch
    kalman_smoothed_by_area_catch_plot = apply(rvalues$year_smoothed_total_catch_kalman[, 1:2], c(1), sum) + unaccounted_total_catch
    
    par(xpd = NA, mar = c(4, 5, 3, 1))
    
    x_min = min(years)
    x_max = max(years)
    y_min = 0
    y_ticks = pretty(c(0, official_catch_num[, s]), n = 4)
    y_max = max(y_ticks)

    plot(NA, axes = F, xaxs = "i", yaxs = "i", xlim = c(x_min, x_max), ylim = c(y_min, y_max), xlab = "", ylab = "")
    axis(side = 1, at = seq(from = 2005, to = 2020, by = 5))
    lines(x = c(x_min, x_max), y = c(y_min, y_min))
    axis(side = 2, at = y_ticks[y_ticks >= 0], las = 2)
    lines(x = c(x_min, x_min), y = c(y_min, y_max))
    mtext("Year", side = 1, line = 3)
    mtext("Catch Number per Year", side = 2, line = 4)
    mtext(paste0(species_df$hawaiian_name[s], " (Year-Smoothed)"), side = 3, line = 1)
    
    selected = c()
    if(1 %in% input$lines) {
      lines(years, ma_smoothed_catch_plot, col = cols[2])
      selected = c(1)
    }
    if(2 %in% input$lines) {
      lines(years, ma_smoothed_by_area_catch_plot, col = cols[2], lty = "dashed")
      selected = c(selected, 2)
    }
    if(3 %in% input$lines) {
      lines(years, kalman_smoothed_catch_plot, col = cols[3])
      selected = c(selected, 3)
    }
    if(4 %in% input$lines) {
      lines(years, kalman_smoothed_by_area_catch_plot, col = cols[3], lty = "dashed")
      selected = c(selected, 4)
    }
    if(5 %in% input$lines) {
      lines(years, official_catch_num[, s], col = cols[1])
      selected = c(selected, 5)
    }
    
    if(length(input$lines) > 0) {
      legend(x = "topright", col = c(cols[2], cols[2], cols[3], cols[3], cols[1])[selected], lty = c("solid", "dashed", "solid", "dashed", "solid")[selected], legend = c("Moving Average", "Moving Average by area", "Kalman", "Kalman by area", "Original")[selected], bty = "n")
    }
  })
  
  
  
  output$plot_w = renderPlot({
    s = as.numeric(input$species)
    
    unaccounted_total_catch = apply(total_catch_num[, , s, 2,], c(1), sum, na.rm = T)
    ma_smoothed_catch_plot = rvalues$wave_smoothed_total_catch_ma[, 3] + unaccounted_total_catch
    ma_smoothed_by_area_catch_plot = apply(rvalues$wave_smoothed_total_catch_ma[, 1:2], c(1), sum) + unaccounted_total_catch
    kalman_smoothed_catch_plot = rvalues$wave_smoothed_total_catch_kalman[, 3] + unaccounted_total_catch
    kalman_smoothed_by_area_catch_plot = apply(rvalues$wave_smoothed_total_catch_kalman[, 1:2], c(1), sum) + unaccounted_total_catch
    
    par(xpd = NA, mar = c(4, 5, 3, 1))
    
    x_min = min(years)
    x_max = max(years)
    y_min = 0
    y_ticks = pretty(c(0, official_catch_num[, s]), n = 4)
    y_max = max(y_ticks)
    
    plot(NA, axes = F, xaxs = "i", yaxs = "i", xlim = c(x_min, x_max), ylim = c(y_min, y_max), xlab = "", ylab = "")
    axis(side = 1, at = seq(from = 2005, to = 2020, by = 5))
    lines(x = c(x_min, x_max), y = c(y_min, y_min))
    axis(side = 2, at = y_ticks[y_ticks >= 0], las = 2)
    lines(x = c(x_min, x_min), y = c(y_min, y_max))
    mtext("Year", side = 1, line = 3)
    mtext("Catch Number per Year", side = 2, line = 4)
    mtext(paste0(species_df$hawaiian_name[s], " (Wave-Smoothed)"), side = 3, line = 1)
    
    selected = c()
    if(1 %in% input$lines) {
      lines(years, ma_smoothed_catch_plot, col = cols[2])
      selected = c(1)
    }
    if(2 %in% input$lines) {
      lines(years, ma_smoothed_by_area_catch_plot, col = cols[2], lty = "dashed")
      selected = c(selected, 2)
    }
    if(3 %in% input$lines) {
      lines(years, kalman_smoothed_catch_plot, col = cols[3])
      selected = c(selected, 3)
    }
    if(4 %in% input$lines) {
      lines(years, kalman_smoothed_by_area_catch_plot, col = cols[3], lty = "dashed")
      selected = c(selected, 4)
    }
    if(5 %in% input$lines) {
      lines(years, official_catch_num[, s], col = cols[1])
      selected = c(selected, 5)
    }
    
    if(length(input$lines) > 0) {
      legend(x = "topright", col = c(cols[2], cols[2], cols[3], cols[3], cols[1])[selected], lty = c("solid", "dashed", "solid", "dashed", "solid")[selected], legend = c("Moving Average", "Moving Average by area", "Kalman", "Kalman by area", "Original")[selected], bty = "n")
    }
  })
}

shinyApp(ui = ui, server = server)