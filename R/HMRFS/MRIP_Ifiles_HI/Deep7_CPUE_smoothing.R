library(shiny)
library(zoo)
library(KFAS)

ui = fluidPage(
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
      numericInput("ma", "Moving Average # of Years +/-", value = 1, step = 1),
      numericInput("q", "Kalman Q", value = 0.01), # Random walk
      numericInput("h", "Kalman H", value = 0.01), # Observation process
      checkboxGroupInput("lines", "Series to Plot", choices = list("Moving Average" = 1,
                                                                   "Moving Average by Area" = 2,
                                                                   "Kalman" = 3,
                                                                   "Kalman by Area" = 4,
                                                                   "Official Values" = 5), selected = 1:5)
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
)

server = function(input, output) {
  rvalues = reactiveValues()
  rvalues$year_smoothed_total_catch_ma = array(0, dim = c(n_years, 3))
  rvalues$year_smoothed_total_catch_kalman = array(0, dim = c(n_years, 3))
  
  observeEvent(input$species, {
    s = as.numeric(input$species)
    
    total_caught = observed_total_caught[, , s, , , ] + unavailable_total_caught[, , s, , , ]
    
    rvalues$year_smoothed_total_catch_ma[, 3] = rollapplyr(apply(total_caught[, , 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum), 1 + 2 * input$ma, FUN = mean, partial = T, align = "center") * apply(effort[, , 1, 1:2], c(1), sum)
    
    for(a in 1:(n_areas - 1)) {
      rvalues$year_smoothed_total_catch_ma[, a] = rollapplyr(apply(total_caught[, , 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum), 1 + 2 * input$ma, FUN = mean, partial = T, align = "center") * apply(effort[, , 1, a], c(1), sum)
    }
    
    cpues = apply(total_caught[, , 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum)
    rvalues$year_smoothed_total_catch_kalman[, 3] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q), H = input$h))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, 1:2]))
    
    for(a in 1:(n_areas - 1)) {
      cpues = apply(total_caught[, , 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum)
      rvalues$year_smoothed_total_catch_kalman[, a] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q), H = input$h))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, a]))
    }
  })
  
  observeEvent(input$ma, {
    s = as.numeric(input$species)
    
    total_caught = observed_total_caught[, , s, , , ] + unavailable_total_caught[, , s, , , ]

    rvalues$year_smoothed_total_catch_ma[, 3] = rollapplyr(apply(total_caught[, , 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum), 1 + 2 * input$ma, FUN = mean, partial = T, align = "center") * apply(effort[, , 1, 1:2], c(1), sum)
    
    for(a in 1:(n_areas - 1)) {
      rvalues$year_smoothed_total_catch_ma[, a] = rollapplyr(apply(total_caught[, , 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum), 1 + 2 * input$ma, FUN = mean, partial = T, align = "center") * apply(effort[, , 1, a], c(1), sum)
    }
  })
  
  observeEvent(input$q, {
    s = as.numeric(input$species)
    
    total_caught = observed_total_caught[, , s, , , ] + unavailable_total_caught[, , s, , , ]
    
    cpues = apply(total_caught[, , 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum)
    rvalues$year_smoothed_total_catch_kalman[, 3] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q), H = input$h))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, 1:2]))
    
    for(a in 1:(n_areas - 1)) {
      cpues = apply(total_caught[, , 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum)
      rvalues$year_smoothed_total_catch_kalman[, a] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q), H = input$h))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, a]))
    }
  })
  
  observeEvent(input$h, {
    s = as.numeric(input$species)
    
    total_caught = observed_total_caught[, , s, , , ] + unavailable_total_caught[, , s, , , ]
    
    cpues = apply(total_caught[, , 1, 1:2, ], c(1), sum) / apply(anglers_by_trip[, , 1, 1:2, ], c(1), sum)
    rvalues$year_smoothed_total_catch_kalman[, 3] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q), H = input$h))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, 1:2]))
    
    for(a in 1:(n_areas - 1)) {
      cpues = apply(total_caught[, , 1, a, ], c(1), sum) / apply(anglers_by_trip[, , 1, a, ], c(1), sum)
      rvalues$year_smoothed_total_catch_kalman[, a] = KFS(SSModel(cpues ~ SSMtrend(1, Q = input$q), H = input$h))$alphahat * sapply(1:n_years, function(f) sum(effort[f, , 1, a]))
    }
  })
  
  output$plot = renderPlot({
    s = as.numeric(input$species)
    
    unaccounted_total_catch = apply(total_catch_num[, , s, 2,], c(1), sum, na.rm = T)
    ma_smoothed_catch_plot = rvalues$year_smoothed_total_catch_ma[, 3] + unaccounted_total_catch
    ma_smoothed_by_area_catch_plot = apply(rvalues$year_smoothed_total_catch_ma[, 1:2], c(1), sum) + unaccounted_total_catch
    kalman_smoothed_catch_plot = rvalues$year_smoothed_total_catch_kalman[, 3] + unaccounted_total_catch
    kalman_smoothed_by_area_catch_plot = apply(rvalues$year_smoothed_total_catch_kalman[, 1:2], c(1), sum) + unaccounted_total_catch
    
    par(xpd = NA, mar = c(4, 5, 3, 1))
    
    cols = c("black", "deeppink", "deepskyblue", "blueviolet", "lawngreen")
    
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
    mtext("Catch Number", side = 2, line = 4)
    mtext(species_df$hawaiian_name[s], side = 3, line = 1)
    
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