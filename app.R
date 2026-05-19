library(shiny)
library(ggplot2)
library(EnvStats) # Required for Pareto functions

ui <- fluidPage(
  titlePanel("Welch's T-Test vs. Wilcoxon Rank-Sum Test"),
  sidebarLayout(
    sidebarPanel(
      # Section 1: Inputs
      h4("Configuration"),
      selectInput("dist_type", "Distribution Family:", 
                  choices = c("Normal", "Pareto", "Normal (Floor at 0)")),
      
      # --- NORMAL SLIDERS ---
      conditionalPanel(
        condition = "input.dist_type == 'Normal'",
        hr(),
        tags$b("Population 1 (Data)"),
        sliderInput("mean1_n", "Mean:", min = -2, max = 5, value = 2, step = 0.1),
        sliderInput("sd1_n", "SD:", min = 0.5, max = 5, value = 1, step = 0.1),
        hr(),
        tags$b("Population 2 (Control)"),
        sliderInput("mean2_n", "Mean:", min = -2, max = 5, value = 0, step = 0.1),
        sliderInput("sd2_n", "SD:", min = 0.5, max = 5, value = 1, step = 0.1)
      ),
      
      # --- NORMAL (FLOOR) SLIDERS ---
      conditionalPanel(
        condition = "input.dist_type == 'Normal (Floor at 0)'",
        hr(),
        tags$b("Population 1"),
        sliderInput("mean1_f", "Mean (Before Floor):", min = -2, max = 5, value = 2, step = 0.1),
        sliderInput("sd1_f", "SD:", min = 0.5, max = 5, value = 1, step = 0.1),
        hr(),
        tags$b("Population 2"),
        sliderInput("mean2_f", "Mean (Before Floor):", min = -2, max = 5, value = 0, step = 0.1),
        sliderInput("sd2_f", "SD:", min = 0.5, max = 5, value = 1, step = 0.1)
      ),
      
      # --- PARETO SLIDERS ---
      conditionalPanel(
        condition = "input.dist_type == 'Pareto'",
        hr(),
        tags$b("Population 1 (Data)"),
        sliderInput("alpha1_p", "Shape (Alpha):", min = 1.1, max = 5, value = 3, step = 0.1),
        sliderInput("scale1_p", "Scale (xm):", min = 1, max = 5, value = 1, step = 0.1),
        hr(),
        tags$b("Population 2 (Control)"),
        sliderInput("alpha2_p", "Shape (Alpha):", min = 1.1, max = 5, value = 3, step = 0.1),
        sliderInput("scale2_p", "Scale (xm):", min = 1, max = 5, value = 1, step = 0.1)
      ),
      
      hr(),
      sliderInput("n", "Sample Size (n per group):", min = 2, max = 500, value = 30),
      actionButton("resample", "Run Simulation", class = "btn-primary", width = "100%"),
      
      # Section 2: Statistical Guide
      hr(),
      h4("Statistical Guide"),
      wellPanel(
        style = "background-color: #f8f9fa; border: 1px solid #ddd; padding: 10px;",
        tags$small(
          conditionalPanel(
            condition = "input.dist_type == 'Normal'",
            tags$p(tags$b("Normal:"), " The 'Bell Curve.' Symmetric and predictable. T-tests are theoretically optimized for this shape."),
            hr()
          ),
          conditionalPanel(
            condition = "input.dist_type == 'Pareto'",
            tags$p(tags$b("Pareto:"), " Highly skewed with a 'heavy' right tail. Often causes false alarms (Type I error) because outliers inflate the mean and variance disproportionately."),
            hr()
          ),
          conditionalPanel(
            condition = "input.dist_type == 'Normal (Floor at 0)'",
            tags$p(tags$b("Normal (Floor):"), " Simulates data that cannot be negative (e.g., income). The 'spike' at zero introduces skewness that can bias the T-statistic."),
            hr()
          ),
          tags$p(tags$b("Power:"), " Probability of correctly catching a real effect (Rejecting H0 when the groups are actually different). Goal: ≥ 0.80."),
          tags$p(tags$b("Type I Error:"), " Probability of a 'false alarm' (Rejecting H0 when the groups are identical). Goal: ≤ 0.06. If red, the test is not robust for this distribution.")
        )
      )
    ),
    
    mainPanel(
      # --- NEW DESCRIPTION BOX ---
      wellPanel(
        style = "background-color: #ffffff; border: 1px solid #ddd;",
        tags$p(tags$b("Description: "), "This R-Shiny App is designed to visualize instances in which the normality assumption is violated in three different types of distributions (Normal, Normal with floor effect, and Pareto) and how Welch’s T-test and Wilcoxon Rank-Sum Test perform under these conditions. To use, select the desired distribution and choose appropriate parameters for both populations. Results will be displayed below, additional information is found in the Statistical Guide.")
      ),
      
      # --- EXISTING PLOTS AND RESULTS ---
      fluidRow(
        column(12, plotOutput("popPlot", height = "260px")),
        column(12, plotOutput("tPlot", height = "260px"))
      ),
      fluidRow(
        column(12, 
               wellPanel(
                 fluidRow(
                   column(5, h4("Current Sample Results"), htmlOutput("results")),
                   column(7, h4("Simulation (1,000 runs)"), htmlOutput("sim_results"))
                 )
               )
        )
      )
    )
  )
)

server <- function(input, output) {
  
  test_data <- eventReactive(input$resample, {
    n <- input$n
    dist <- input$dist_type
    
    if(dist == "Normal") { 
      p1_a <- input$mean1_n; p1_b <- input$sd1_n
      p2_a <- input$mean2_n; p2_b <- input$sd2_n
    } else if(dist == "Normal (Floor at 0)") { 
      p1_a <- input$mean1_f; p1_b <- input$sd1_f
      p2_a <- input$mean2_f; p2_b <- input$sd2_f
    } else { 
      p1_a <- input$alpha1_p; p1_b <- input$scale1_p
      p2_a <- input$alpha2_p; p2_b <- input$scale2_p
    }
    
    get_samp <- function(N, D, P_A, P_B) {
      if(D == "Normal") return(rnorm(N, P_A, P_B))
      if(D == "Pareto") return(rpareto(N, location = P_B, shape = P_A))
      if(D == "Normal (Floor at 0)") return(pmax(0, rnorm(N, P_A, P_B))) 
    }
    
    samp1 <- get_samp(n, dist, p1_a, p1_b)
    samp2 <- get_samp(n, dist, p2_a, p2_b) 
    
    t_res <- t.test(samp1, samp2)
    w_res <- wilcox.test(samp1, samp2, exact = FALSE)
    
    sim_results <- replicate(1000, {
      # Power Sim
      s_pow_1 <- get_samp(n, dist, p1_a, p1_b)
      s_pow_2 <- get_samp(n, dist, p2_a, p2_b)
      
      p_t_pow <- t.test(s_pow_1, s_pow_2)$p.value < 0.05
      p_w_pow <- wilcox.test(s_pow_1, s_pow_2, exact = FALSE)$p.value < 0.05 
      
      # Type I Error Sim
      s_t1_1 <- get_samp(n, dist, p1_a, p1_b)
      s_t1_2 <- get_samp(n, dist, p1_a, p1_b)
      
      p_t_t1 <- t.test(s_t1_1, s_t1_2)$p.value < 0.05
      p_w_t1 <- wilcox.test(s_t1_1, s_t1_2, exact = FALSE)$p.value < 0.05
      
      c(p_t_pow, p_w_pow, p_t_t1, p_w_t1)
    })
    
    list(samp1 = samp1, samp2 = samp2, 
         t_stat = t_res$statistic, df = t_res$parameter,
         p_t = t_res$p.value, p_w = w_res$p.value,
         t_pow = mean(sim_results[1,]), w_pow = mean(sim_results[2,]),
         t_t1 = mean(sim_results[3,]), w_t1 = mean(sim_results[4,]))
  }, ignoreNULL = FALSE)
  
  output$popPlot <- renderPlot({
    d <- test_data()
    df_plot <- data.frame(
      val = c(d$samp1, d$samp2),
      Group = factor(rep(c("Sample 1", "Sample 2"), each = input$n))
    )
    
    ggplot(df_plot, aes(x = val, fill = Group)) +
      geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.5, bins = 30, color = "white") +
      geom_vline(xintercept = mean(d$samp1), color = "steelblue", linewidth = 1, linetype = "dashed") +
      geom_vline(xintercept = mean(d$samp2), color = "indianred", linewidth = 1, linetype = "dashed") +
      scale_fill_manual(values = c("Sample 1 (Data)" = "steelblue", "Sample 2 (Control)" = "indianred")) +
      theme_minimal() + 
      theme(legend.position = "top") +
      labs(title = "Two-Sample Distribution Comparison", subtitle = "Dashed lines = Group Sample Means", x = "Value", y = "Density")
  })
  
  output$tPlot <- renderPlot({
    d <- test_data()
    df_val <- d$df 
    t_crit <- qt(0.975, df_val)
    x_ax <- seq(-6, 6, length.out = 200)
    
    ggplot(data.frame(x = x_ax, y = dt(x_ax, df_val)), aes(x, y)) +
      geom_line() +
      geom_area(data = subset(data.frame(x=x_ax, y=dt(x_ax, df_val)), x > t_crit), fill = "red", alpha = 0.3) +
      geom_area(data = subset(data.frame(x=x_ax, y=dt(x_ax, df_val)), x < -t_crit), fill = "red", alpha = 0.3) +
      geom_vline(xintercept = d$t_stat, linewidth = 1.2, color = "black") + 
      theme_minimal() + 
      labs(title = paste0("T-Distribution (Welch df = ", round(df_val, 1), ")"), subtitle = "Solid line = Observed T-statistic", x = "t", y = "Density")
  })
  
  output$results <- renderUI({
    d <- test_data()
    HTML(paste0("<div style='font-size: 1.1em;'>",
                "<b>Mean 1 (Data):</b> ", round(mean(d$samp1), 4), "<br>",
                "<b>Mean 2 (Control):</b> ", round(mean(d$samp2), 4), "<br>",
                "<b>T-Statistic:</b> ", round(d$t_stat, 4), "<br>",
                "<hr style='margin: 8px 0;'>",
                "<b>Welch P-Value:</b> ", round(d$p_t, 4), "<br>",
                "<b>Decision:</b> ", ifelse(d$p_t < 0.05, "Reject Null", "Fail to Reject Null"), "<br><br>",
                "<b>Wilcoxon P-Value:</b> ", round(d$p_w, 4), "<br>",
                "<b>Decision:</b> ", ifelse(d$p_w < 0.05, "Reject Null", "Fail to Reject Null"),
                "</div>"))
  })
  
  output$sim_results <- renderUI({
    d <- test_data()
    
    get_col <- function(val, type) {
      if(type == "pow") return(if(val >= 0.8) "green" else "red")
      return(if(val <= 0.06) "green" else "red")
    }
    
    HTML(paste0(
      "<div style='font-size: 1.1em;'>",
      "<u><b>Welch's T-Test</b></u><br>",
      "<span style='color:", get_col(d$t_pow, "pow"), "'><b>Power: ", round(d$t_pow, 3), "</b> (Goal ≥ 0.8)</span><br>",
      "<span style='color:", get_col(d$t_t1, "t1"), "'><b>Type I Error: ", round(d$t_t1, 3), "</b> (Goal ≤ 0.06)</span><br><br>",
      "<u><b>Wilcoxon Rank-Sum Test</b></u><br>",
      "<span style='color:", get_col(d$w_pow, "pow"), "'><b>Power: ", round(d$w_pow, 3), "</b> (Goal ≥ 0.8)</span><br>",
      "<span style='color:", get_col(d$w_t1, "t1"), "'><b>Type I Error: ", round(d$w_t1, 3), "</b> (Goal ≤ 0.06)</span>",
      "</div>"
    ))
  })
}

shinyApp(ui, server)

