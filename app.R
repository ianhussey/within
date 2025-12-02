
# todo: 
# pivot p values in output
# flag input
# add consistency checking among d_z, t, and p values
# add notes that if p is approximated it will provide a limit on r but not a precise value.
# what's happening when the table doesn't display with extreme values? is r outside [-1, 1]? then display a message that r is impossible, but i think its more complicated than this.

# references for plausible r being probably [.5, .75] or more generously [.25, .90] but unlikely to be beyond this.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6331475/
# https://www.researchgate.net/publication/247720873_Estimating_Effect_Sizes_From_Pretest-Posttest-Control_Group_Designs
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6998624/
# https://www.ncbi.nlm.nih.gov/books/NBK154408/
# https://www.ncbi.nlm.nih.gov/books/NBK115798/
# https://www.campbellcollaboration.org/calculator/equations


# libraries ----
suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinydashboard)
  library(shinyWidgets)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(DT)
  library(faux)
  library(scales)
  library(cowplot)
})

# setup ----


# functions ----
source("scripts/func.R") # helper functions

# user interface ----

## UI ----
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Within-subject t-test forensics",
                  titleWidth = "calc(100% - 44px)" # puts sidebar toggle on right
  ),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    shinyjs::useShinyjs(),
    tags$head(
      # links to files in www/
      tags$link(rel = "stylesheet", type = "text/css", href = "basic_template.css"),
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
      tags$script(src = "custom.js")
    ),
    column(width = 4,
           ## data paramaters ----
           box(id = "input_box",
               title = "Reported Data Parameters", width = 12, solidHeader = TRUE,
               p("Check the possible p-values and t-values for the given N, means and SDs for correlations between -1 and +1."),
               numericInput("n", "N", 20, min = 1, step = .01),
               numericInput("m1", "Mean 1", 0.2, step = .01),
               numericInput("sd1", "SD 1", 1, min = 0.0001, step = .01),
               numericInput("m2", "Mean 2", 0, step = .01),
               numericInput("sd2", "SD 2", 1, min = 0.0001, step = .01)
           ),
           
           ## reported test parameters ----
           box(id = "test_box",
               title = "Reported Test Parameters", width = 12, solidHeader = TRUE,
               p("Given a t-value or p-value, get the possible correlations for one- and two-tailed tests. If both are supplied, ensure they produce consistent results."),
               numericInput("reported_t", "t-value", NULL, step = .01),
               numericInput("reported_p", "p-value", NULL, min = 0, max = 1, step = .01)
           ),
           tags$a(href="https://github.com/ianhussey/within", "Original code by Lisa DeBruine, modifications by Ian Hussey. Code & Citation.")
    ),
    column(width = 8,
           box(title = "t-values", width = 6, solidHeader = TRUE,
               plotOutput("t_plot")
           ),
           box(title = "p-values", width = 6, solidHeader = TRUE,
               plotOutput("p_plot")
           ),
           DTOutput("r_table")
    )
  )
)


# server ----
server <- function(input, output, session) {
  ## r_table ----
  r_table <- reactive({
    check_values(
      n = input$n,
      m1 = input$m1,
      sd1 = input$sd1,
      m2 = input$m2,
      sd2 = input$sd2
    )
  })
  
  ## t_plot ----
  output$t_plot <- renderPlot({
    req(r_table())
    t_plot(r_table(), input$reported_t)
  })
  
  ## p_plot ----
  output$p_plot <- renderPlot({
    req(r_table())
    p_plot(r_table(), input$reported_p)
  })
  
  less_r <- reactive({
    solve_r(param = "p",
            reported = input$reported_p,
            n = input$n,
            m1 = input$m1,
            m2 = input$m2,
            sd1 = input$sd1,
            sd2 = input$sd2,
            alternative = "less")
  })
  
  greater_r <- reactive({
    solve_r(param = "p",
            reported = input$reported_p,
            n = input$n,
            m1 = input$m1,
            m2 = input$m2,
            sd1 = input$sd1,
            sd2 = input$sd2,
            alternative = "greater")
  })
  
  two.sided_r <- reactive({
    solve_r(param = "p",
            reported = input$reported_p,
            n = input$n,
            m1 = input$m1,
            m2 = input$m2,
            sd1 = input$sd1,
            sd2 = input$sd2,
            alternative = "two.sided")
  })
  
  ## p_table ----
  p_table <- reactive({
    if (is.na(input$reported_p)) return(NULL)
    message("making p_table")
    
    r    <- list(less_r(), greater_r(), two.sided_r())
    alt  <- c("less", "greater", "two.sided")
    nosol <- sapply(r, is.null)
    
    # if all solutions are NULL, nothing to show
    if (all(nosol)) return(NULL)
    
    r_vals  <- unlist(r)
    alt_ok  <- alt[!nosol]
    
    df_p <- tidyr::crossing(
      n   = input$n,
      m1  = input$m1,
      m2  = input$m2,
      sd1 = input$sd1,
      sd2 = input$sd2,
      r   = r_vals
    ) %>%
      # match each r solution to its alternative
      dplyr::mutate(
        alternative = alt_ok[seq_along(r)]
      ) %>%
      # within_t must return diff_mean, diff_sd, t, p
      dplyr::bind_cols(purrr::pmap_df(., within_t)) %>%
      # add effect sizes
      dplyr::mutate(
        d_z  = t / sqrt(n),
        d_av = diff_mean / ((sd1 + sd2) / 2)
      ) %>%
      dplyr::select(r, diff_sd, t, p, d_z, d_av, alternative) %>%
      dplyr::mutate_if(is.numeric, round, 3)
    
    return(df_p)
  })
  
  ## t_table ----
  t_table <- reactive({
    if (is.na(input$reported_t)) return(NULL)
    message("making t_table")
    
    r <- solve_r(
      param     = "t",
      reported  = input$reported_t,
      n         = input$n,
      m1        = input$m1,
      m2        = input$m2,
      sd1       = input$sd1,
      sd2       = input$sd2,
      alternative = "two.sided"
    )
    
    if (is.null(r)) return(NULL)
    
    alt <- c("less", "greater", "two.sided")
    
    df_t <- tidyr::crossing(
      n   = input$n,
      m1  = input$m1,
      m2  = input$m2,
      sd1 = input$sd1,
      sd2 = input$sd2,
      r   = r,
      alternative = alt
    ) %>%
      dplyr::bind_cols(purrr::pmap_df(., within_t)) %>%
      dplyr::mutate(
        d_z  = t / sqrt(n),
        d_av = diff_mean / ((sd1 + sd2) / 2)
      ) %>%
      dplyr::select(r, diff_sd, t, p, d_z, d_av, alternative) %>%
      dplyr::mutate_if(is.numeric, round, 3)
    
    return(df_t)
  })
  
  ## r_table ----
  output$r_table <- renderDT({
    message("r_table")
    dplyr::bind_rows(p_table(), t_table()) %>%
      pivot_wider(names_from = "alternative",
                  names_prefix = "p_",
                  values_from = "p") %>%
      select(r, diff_sd, t, d_z, d_av, p_greater, p_two.sided, p_less)
  }, options = list(dom = 't'))
}

shinyApp(ui, server)
