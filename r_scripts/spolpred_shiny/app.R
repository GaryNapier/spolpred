# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

data_url <- "https://raw.githubusercontent.com/GaryNapier/spolpred/master/results/spol_lin_levels_table_github.csv"

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Spoligotype and SNP-based lineage associations"),
  
  # mainPanel(
  tags$br(),
  img(src="LSHTM-logo-bw.jpeg", width='400px'),
  tags$br(),
  tags$br(),
  img(src="tb-profiler-logo.png", width='400px'),
  
  h4(tags$div(
    "Data source:",
    data_url,
    tags$br(), tags$br()
  )),
  
  h3(tags$div(
    "Data description:",
    tags$br()
  )),
  
  h4(tags$div(
    "This is a table of associations between spoligotypes and SNP-based lineages defined in Napier (2020).",
    tags$br(),
    "Lineages are in a numerical hierarchy - 1, 1.1, 1.1.1 etc. - which have been given a 'Level' - 1, 2, 3 etc",
    tags$br(),
    "The important column is 'Correlation' which measures how closely associated a given spoligotype is to a lineage, at that level.",
    tags$br(),
    "A score of 1 means the spoligotype is uniquely found in the lineage-level.",
    tags$br(),
    tags$br(),
    "Example:",
    tags$br(),
    tags$br(),
    tableOutput("example"),
    tags$br(),
    "Spoligotype 110...111 (SIT 19) is correlated 1 with lineage 1 (level 1). i.e. this spoligotype is only found
    in this lineage.",
    tags$br(),
    tags$br(),
    "Columns:",
    tags$br(),
    tags$br(),
    "- Level: Level of the lineage hierarchy (e.g. 1 = level 1, 1.1 = level 2, 1.1.1 = level 3 etc.)",
    tags$br(),
    tags$br(),
    "- Lineage: Lineages in Napier (2020).",
    tags$br(),
    tags$br(),
    "- Spoligotype: Spoligotype predicted by TB-Profiler (Phelan 2019).",
    tags$br(),
    tags$br(),
    "- SIT: SIT number from SITVIT2; database of genotyping markers for Mtb.",
    tags$br(),
    tags$br(),
    "- Family: Spoligotype family - see Brudley (2006).",
    tags$br(),
    tags$br(),
    "- Correlation: Scale 0-1 of how associated the spoligotype is to the lineage-level. 1 = spoligotype is uniquely found in the lineage-level. Comparable to Fst score.",
    tags$br(),
    tags$br(),
    "- N: Number of samples with this spoligotype belonging to the lineage-level.",
    tags$br(),
    tags$br(),
    "- % in level-lineage: Percentage of samples with this spoligotype at the level-lineage.",
    tags$br()
  )),
  
  h4(tags$div(
    tags$br(),
    tags$br(),
    "Comprehensive MTBC tree from Napier 2020"
  )),
  
  tags$br(),
  img(src="mtbc_tree.png", width='1000px'),
  tags$br(),
  
  h5(tags$div(
    "References: ",
    tags$br(),
    tags$br(),
    "Napier, G., Campino, S., Merid, Y. et al. Robust barcoding and identification of Mycobacterium tuberculosis lineages for epidemiological and clinical studies. Genome Med 12, 114 (2020). https://doi.org/10.1186/s13073-020-00817-3",
    tags$br(),
    tags$br(),
    "Phelan, J. E. et al. Integrating informatics tools and portable sequencing technology for rapid detection of resistance to anti-tuberculous drugs. Genome Med. 11, 41 (2019).",
    tags$br(),
    tags$br(),
    "SITVIT2: http://www.pasteur-guadeloupe.fr:8081/SITVIT2/",
    tags$br(),
    tags$br(),
    "Brudey, K., Driscoll, J.R., Rigouts, L. et al. Mycobacterium tuberculosis complex genetic diversity: mining the fourth international spoligotyping database (SpolDB4) for classification, population genetics and epidemiology. BMC Microbiol 6, 23 (2006). https://doi.org/10.1186/1471-2180-6-23",
    tags$br(),
    tags$br()
  )),
  
  sidebarLayout(
    
    sidebarPanel(
      helpText("Filter results table (enter comma-separated values for multiple)"),
      textInput(inputId = "level", 
                label = "Level", 
                value = NULL),
      textInput(inputId = "lineage",
                label = "Lineage",
                value = NULL),
      textInput(inputId = "spoligotype",
                label = "Spoligotype",
                value = NULL),
      textInput(inputId = "SIT",
                label = "SIT",
                value = NULL),
      textInput(inputId = "family",
                label = "Family",
                value = NULL),
      submitButton("Submit")
    ), # sidebarPanel
    
    mainPanel(
      tableOutput("spol_out"),
      downloadButton('downloadData', 'Download')
    )
  ) # sidebarLayout
)

server <- function(input, output) {
  
  # data_url <- "https://raw.githubusercontent.com/GaryNapier/spolpred/master/results/spol_lin_levels_table_github.csv"
  spol_lin_levels_table <- read.csv(data_url, colClasses = c("Spoligotype" = "character"))
  
  names(spol_lin_levels_table) <- c("Level",
                                    "Lineage",
                                    "Spoligotype",
                                    "SIT",
                                    "Family",
                                    "Correlation",
                                    "N",
                                    "% in level-lineage")
  
  output$example <- renderTable({
    df <- data.frame(
      Level = "1",
      Lineage = "1",
      Spoligotype = "1101111111111111111001111111000010111111111",
      SIT = "19",
      Family = "EAI2-Manila",
      Correlation = 1.00,
      N = "402",
      pc = 17.20
    )
    df <- dplyr::rename(df, "% in level-lineage" = pc)
    df
  })
  
  # Subset the table for the main results

  level_values <- reactive({
    if (isTruthy(input$level)) {
      return(input$level)
    } else {
      return(unique(spol_lin_levels_table$Level))
    }
  })
  
  lineage_values <- reactive({
    if (isTruthy(input$lineage)) {
      return(input$lineage)
    } else {
      return(unique(spol_lin_levels_table$Lineage))
    }
  })
  
  spol_values <- reactive({
    if (isTruthy(input$spoligotype)) {
      return(input$spoligotype)
    } else {
      return(unique(spol_lin_levels_table$Spoligotype))
    }
  })
  
  sit_values <- reactive({
    if (isTruthy(input$SIT)) {
      return(input$SIT)
    } else {
      return(unique(spol_lin_levels_table$SIT))
    }
  })
  
  family_values <- reactive({
    if (isTruthy(input$family)) {
      return(input$family)
    } else {
      return(unique(spol_lin_levels_table$Family))
    }
  })
  
  
  filter_df <- reactive({
    return(spol_lin_levels_table %>% 
             filter(Level %in% level_values(),
                    Lineage %in% lineage_values(),
                    Spoligotype %in% spol_values(),
                    SIT %in% sit_values(),
                    Family %in% family_values()
                    ) # filter
    ) # return
  }) # reactive
  
  output$spol_out <- renderTable({filter_df()}) # renderTable
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("spol_correlation_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(filter_df(), file, row.names = F) # call to the reactive again
    }
  )
  
} # server()

# Run the application
shinyApp(ui = ui, server = server)



