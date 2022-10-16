#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Get data

data_url <- "https://raw.githubusercontent.com/GaryNapier/spolpred/master/results/spol_lin_levels_table_github.csv"
spol_lin_levels_table <- read.csv(data_url, colClasses = c("Spoligotype" = "character"))

# Clean
# names(spol_lin_levels_table) <- c("Level",
#                                   "Lineage", 
#                                   "Spoligotype", 
#                                   "SIT", 
#                                   "Family", 
#                                   "Spol-lin-prob", 
#                                   "Freq",
#                                   "Lin-level_pc")

names(spol_lin_levels_table) <- c("Level", 
                                  "Lineage", 
                                  "Spoligotype", 
                                  "SIT", 
                                  "Family", 
                                  "Correlation",
                                  "N",
                                  "% in level-lineage")

# Define UI
ui <- fluidPage(
  
    mainPanel(
      tags$br(),
      img(src="LSHTM-logo-bw.jpeg", width='400px'),
      tags$br(),
      tags$br(),
      img(src="tb-profiler-logo.png", width='400px'),

   # Application title
   titlePanel("Spoligotype and lineage associations"),
   
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
     "Lineages are taken from Napier (2020). They are arranged in a numerical hierarchy (e.g. 1, 1.1, 1.1.1 etc.) - see tree figure below",
     tags$br(), 
     tags$br(), 
     "Columns:", 
     tags$br(), 
     tags$br(), 
     "- Level: Level of the lineage hierarchy (e.g. 1 = level 1, 1.1 = level 2, 1.1.1 = level 3 etc.)",
     tags$br(), 
     "- Lineage: Lineages in Napier (2020).",
     tags$br(), 
     "- Spoligotype: Spoligotype predicted by TB-Profiler (Phelan 2019).",
     tags$br(), 
     "- SIT: SIT number from SITVIT2; database of genotyping markers for Mtb.", 
     tags$br(), 
     "- Family: Spoligotype family - see Brudley (2006).", 
     tags$br(), 
     "- Correlation: Scale 0-1 of how associated the spoligotype is to the lineage-level. 1 = spoligotype is uniquely found in the lineage-level. Comparable to Fst score.", 
     tags$br(), 
     "- N: Number of samples with this spoligotype belonging to the lineage-level.",
     tags$br(), 
     "- % in level-lineage: Percentage of samples with this spoligotype at the level-lineage.", 
     tags$br()
   )),
   
   h5(tags$div(
     "References: ", 
     tags$br(), 
     "Napier, G., Campino, S., Merid, Y. et al. Robust barcoding and identification of Mycobacterium tuberculosis lineages for epidemiological and clinical studies. Genome Med 12, 114 (2020). https://doi.org/10.1186/s13073-020-00817-3", 
     tags$br(), 
     "Phelan, J. E. et al. Integrating informatics tools and portable sequencing technology for rapid detection of resistance to anti-tuberculous drugs. Genome Med. 11, 41 (2019).", 
     tags$br(), 
     "SITVIT2: http://www.pasteur-guadeloupe.fr:8081/SITVIT2/", 
     tags$br(), 
     "Brudey, K., Driscoll, J.R., Rigouts, L. et al. Mycobacterium tuberculosis complex genetic diversity: mining the fourth international spoligotyping database (SpolDB4) for classification, population genetics and epidemiology. BMC Microbiol 6, 23 (2006). https://doi.org/10.1186/1471-2180-6-23", 
     tags$br()
   )),
   
   tags$br(),
   img(src="mtbc_tree.png", width='1000px'),
   tags$br(),
   
   
   h4(tags$div(
     "Usage: Input spoligotype into box, e.g. 1111111111111111111111111111000010111110111"
   )), 
   
   # Sidebar with texbox input
   sidebarLayout(
      sidebarPanel(
        textInput(inputId = "spoligotype", 
                  label = "Input spoligotype(s) (comma separated)"),
        actionButton("textSearchButton", "See table")
      ),
      
      # Show a table
      fluidRow(
        column(12, align="center", tableOutput("spol_out"))
      ) # fluidRow
   ) # sidebarLayout

   
) # mainPanel
) # fluidPage

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   # output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      # x    <- faithful[, 2] 
      # bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      # hist(x, breaks = bins, col = 'darkgray', border = 'white')
   # })
  
  spol_df <- eventReactive(input$textSearchButton, {
    # outputFunc(input$textBox, df)
    req(input$spoligotype)
    spols <- unlist(strsplit(input$spoligotype, ","))
    # return(data.frame(Code = spols,
    #                   Description = "Something here",
    #                   Value = "Some value"))
    return(subset(spol_lin_levels_table, Spoligotype %in% spols))
  })
  output$spol_out <- renderTable({
    spol_df()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)


# Test spol:
# 1111111111111111111111111111000010111110111











