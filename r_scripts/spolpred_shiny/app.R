# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ape)
library(phytools)
library(ggplot2)
library(ggtree)
library(gridExtra)

# Functions ----

# Split input on comma (default) to provide vector - allows multiple selection in textbox
to_vect <- function(x, split_on = ","){
  unlist(strsplit(x, split_on))
}

# Data source ----

data_url <- "https://raw.githubusercontent.com/GaryNapier/spolpred/master/results/spol_43_table_github.csv"

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
    "Lineages are in a numerical hierarchy - 1, 1.1, 1.1.1 etc. - which have been given a 'Level' - 1, 2, 3 etc - e.g.:",
    tags$br(),
    tags$br(),
    plotOutput("levels_plot"), 
    tags$br(),
    tags$br(),
    tags$br(),
    tags$br(),
    tags$br(),
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
    tags$br(),
    style = 'width:1000px', align = "justify"
  )),
  
  h3(tags$div(
    "Main results:",
    tags$br()
  )),
  
  h4(tags$div(
    "The plots below show the relationships between the SNP-based lineages and the spoligotypes.
    The variation in '% in level-lineage' and 'Correlation' columns are plotted by lineage level and grouped 
    by main lineage. 
    The main message is that variation in both measures increases at the lower lineage levels, 
    reflecting noise in the spoligotypes.",
    style = 'width:1000px', align = "justify",
    tags$br()
  )),
  
  tags$br(),
  img(src="spol_43_boxplots.png", width='1000px'),
  tags$br(),
  
  tags$br(),
  img(src="spol_68_boxplots.png", width='1000px'),
  tags$br(),
  
  h5(tags$div(
    tags$br(),
    tags$br(),
    "Top panels: The maximum percentage of samples with a spoligotype at that lineage-level 
    ('% in level-lineage' column) increases as the lineage level decreases (1-4). 
    At the same time, there is more variation in these percentages at the lower lineage levels.
    For example, for lineage 1, level 1, the highest percentage of samples is 17.2, for spoligotype SIT-19 (n=402).
    Then, at lowest level (4), 100% of samples in lineage 1.1 have spoligotype SIT-939, 
    however, there are only n=14 samples in this lineage-level. 
    Bottom panels: At the lower lineage levels (1-4), there remain spoligotypes with correlation 1, 
    but the variation in correlation scores increases, reflecting noise in the spoligotypes,
    Red diamond = maximum values.", 
    style = 'width:1000px', align = "justify"
  )),
  
  h3(tags$div(
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
    tags$br(), 
    style = 'width:1000px', align = "justify"
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
      
      tabsetPanel(
        # tabPanel("Expression values", tableOutput("mainTable")),
        # tabPanel("ID filtering", tableOutput("table"))
        # tableOutput("spol_out"),
        # downloadButton('downloadData', 'Download')
        tabPanel("43-spacer spol", tableOutput("spol_out_43"), 
                 downloadButton('downloadData_43', 'Download')),
        tabPanel("68-spacer spol", tableOutput("spol_out_68"), 
                 downloadButton('downloadData_68', 'Download'))
      ) # tabsetPanel
    )
  ) # sidebarLayout
)

server <- function(input, output) {
  
  # Levels illustration ----
  
  n <- 20
  title_sz <- 24
  set.seed(1)
  
  tree <- midpoint.root(rcoal(n = n, tip.label = LETTERS[1:n]))
  tree$tip.label <- sort(tree$tip.label)
  
  dat_lv1 <- data.frame(id = c(22, 27), Lineage = c("Lineage 1", "Lineage 2"))
  lv1 <- ggtree(tree)+
    geom_tiplab()+
    # geom_text(aes(label=node), hjust=2, vjust = -1)+
    geom_hilight(data = dat_lv1, 
                 mapping = aes(node = id, fill = Lineage))+
    # scale_fill_manual(values=c("blue", "darkred"))+
    ggtitle("Level 1")+
    theme(plot.title = element_text(size = title_sz))
  
  dat_lv2 <- data.frame(id = c(38, 28), Lineage = c("Lineage 2.1", "Lineage 2.2"))
  lv2 <- ggtree(tree)+
    geom_tiplab()+
    # geom_text(aes(label=node), hjust=2, vjust = -1)+
    geom_hilight(data = dat_lv2, 
                 mapping = aes(node = id, fill = Lineage))+
    ggtitle("Level 2")+
    theme(plot.title = element_text(size = title_sz))
  
  dat_lv3 <- data.frame(id = c(36, 29), Lineage = c("Lineage 2.2.1", "Lineage 2.2.2"))
  lv3 <- ggtree(tree)+
    geom_tiplab()+
    # geom_text(aes(label=node), hjust=2, vjust = -1)+
    geom_hilight(data = dat_lv3, 
                 mapping = aes(node = id, fill = Lineage))+
    ggtitle("Level 3")+
    theme(plot.title = element_text(size = title_sz))
  
  dat_lv4 <- data.frame(id = c(35, 30), Lineage = c("Lineage 2.2.2.1", "Lineage 2.2.2.2"))
  lv4 <- ggtree(tree)+
    geom_tiplab()+
    # geom_text(aes(label=node), hjust=2, vjust = -1)+
    geom_hilight(data = dat_lv4, 
                 mapping = aes(node = id, fill = Lineage))+
    ggtitle("Level 4")+
    theme(plot.title = element_text(size = title_sz))
  
  output$levels_plot <- renderPlot({
    grid.arrange(lv1, lv2, lv3, lv4, nrow = 2)
  }, height = 500, width = 700)
  
  
  # Spol table ----
  
  spol_data_43_url <- "https://raw.githubusercontent.com/GaryNapier/spolpred/master/results/spol_lin_levels_table_github.csv"
  spol_data_68_url <- "https://raw.githubusercontent.com/GaryNapier/spolpred/master/results/spol_lin_levels_table_68_github.csv"
  
  spol_43_table <- read.csv(spol_data_43_url, colClasses = c("Spoligotype" = "character"))
  spol_68_table <- read.csv(spol_data_68_url, colClasses = c("Spoligotype" = "character"))
  
  names(spol_43_table) <- c("Level",
                            "Lineage",
                            "Spoligotype",
                            "SIT",
                            "Family",
                            "Correlation",
                            "N",
                            "% in level-lineage")
  
  names(spol_68_table) <- c("Level",
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

  # 68-spacer
  
  level_values <- reactive({
    if (isTruthy(input$level)) {
      return(to_vect(input$level))
    } else {
      return(unique(spol_43_table$Level))
    }
  })
  
  lineage_values <- reactive({
    if (isTruthy(input$lineage)) {
      return(to_vect(input$lineage))
    } else {
      return(unique(spol_43_table$Lineage))
    }
  })
  
  spol_values <- reactive({
    if (isTruthy(input$spoligotype)) {
      return(to_vect(input$spoligotype))
    } else {
      return(unique(spol_43_table$Spoligotype))
    }
  })
  
  sit_values <- reactive({
    if (isTruthy(input$SIT)) {
      return(to_vect(input$SIT))
    } else {
      return(unique(spol_43_table$SIT))
    }
  })
  
  family_values <- reactive({
    if (isTruthy(input$family)) {
      return(to_vect(input$family))
    } else {
      return(unique(spol_43_table$Family))
    }
  })
  
  level_values <- reactive({
    if (isTruthy(input$level)) {
      return(to_vect(input$level))
    } else {
      return(unique(spol_43_table$Level))
    }
  })
  
  # 68-spacer
  
  level_values_68 <- reactive({
    if (isTruthy(input$level)) {
      return(to_vect(input$level))
    } else {
      return(unique(spol_68_table$Level))
    }
  })
  
  lineage_values_68 <- reactive({
    if (isTruthy(input$lineage)) {
      return(to_vect(input$lineage))
    } else {
      return(unique(spol_68_table$Lineage))
    }
  })
  
  spol_values_68 <- reactive({
    if (isTruthy(input$spoligotype)) {
      return(to_vect(input$spoligotype))
    } else {
      return(unique(spol_68_table$Spoligotype))
    }
  })
  
  sit_values_68 <- reactive({
    if (isTruthy(input$SIT)) {
      return(to_vect(input$SIT))
    } else {
      return(unique(spol_68_table$SIT))
    }
  })
  
  family_values_68 <- reactive({
    if (isTruthy(input$family)) {
      return(to_vect(input$family))
    } else {
      return(unique(spol_68_table$Family))
    }
  })
  
  
  filter_df_43 <- reactive({
    return(spol_43_table %>% 
             filter(Level %in% level_values(),
                    Lineage %in% lineage_values(),
                    Spoligotype %in% spol_values(),
                    SIT %in% sit_values(),
                    Family %in% family_values()
                    ) # filter
    ) # return
  }) # reactive
  
  filter_df_68 <- reactive({
    return(spol_68_table %>% 
             filter(Level %in% level_values_68(),
                    Lineage %in% lineage_values_68(),
                    Spoligotype %in% spol_values_68(),
                    SIT %in% sit_values_68(),
                    Family %in% family_values_68()
             ) # filter
    ) # return
  }) # reactive
  
  output$spol_out_43 <- renderTable({filter_df_43()})
  output$spol_out_68 <- renderTable({filter_df_68()}) 
  
  output$downloadData_43 <- downloadHandler(
    filename = function() {
      paste("spol_correlation_43_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(filter_df_43(), file, row.names = F) # call to the reactive again
    }
  )
  
  output$downloadData_68 <- downloadHandler(
    filename = function() {
      paste("spol_correlation_68_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(filter_df_68(), file, row.names = F) # call to the reactive again
    }
  )
  
} # server()

# Run the application
shinyApp(ui = ui, server = server)



