# Application Shiny v4
# Créateur: Mathieu Zallio 
# Email: mathieu.zallio@univ-rouen.fr
# Université affiliée : Université de Rouen Normandie 

library(shiny)
library(shinydashboard)
library(plotly)
library(shinyalert)
library(htmlwidgets)
library(ggplot2)
library(shinyjs)



shinyUI(
  dashboardPage(skin="blue", 
    dashboardHeader(title = "My Dashboard"),
    
    #Personnalisation de la barre latérale 
    dashboardSidebar(
      #Ici on va lister toutes les sections, sliders, inputs etc que l'on souhaite avoir sur la barre latérale
      sidebarMenu(
        menuItem("Home", tabName = "Home", icon = icon("home")),
        
        #On laisse à l'utilisateur la possibilité de upload un fichier (dans notre cas on demande un fichier CSV)
        fileInput("fileUpload", "Select a CSV file:"),
        
        #On laisse à l'utilisateur la possibilité de choisir un organisme souhaité que l'on trie par ordre alphabétique.
        selectInput("organisme", "Select an organism name:",
                  choices = sort(
                    c("Mus Musculus","Homo Sapiens", "Gallus gallus domesticus","Pan troglodytes")
                    )),
        
        br(),
        br(),
        
        #Si on veut rajouter une icône à côté du titre d'une la section on utilise l'option icon= icon("élément")
        #Il y a plusieurs icon disponible sur cette librairie : https://getbootstrap.com/docs/3.4/components/#glyphicons 
        # Il y a aussi une autre librairie https://fontawesome.com/icons/categories qui s'appelle font-awesome
        # Pour les utiliser on doit indiquer la librairie utilisée suivie de l'icon
        menuItem("Whole Data Inspection", tabName = "WholeDataInspect", icon = icon("database", lib = "font-awesome")),
        menuItem("Go Term Enrichment", tabName = "GoTermEnrichment", icon = icon("sitemap", lib = "font-awesome"), 
                 badgeLabel = "News", badgeColor = "green"), 
        menuItem("Pathway Enrichment", tabName = "PathwayEnrichment", icon = icon("chart-pie", lib = "font-awesome"), 
                 badgeLabel = "News", badgeColor = "green"), 
        menuItem("About", tabName = "About", icon = icon("th", lib = "glyphicon"))
        )
    ),
    
    # Personnalisation du corps de la page que l'on souhaite
    dashboardBody(
      # Permet d'activer ou désactiver des boutons d'actions 
      useShinyjs(),
      tags$head(
        tags$style(HTML("#parametersBox { font-size: 20px; }"))  # Augmente la taille de la police pour le cadre "Parameters"
      ),
      # Pour personnaliser une page en particulier, on va préciser le tabName associé à la page que l'on a renseigné au dessus
      tabItems(
        tabItem(tabName = "Home",
                
                tags$h2(style = "text-align: center; text-decoration: underline; font-size: 50px;", "Welcome to my Shiny Application"),
                
                br(),
                br(),
                div(
                  style = "max-width: 800px; margin: 0 auto; font-size: 20px;",
                  tags$p("Welcome to this Shiny application. This home page is designed to give you a quick overview of what you can do with the application."),
                
                  tags$h3("How to use the CSV file:"),
                  tags$ul(
                    tags$li(
                      "First, click on the item", shiny::strong("Whole Data Inspection"), "to begin."
                    ),
                    tags$li("Select a CSV file to analyze using the 'Select a CSV file' button."),
                    tags$li("Ensure that the CSV file has the necessary columns: GeneName, ID, baseMean, log2FC, pval, padj.")
                  ),
                
                  tags$h3("Volcano Plot Options:"),
                  tags$ul(
                    tags$li("Use the sliders to set filtering thresholds for p-value and log2 FoldChange."),
                    tags$li("Explore the generated Volcano Plot to visualize upexpressed, downexpressed genes, and filtered data."),
                    tags$li("The table of filtered data is also displayed below the Volcano Plot.")
                  ),
                
                  tags$h3("Downloading Data:"),
                  tags$ul(
                    tags$li("Use the 'Download Volcano Plot & Data' button to download the Volcano Plot and filtered data."),
                    tags$li("The downloaded ZIP file contains an HTML file with the Volcano Plot and a CSV file with the filtered data.")
                  )
                )
        ),
        tabItem(tabName = "WholeDataInspect",
                tags$h2(style = "text-align: center; text-decoration: underline;font-size: 50px;", "Whole Data Inspection"),
                br(),
                
                # Création d'une ligne imaginaire qui va permettre d'aligner les différents éléments que l'on va écrire
                fluidRow(
                  box(
                    title = "Volcano Plot", status = "primary", solidHeader = TRUE,
                    height = 500,
                    plotlyOutput("volcanoPlot")
                  ),
                  
                  # Création d'une box contenant les sliders de filtrage avec un bouton de téléchargement 
                  box(
                    title = "Options", status = "primary", solidHeader = TRUE, 
                    height = 300, 
                    sliderInput("p_value", "P-value cutoff from input:", 0.01, 1, value=0.05),
                    sliderInput("log2Folchange", "log2 FoldChange cutoff from input:", 0, 5.0, value=1, step= 0.5),
                    downloadButton("downloadButtonID", "Download Volcano Plot & Data",class="pull-left",style="margin-right:15px;")
                  )
                ),
                
                fluidRow(
                  box(
                    title = "Datatable when you browse your CSV", status = "primary", solidHeader = TRUE, width = 12,
                    dataTableOutput("dataTable")
                  )
                )
        ),
        tabItem(tabName = "GoTermEnrichment",
              tags$h2("Go Term Enrichment",style = "text-align: center; text-decoration: underline;font-size: 50px;"),
              br(),
              fluidRow(
                box(title = "Parameters", id = "GoTermParametersBox", width = 8, height = 330, solidHeader = TRUE, status = 'primary',
                    fluidRow(
                      column(4, 
                             checkboxGroupInput("GoTermAnalysisType", "Analyse type", choices = c("ORA", "GSEA"))
                      ),
                      column(4, 
                             checkboxGroupInput("OntologiesSelect", "Ontologies", choices = c("Biological process", "Molecular function", "Cellular component"))
                      ),
                      column(4,
                             br(),
                             img(src = "go-logo.large.png", height = "100%", width = "100%")
                      )
                    ),
                    selectInput("GoTermAdjustmentMethod", "Adjustment method", choices = c("Bonferroni", "Another Method")),
                    
                    ######AJOUT D'UN BOUTON RUN ET CLEAR
                    
                    #le bouton run est désactivé tant que aucun fichier n'est upload 
                    actionButton("GoTermRunButton","Run Analysis", icon = icon("step-forward", lib = "glyphicon"), disabled = TRUE),
                    actionButton("clearGoTerm", "Clear Results", icon = icon("remove", lib = "glyphicon"), class = "btn btn-primary", disabled = TRUE)

                ),
                box(
                  title = "Options (for ORA)", width = 4, height = 330, solidHeader = TRUE, status = 'primary',
                  sliderInput("GoTermSliderP_val", "p_Value cutoff from input:", min = 0.01, max = 1, value = 0.1, step = 0.1),
                  br(),
                  sliderInput("GoTermSliderFoldChange", "log2 FoldChange cutoff from input:", min = 0, max = 5.0, value = 0.5, step = 0.5),
                  selectInput("GoTermDegOptions", NULL, choices = c("Over expressed DEG", "Under expressed DEG", "Both"), selected = "Over expressed DEG")
                )
              ),
              
              uiOutput("GoTermAnalysisBox")
              
      ),
        tabItem(tabName = "PathwayEnrichment",
                tags$h2(style = "text-align: center; text-decoration: underline; font-size: 50px;", "Pathway Enrichment"),
                br(),
                fluidRow(
                  box(title = "Parameters", id = "PathwayParametersBox", width = 8, height = 330, solidHeader = TRUE, status = 'primary',
                      fluidRow(
                        column(4, 
                               checkboxGroupInput("PathwayAnalysisType", "Analyse type", choices = c("ORA", "GSEA"))
                        ),
                        column(4, 
                               checkboxGroupInput("PathwayDatabase", "Database", choices = c("KEGG", "Reactome"))
                        ),
                        column(4,
                               img(src = "kegg128.gif", height = "40%", width = "40%"),
                               img(src = "Reactome_Imagotype_Positive_100mm.png", height = "66%", width = "66%")
                        )
                      ),
                      selectInput("PathwayAdjustmentMethod", "Adjustment method", choices = c("Bonferroni", "Another Method")),
                      
                      ######AJOUT D'UN BOUTON RUN ET CLEAR
                      actionButton("PathwayRunButton","Run Analysis", icon = icon("step-forward", lib = "glyphicon"), disabled = TRUE),
                      actionButton("clearPathway", "Clear Results", icon = icon("remove", lib = "glyphicon"), class = "btn btn-primary", disabled = TRUE)
                  ),
                  box(
                    title = "Options (for ORA)", width = 4, height = 330, solidHeader = TRUE, status = 'primary',
                    sliderInput("PathwaySliderP_val", "p_Value cutoff from input:", min = 0.01, max = 1, value = 0.1, step = 0.1),
                    br(),
                    sliderInput("PathwaySliderFoldChange", "log2 FoldChange cutoff from input:", min = 0, max = 5.0, value = 0.5, step = 0.5),
                    selectInput("PathwayDegOptions", NULL, choices = c("Over expressed DEG", "Under expressed DEG", "Both"), selected = "Over expressed DEG")
                  )
                ),
                
                uiOutput("PathwayAnalysisBox")
                
        ),
        tabItem(tabName = "About",
                tags$h2(style = "text-align: center; text-decoration: underline; font-size: 50px;", "About"),
                br(),
                br(),
                div(
                  style = "font-size: 18px;",
                  tags$p(shiny::strong("Version :"), "3.0"),
                  tags$p(shiny::strong("Creators :"), "Baptiste Herlemont / Adam Schumacher / Mathieu Zallio"),
                  tags$p(shiny::strong("Email :"), "baptiste.herlemont@univ-rouen.fr / adam.schumacher@univ-rouen.fr / mathieu.zallio@univ-rouen.fr"),
                  tags$p(shiny::strong("University :"), "Université de Rouen Normandie ")
                  )
                )
        )
      )
    )
  )




