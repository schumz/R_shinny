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
library(enrichplot)
library(GOSemSim)
library(genekitr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(DOSE)
library(GOSemSim)


#################################################################################################################################################################################################################"

ontologies_associations <- list(
  "Biological process" = "BP",
  "Molecular function" = "MF",
  "Cellular component" = "CC"
)


#################################################################################################################################################################################################################"
GO_ora_analysis <- function(filtered_data, ontologies_terms, org_db = org.Mm.eg.db, simplifyCutoff = 0.7) {
  
  # Mappage ENSEMBL à Entrez
  entrez_ids <- mapIds(org_db,
                       keys = filtered_data$ID,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  
  # Initialiser une liste pour stocker les résultats simplifiés de chaque ontologie
  simplified_results_list <- list()
  
  for (term in ontologies_terms) {
    # Utiliser la liste nommée pour obtenir l'abréviation correspondante
    ont <- ontologies_associations[[term]]
    message = sprintf("ORA pour %s en cours",ont)
    print(message)
    
    
    tryCatch({
      go_results <- enrichGO(gene = entrez_ids, OrgDb = org_db, keyType = "ENTREZID", ont = ont)
      go_results_simplified <- simplify(go_results, cutoff = simplifyCutoff)
      
      # Ajouter les résultats simplifiés à la liste
      simplified_results_list[[ont]] <- go_results_simplified
    }, error = function(e) {
      cat("Error in ORA for", ont, ": ", e$message, "\n")
    })
  }
  
  # Retourner la liste des résultats simplifiés
  return(simplified_results_list)
}


GO_gsea_analysis <- function(data, ontologies_terms, org_db = org.Mm.eg.db, eps = 1e-300, pAdjustMethod = "BH",simplifyCutoff = 0.7) {
  data_sorted <- data %>%
    filter(!duplicated(ID)) %>%
    filter(!is.na(log2FC) & !is.na(ID) & !is.na(padj)) %>%
    arrange(desc(log2FC))
  
  # Mappage ENSEMBL à Entrez
  entrez_ids <- mapIds(org_db,
                       keys = data_sorted$ID,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  
  # Ajout de la colonne EntrezID aux données triées
  data_sorted <- data_sorted %>%
    mutate(EntrezID = entrez_ids)
  
  # Préparation de la liste de gènes pour GSEA
  gene_GSEA_sorted <- setNames(data_sorted$log2FC, data_sorted$EntrezID)
  
  simplified_results_list <- list()
  
  for (term in ontologies_terms) {
    ont <- ontologies_associations[[term]]
    message = sprintf("GSEA en cours pour %s en cours",ont)
    print(message)
    
    tryCatch({
      gse_go_result <- gseGO(geneList = gene_GSEA_sorted, ont = ont, keyType = "ENTREZID", OrgDb = org_db, eps = eps, pAdjustMethod = pAdjustMethod)
      gse_go_results_simplified <- simplify(gse_go_result, cutoff = simplifyCutoff)
      
      simplified_results_list[[ont]] <- gse_go_results_simplified
    }, error = function(e) {
      cat("Error in GSEA for", ont, ": ", e$message, "\n")
    })
  }
  
  return(simplified_results_list)
}




kegg_ora_analysis <- function(filtered_data, organism = "mmu", key_type = "ENSEMBL") {
  
  gene_ORA <- as.character(filtered_data$ID)
  
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = gene_ORA,
                       column = "ENTREZID",
                       keytype = key_type,
                       multiVals = "first")
  
  
  tryCatch({
    kegg_results <- enrichKEGG(gene = entrez_ids, organism = organism, keyType = "ncbi-geneid")
  }, error = function(e) {
    cat("Error in KEGG ORA:", e$message, "\n")
  })
  return(kegg_results)
}


kegg_gsea_analysis <- function(data, organism = "mmu", key_type = "ENSEMBL", eps = 1e-300, pAdjustMethod = "BH") {
  
  data_sorted <- data %>%
    filter(!duplicated(ID)) %>%
    filter(!is.na(log2FC) & !is.na(ID) & !is.na(padj)) %>%
    arrange(desc(log2FC))
  gene_GSEA <- setNames(data_sorted$log2FC, data_sorted$ID)
  
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = names(gene_GSEA),
                       column = "ENTREZID",
                       keytype = key_type,
                       multiVals = "first")
  
  # Créez un nouveau vecteur avec les identifiants Entrez mappés, mais conservez l'ordre de gene_GSEA
  entrez_geneList <- setNames(gene_GSEA[names(gene_GSEA) %in% names(entrez_ids)], entrez_ids[names(gene_GSEA) %in% names(entrez_ids)])
  
  # Assurez-vous que le vecteur est trié en ordre décroissant
  entrez_geneList_sorted <- sort(entrez_geneList, decreasing = TRUE)
  
  tryCatch({
    gse_kegg_result <- gseKEGG(geneList = entrez_geneList_sorted, organism = organism, keyType = "ncbi-geneid", eps = eps, pAdjustMethod = pAdjustMethod)
  }, error = function(e) {
    cat("Error in KEGG GSEA:", e$message, "\n")
  })
}


#################################################################################################################################################################################################################"


shinyServer(function(input, output) {
  raw_data <- reactiveVal()
  #############################################################################################################################################################################################################################################################"
  ######################################################## PARTIE WHOLE DATA INSPECTION: CHECK DE L'INPUT + VOLCANOPLOT + DOWNLOAD + ACTIVATION DES BOUTONS RUN DE GO ET PATHWAY ##############################################################################"
  #############################################################################################################################################################################################################################################################"
  
  #Lorsque l'on va "observer l'évènement " : un fichier a été mis dans le fileUpload, on va analyser le fichier
  observeEvent(input$fileUpload, {
    #ici, on vérifie que le fichier est bien upload, puis on le lit avec des lignes de commandes R classiques
    if (!is.null(input$fileUpload)) {
      
      #Vérification que c'est bien un .csv, si ce n'est pas le cas on émet une alerte
      if (!grepl("\\.csv$", input$fileUpload$name, ignore.case = TRUE)) {
        shinyalert(title = "Error in file format",
                   text ="Please upload a valid CSV file.",
                   type= "error")
        return(NULL)
      }
      
      # Si c'est bien un .csv on va charger les données dans une variable
      data <- read.csv(input$fileUpload$datapath, sep = ";")
      raw_data(data)
      
      #Vérification que le csv possède bien les bonnes colonnes et uniquement ces colonnes, si ce n'est pas le cas on émet une alerte
      required_columns <- c("GeneName", "ID", "baseMean", "log2FC", "pval", "padj")
      if (!all(required_columns %in% names(data)) || !all(names(data) %in% required_columns)) {
        shinyalert(title = "Error in column names", 
                   text = "The CSV file must contain the following columns: GeneName ; ID ; baseMean ; log2FC ; pval ; padj",
                   type= "error")
        return(NULL)
      }
      
      ## On va rendre clickable tout les boutons des onglets Pathway et Go analysis 
      shinyjs::enable("GoTermRunButton")
      shinyjs::enable("clearGoTerm")
      shinyjs::enable("PathwayRunButton")
      shinyjs::enable("clearPathway")
      
      # Si le .csv passe les 2 tests et est donc conforme à nos attentes, on va pouvoir débuter l'analyse
      
      ######################## Filtrage des données en fonction des sliders #############################
      filtered_data <- reactive({
        #on vérifie que les 2 inputs sont bien présents 
        req(input$p_value, input$log2Folchange)
        filtered <- data %>%
          filter(pval < input$p_value & abs(log2FC) > input$log2Folchange)
        
        return(filtered)
      })
      ##################################################################################################
      
      
      ######################## Création du volcano plot en fonction des sliders #############################
      plot_filtered <- reactive({
        #on vérifie que les 2 inputs sont bien présents 
        req(input$p_value, input$log2Folchange)
        
        #Ici on va colorier les points du Volcano plot en fonction de leur pvalue et de leur Log2FC 
        #Même si on utilise les données complète (data), le volcano va refléter les données présentes sur le filtered_data
        #Les points coloriés en rouge ou vert seront présent sur le tableau filtré mais les points gris eux seront supprimés du tableau fitré
        #Le fait de séparer l'apparition des points avec add_trace permet de bien faire aparaître la légende en fonction de la couleur des points
        p <- plot_ly(data = data[data$log2FC > input$log2Folchange & data$pval < input$p_value ,], x = ~log2FC, y = ~-log10(pval), 
                     type = 'scatter',
                     mode = 'markers',
                     name = 'Upexpressed',
                     #cette option permet, lors du passage de la souris, d'afficher les informations associées aux points
                     text = ~paste("Gene: ", GeneName,
                                   "<br>ID: ", ID,
                                   "<br>BaseMean: ", baseMean, 
                                   "<br>Log2FC: ", log2FC, 
                                   "<br>p-value: ", pval, 
                                   "<br>Adjusted p-value: ", padj),
                     marker = list(color = "green", size = 5)) %>%
          add_trace(data = data[data$log2FC < -input$log2Folchange & data$pval < input$p_value,], x = ~log2FC, y = ~-log10(pval), 
                    type = 'scatter',
                    mode = 'markers',
                    name = 'Downexpressed',
                    text = ~paste("Gene: ", GeneName,
                                  "<br>ID: ", ID,
                                  "<br>BaseMean: ", baseMean, 
                                  "<br>Log2FC: ", log2FC, 
                                  "<br>p-value: ", pval, 
                                  "<br>Adjusted p-value: ", padj),
                    marker = list(color = "red", size = 5)) %>%
          add_trace(data = data[data$log2FC < input$log2Folchange & data$pval > input$p_value | data$log2FC > -input$log2Folchange & data$pval > input$p_value | data$log2FC > -input$log2Folchange & data$log2FC < input$log2Folchange ,], x = ~log2FC, y = ~-log10(pval), 
                    type = 'scatter',
                    mode = 'markers',
                    name = 'Filtered Data',
                    text = ~paste("Gene: ", GeneName,
                                  "<br>ID: ", ID,
                                  "<br>BaseMean: ", baseMean, 
                                  "<br>Log2FC: ", log2FC, 
                                  "<br>p-value: ", pval, 
                                  "<br>Adjusted p-value: ", padj),
                    marker = list(color = "grey", size = 5)) %>%
          layout(title = '-log(pval) vs log2FC',
                 xaxis = list(title = 'log2 Fold Change'),
                 yaxis = list(title = '-log10(p-value)') )
        
        return(p)
      })
      ############################################################################################################

      ######################## Génération des Outputs #############################
      #On va réutiliser les variables réactives crées auparavant 
      
      #Affichage du tableau de données filtrés
      output$dataTable <- renderDataTable({ filtered_data()})
      
      #Affichage du volcano plot reflétant les données filtrées
      output$volcanoPlot <- renderPlotly({plot_filtered()})
      
      #############################################################################
      
      
      ######################## Mise en action du bouton Download #############################
      
      #On va sauvegarder nos 2 éléments importants: le volcano plot et le tableau associé.
      #On va donc créer une archive zip avec les 2 fichiers (au format html et csv)
      #L'archive sera nommée de la façon suivante: Downloaded_date_with_x_log2FC_cutoff_y_pval_cutoff.zip
      # Avec x=valeur indiquée par le slider log2FC, y=valeur indiquée par le slider pval et date=date du téléchargement  
      
      output$downloadButtonID <- downloadHandler(
        filename = function() {
          paste("Downloaded_",Sys.Date(),"_with_",input$log2Folchange,"_log2FC_cutoff_",input$p_value,"_pval_cutoff.zip", sep = "")
        },
        content = function(file) {
          
          # Sauvegarde en fichier html du volcao plot
          NameVolcano <- paste("Volcano_plot_",input$log2Folchange,"_log2FC_cutoff_", input$p_value, "_pval_cutoff.html", sep = "")
          saveWidget(plot_filtered(), NameVolcano , selfcontained = TRUE)
          
          # Sauvegarde en fichier CSV des données filtrées
          NameFilterData <- paste("Filtered_data_",input$log2Folchange,"_log2FC_cutoff_", input$p_value, "_pval_cutoff.csv", sep = "")
          write.csv(filtered_data(), NameFilterData, row.names = FALSE)
          
          # Création de l'archive ZIP
          zip(file, files = c( NameVolcano , NameFilterData ))
        }
      )
      
      #La fonctionnalité de téléchargement pourra évoluer par la suite 
      #Par exemple, laisser le choix à l'utilisateur de quel élément (html ou csv) il souhaite télécharger
      #À ce stade, j'ai préféré télécharger les 2 éléments en même temps
      
      #######################################################################################
      
    }
  })
  #############################################################################################################################################################################################################################################################"
  #############################################################################################################################################################################################################################################################"
  #############################################################################################################################################################################################################################################################"
  
  
  ##### SERA SUPPRIME
  ###### FIGURES TESTS POUR LES PLACER :
  output$barplot1 <- renderPlot({
    # Données exemple pour le premier barplot
    data1 <- data.frame(
      category = c("A", "B", "C"),
      value = c(2, 3, 5)
    )
    ggplot(data1, aes(x = category, y = value, fill = category)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = "Barplot 1", x = "Category", y = "Value")
  })
  
  output$barplot2 <- renderPlot({
    # Données exemple pour le deuxième barplot
    data2 <- data.frame(
      category = c("D", "E", "F"),
      value = c(4, 1, 7)
    )
    ggplot(data2, aes(x = category, y = value, fill = category)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(title = "Barplot 2", x = "Category", y = "Value")
  })
  
  # Plot 1
  output$plot1 <- renderPlot({
    ggplot(mtcars, aes(x = mpg, y = disp)) + geom_point() + ggtitle("Plot 1: MPG vs Displacement")
  })
  
  # Plot 2
  output$plot2 <- renderPlot({
    ggplot(mtcars, aes(x = factor(cyl), fill = factor(gear))) + geom_bar(position = "dodge") + ggtitle("Plot 2: Cylinder Count by Gear Type")
  })
  
  # Plot 3
  output$plot3 <- renderPlot({
    ggplot(mtcars, aes(x = wt, y = qsec)) + geom_line() + ggtitle("Plot 3: Weight vs Quarter Mile Time")
  })
  
  # Plot 4
  output$plot4 <- renderPlot({
    ggplot(mtcars, aes(x = hp)) + geom_histogram(binwidth = 20) + ggtitle("Plot 4: Distribution of Horsepower")
  })
  
  # Plot 5
  output$plot5 <- renderPlot({
    ggplot(mtcars, aes(x = drat, y = carb)) + geom_point(aes(color = factor(am))) + ggtitle("Plot 5: Rear Axle Ratio vs Number of Carburetors")
  })
  
  # Plot 6
  output$plot6 <- renderPlot({
    ggplot(mtcars, aes(x = gear, y = mpg)) + geom_boxplot() + ggtitle("Plot 6: MPG by Gear Count")
  })
  
  # Plot 7
  output$plot7 <- renderPlot({
    ggplot(mtcars, aes(x = vs, y = mpg, color = factor(vs))) + geom_jitter() + ggtitle("Plot 7: MPG by Engine Shape")
  })
  
  # Plot 8
  output$plot8 <- renderPlot({
    ggplot(mtcars, aes(x = factor(cyl), y = mpg)) + geom_violin() + ggtitle("Plot 8: MPG Distribution by Cylinder Count")
  })
  
  
  
  
  
  #############################################################################################################################################################################################################################################################"
  ######################################################## PARTIE GO TERM ANALYSIS: CHECK DU COCHAGE POUR ANALYSE + APPARITION DES BOXS (AVEC MAUVAIS PLOTS) ##############################################################################"
  #############################################################################################################################################################################################################################################################"
  GoAnalysisReady <- reactiveVal(FALSE) ## Création d'une variable invisible permettant d'activer l'analyse ou non selon les conditions 
  
  observeEvent(input$GoTermRunButton, {
    # Vérifier si les sélections ont été faites dans "Analyse type" et "Ontologies"
    if (length(input$GoTermAnalysisType) == 0 || length(input$OntologiesSelect) == 0) {
      shinyalert::shinyalert(title = "Selection Required",
                             text = "Please select at least one option in both 'Analyse Type' and 'Ontologies'.",
                             type = "warning")
      GoAnalysisReady(FALSE)  # Indiquer que l'analyse n'est pas prête
      return()
    }
    
    # Si les conditions sont remplies
    GoAnalysisReady(TRUE)  # Indiquer que l'analyse est prête
  })
  
 
  #### Pas de height pour les box au moins elles vont s'adapter à la taille de la figure  quand on les intègrera 
  output$GoTermAnalysisBox <- renderUI({
    if (GoAnalysisReady()) {
      
      GO_choix_analyse <- input$GoTermAnalysisType
      ontologies <- input$OntologiesSelect
      
      GO_seuil_pval <- input$GoTermSliderP_val
      GO_seuil_fc <- input$GoTermSliderFoldChange
      GO_ora_option <- input$GoTermDegOptions
      
      
      
      if ("ORA" %in% GO_choix_analyse) {
       print("c'est good")
       data_test <- raw_data() 
       filtered_data_ora <- subset(data_test,padj <= seuil_pval & abs(log2FC) >= GO_seuil_fc)
       if (GO_ora_option == "Over expressed DEG") {
         # Gardez uniquement les lignes avec log2FC > 0 pour les gènes sur-exprimés
         filtered_data_ora <- filtered_data_ora[filtered_data_ora$log2FC > 0, ]
       } else if (GO_ora_option == "Under expressed DEG") {
         # Gardez uniquement les lignes avec log2FC < 0 pour les gènes sous-exprimés
         filtered_data_ora <- filtered_data_ora[filtered_data_ora$log2FC < 0, ]
       }
       res_GO_ora_globaux = GO_ora_analysis(filtered_data_ora,ontologies,org_db = org.Mm.eg.db,simplifyCutoff = 0.7)
      }
      
      if ("GSEA" %in% GO_choix_analyse) {
        
       data_test <- raw_data() 
       res_GO_gsea_globaux = GO_gsea_analysis(data_test,ontologies,org_db = org.Mm.eg.db,simplifyCutoff = 0.7)
       print("c'est good pour la GSEA")
       print(length(res_gsea_globaux))
       }
      
      
      
      
      if ("ORA" %in% GO_choix_analyse && "GSEA" %in% GO_choix_analyse) {
        fluidRow(
          box(title = "BarPlot", id = "GoTermOraGseaBarplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              fluidRow(
                column(6, plotOutput("barplot1")), # Pour le premier plot
                column(6, plotOutput("barplot2"))  # Pour le deuxième plot
              )
              ),
          box(title = "CnetPlot", id = "GoTermOraGseaCnetplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              fluidRow(
                column(6, plotOutput("plot1")), 
                column(6, plotOutput("plot2"))  
              )
              ),
          box(title = "EmaPlot", id = "GoTermOraGseaEmaplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              fluidRow(
                column(6, plotOutput("plot3")), 
                column(6, plotOutput("plot4")) 
              )
              ),
          box(title = "DotPlot", id = "GoTermOraGseaDotplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              fluidRow(
                column(6, plotOutput("plot5")),
                column(6, plotOutput("plot6"))  
              )
              ),
          box(title = "GoPlot & GseaPlot", id = "GoTermGoGseaplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              fluidRow(
                column(6, plotOutput("plot7")),
                column(6, plotOutput("plot8"))  
              )
              ),
          shinyjs::disable("GoTermRunButton"),
          shinyjs::disable("GoTermAnalysisType"),
          shinyjs::disable("OntologiesSelect")
        ) }
      else if ("ORA" %in% GO_choix_analyse) { 
        fluidRow(
          box(title = "BarPlot", id = "GoTermBarplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              # Contenu de la boîte des résultats ici
              plotOutput("barplot1")
            ),
          box(title = "CnetPlot", id = "GoTermOraCnetplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot1")
            ),
          box(title = "EmaPlot", id = "GoTermOraEmaplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot3")
            ),
          box(title = "DotPlot", id = "GoTermOraDotplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot5")
            ),
          box(title = "GoPlot ", id = "GoTermGoplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot7")
            ),
          shinyjs::disable("GoTermRunButton"),
          shinyjs::disable("GoTermAnalysisType"),
          shinyjs::disable("OntologiesSelect")
        ) } 
      else if ("GSEA" %in% GO_choix_analyse) { 
        fluidRow( 
          box(title = "BarPlot", id = "GoTermGseaBarplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("barplot2")
          ),
          box(title = "CnetPlot", id = "GoTermGseaCnetplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot2")
          ),
          box(title = "EmaPlot", id = "GoTermGseaEmaplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot4")
          ),
          box(title = "DotPlot", id = "GoTermGseaDotplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot6")
          ),
          box(title = "GseaPlot ", id = "GoTermGseaplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot8")
          ),
          shinyjs::disable("GoTermRunButton"),
          shinyjs::disable("GoTermAnalysisType"),
          shinyjs::disable("OntologiesSelect")
        )}
    }
  })
  
  observeEvent(input$clearGoTerm, {
    # Réinitialiser GoAnalysisReady à FALSE pour "effacer" les box
    GoAnalysisReady(FALSE)
    
    # Réactiver le bouton Run si nécessaire
    shinyjs::enable("GoTermRunButton")
    shinyjs::enable("GoTermAnalysisType")
    shinyjs::enable("OntologiesSelect")
  })
  
  #############################################################################################################################################################################################################################################################"
  #############################################################################################################################################################################################################################################################"
  #############################################################################################################################################################################################################################################################"
 
  
  #############################################################################################################################################################################################################################################################"
  ######################################################## PARTIE PATHWAY ANALYSIS: CHECK DU COCHAGE POUR ANALYSE + APPARITION DES BOXS (VIDES POUR L'INSTANT) ##############################################################################"
  #############################################################################################################################################################################################################################################################"
  PathwayAnalysisReady <- reactiveVal(FALSE) ## Création d'une variable invisible permettant d'activer l'analyse ou non selon les conditions 
  
  observeEvent(input$PathwayRunButton, {
    # Vérifier si les sélections ont été faites dans "Analyse type" et "Ontologies"
    if (length(input$PathwayAnalysisType) == 0 || length(input$PathwayDatabase) == 0) {
      shinyalert::shinyalert(title = "Selection Required",
                             text = "Please select at least one option in both 'Analyse Type' and 'Database'.",
                             type = "warning")
      PathwayAnalysisReady(FALSE)  # Indiquer que l'analyse n'est pas prête
      return()
    }
    
    # Si les conditions sont remplies
    PathwayAnalysisReady(TRUE)  # Indiquer que l'analyse est prête
  })
  
  
  #### Pas de height pour les box au moins elles vont s'adapter à la taille de la figure  quand on les intègrera 
  output$PathwayAnalysisBox <- renderUI({
    if (PathwayAnalysisReady()) {
      
      
      pathways_choix_analyse <- input$PathwayAnalysisType
      databases <- input$PathwayDatabase
      
      pathways_seuil_pval <- input$PathwaySliderP_val
      pathways_seuil_fc <- input$PathwaySliderFoldChange
      pathways_ora_option <- input$PathwayDegOptions
      
      
      if ("ORA" %in% pathways_choix_analyse) {
        print("c'est good")
        data_test <- raw_data() 
        filtered_data_ora_pathways <- subset(data_test,padj <= seuil_pval & abs(log2FC) >= GO_seuil_fc)
        if (pathways_ora_option == "Over expressed DEG") {
          # Gardez uniquement les lignes avec log2FC > 0 pour les gènes sur-exprimés
          filtered_data_ora <- filtered_data_ora[filtered_data_ora$log2FC > 0, ]
        } else if (pathways_ora_option == "Under expressed DEG") {
          # Gardez uniquement les lignes avec log2FC < 0 pour les gènes sous-exprimés
          filtered_data_ora <- filtered_data_ora[filtered_data_ora$log2FC < 0, ]
        }
        res_kegg_ora = kegg_ora_analysis(filtered_data, organism = "mmu", key_type = "ENSEMBL")
      }
      
      if ("GSEA" %in% pathways_choix_analyse) {
        
        data_test <- raw_data() 
        res_kegg_gsea = kegg_gsea_analysis(data_test,organism = "mmu", key_type = "ENSEMBL", eps = 1e-300, pAdjustMethod = "BH")
        print("c'est good pour la GSEA")
      }
      
      
      
      
      
      
      
      
      
      if ("ORA" %in% input$PathwayAnalysisType && "GSEA" %in% input$PathwayAnalysisType) {
        fluidRow(
          # Remplacez ceci par le contenu que vous souhaitez afficher dans la box
          box(title = "BarPlot", id = "PathwayOraGseaBarplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              fluidRow(
                column(6, plotOutput("barplot1")), # Pour le premier plot
                column(6, plotOutput("barplot2"))  # Pour le deuxième plot
              )
          ),
          box(title = "CnetPlot", id = "PathwayOraGseaCnetplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              fluidRow(
                column(6, plotOutput("plot1")), 
                column(6, plotOutput("plot2"))  
              )
          ),
          box(title = "EmaPlot", id = "PathwayOraGseaEmaplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              fluidRow(
                column(6, plotOutput("plot3")), 
                column(6, plotOutput("plot4")) 
              )
          ),
          box(title = "DotPlot", id = "PathwayOraGseaDotplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              fluidRow(
                column(6, plotOutput("plot5")),
                column(6, plotOutput("plot6"))  
              )
          ),
          box(title = "GoPlot & GseaPlot", id = "PathwayGoGseaplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              fluidRow(
                column(6, plotOutput("plot7")),
                column(6, plotOutput("plot8"))  
              )
          ),
          shinyjs::disable("PathwayRunButton"),
          shinyjs::disable("PathwayAnalysisType"),
          shinyjs::disable("PathwayDatabase")
        ) }
      else if ("ORA" %in% input$PathwayAnalysisType) { 
        fluidRow(
          box(title = "BarPlot", id = "PathwayOraBarplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              # Contenu de la boîte des résultats ici
              plotOutput("barplot1")
          ),
          box(title = "CnetPlot", id = "PathwayOraCnetplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot1")
          ),
          box(title = "EmaPlot", id = "PathwayOraEmaplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot3")
          ),
          box(title = "DotPlot", id = "PathwayOraDotplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot5")
          ),
          box(title = "GoPlot ", id = "PathwayGoplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot7")
          ),
          shinyjs::disable("PathwayRunButton"),
          shinyjs::disable("PathwayAnalysisType"),
          shinyjs::disable("PathwayDatabase")
        ) } 
      else if ("GSEA" %in% input$PathwayAnalysisType) { 
        fluidRow( 
          box(title = "BarPlot", id = "PathwayGseaBarplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("barplot2")
          ),
          box(title = "CnetPlot", id = "PathwayGseaCnetplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot2")
          ),
          box(title = "EmaPlot", id = "PathwayGseaEmaplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot4")
          ),
          box(title = "DotPlot", id = "PathwayGseaDotplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot6")
          ),
          box(title = "GseaPlot ", id = "PathwayGseaplotBox", width = 12, solidHeader = TRUE, status = 'primary',
              plotOutput("plot8")
          ),
          shinyjs::disable("PathwayRunButton"),
          shinyjs::disable("PathwayAnalysisType"),
          shinyjs::disable("PathwayDatabase")
        ) }
    }
  })
  
  observeEvent(input$clearPathway, {
    PathwayAnalysisReady(FALSE)
    
    shinyjs::enable("PathwayRunButton")
    shinyjs::enable("PathwayAnalysisType")
    shinyjs::enable("PathwayDatabase")
  })
  
  #############################################################################################################################################################################################################################################################"
  #############################################################################################################################################################################################################################################################"
  #############################################################################################################################################################################################################################################################" 
  
  
})




###Rajouter la liaison avec l'organisme séléctionné == en fait si j'ai bien compris, l'organisme choisi devra être mis dans les différentes fonctions plus précisément dans l'option organism de chacune. Donc il faut bien cahrger toutes les bases de données.
## De plus dans ces mêmes fonctions, il va falloir choisir la pAdjustMethod qui se fera dans le truc de choix (Bonferonni etc) il faut donc changer les fonctions de sorte à laisser la choix à l'utilisateur et que donc l'analyse soit différente selon le choix.
## A discuter avec les boys : 
#                  Commencer à intégrer les fonctions car le choix de la database, ontologies ou le déplacement des sliders ne fait aucune différences actuellement (normal j'ai testé des figures "classiquess" qui n'utilisent pas ces variables) §§§ Je pense que c'est le plus gros point il faudrait 1 sur Go et 1 sur Pathway 
#                  Bien remplir les organismes que l'on peut choisir 
#                  Rajouter des méthodes d'ajustements ? Peut être même juste 1 ou 2 histoire d'être cohérent
#                  faire les modifs de la page home: limite mettre des logos de l'universté, master etc avec un rapide manuel peut être ? 
#                  dans le about, essayer de "styliser" un peu la page: idée == rajouter nos LinkedIn ? mais ça on le fera en dernier c'est plutôt du bonus 
#                  uniformiser les polices d'écriture, tailles etc 








