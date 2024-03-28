# Application Shiny v7
# Créateurs: Baptiste Herlemont
#            Adam Schumacher
#            Mathieu Zallio 
# Email: baptiste.herlemont@univ-rouen.fr
#        adam.schumacher@univ-rouen.fr
#        mathieu.zallio@univ-rouen.fr
# Université affiliée : Université de Rouen Normandie 

library(shiny)
library(shinydashboard)
library(shinycssloaders)
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
library(ReactomePA)
library(stats)

conflicts(detail = TRUE)



#################################################################################################################################################################################################################"

ontologies_associations <- list(
  "Biological process" = "BP",
  "Molecular function" = "MF",
  "Cellular component" = "CC"
)

methodes_ajustement <- list(
  "Holm" = "holm",
  "Hochberg" = "hochberg",
  "Hommel" = "hommel",
  "Bonferroni" = "bonferroni",
  "Benjamini-Hochberg" = "BH",
  "Benjamini-Yekutieli" = "BY",
  "Contrôle du FDR" = "fdr",
  "Aucun ajustement" = "none"
)

#################################################################################################################################################################################################################"
############################################################## PARTIE CREATION DE FONCTIONS UTILISABLES #########################################################################################################"
#################################################################################################################################################################################################################"

        #################################
####### FONCTIONS POUR GO TERM ENRICHMENT ########
        #################################

######################################################################
###################### FONCTION GO ORA ANALYSIS ######################
######################################################################
GO_ora_analysis <- function(filtered_data, ontologies_terms, pajust_method, org_db = org.Mm.eg.db, simplifyCutoff = 0.7) {
  
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
      go_results <- enrichGO(gene = entrez_ids, OrgDb = org_db, keyType = "ENTREZID", ont = ont, pAdjustMethod=pajust_method)
      go_results_simplified <- simplify(go_results, cutoff = simplifyCutoff)
      
      # Ajouter les résultats simplifiés à la liste
      simplified_results_list[[term]] <- go_results_simplified
    }, error = function(e) {
      cat("Error in ORA for", ont, ": ", e$message, "\n")
    })
  }
  
  # Retourner la liste des résultats simplifiés
  return(simplified_results_list)
}


######################################################################
###################### FONCTION GO GSEA ANALYSIS #####################
######################################################################
GO_gsea_analysis <- function(data, ontologies_terms, pajust_method,  org_db = org.Mm.eg.db, eps = 1e-300,simplifyCutoff = 0.7) {
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
      gse_go_result <- gseGO(geneList = gene_GSEA_sorted, ont = ont, keyType = "ENTREZID", OrgDb = org_db, eps = eps, pAdjustMethod = pajust_method)
      gse_go_results_simplified <- simplify(gse_go_result, cutoff = simplifyCutoff)
      
      simplified_results_list[[term]] <- gse_go_results_simplified
    }, error = function(e) {
      cat("Error in GSEA for", ont, ": ", e$message, "\n")
    })
  }
  
  return(simplified_results_list)
}

######################################################################
################### FONCTION GO APPARITION DE BOX ####################
######################################################################

# Fonction auxiliaire pour générer une box DataTable
generateDataTableBox <- function(ontology, analysisType) {
  if(analysisType == "both") {
    titleORA <- paste("Datatable after Analysis ORA")
    outputIdORA <- paste("DatatableORA_", ontology, sep="")
    titleGSEA <- paste("Datatable after Analysis GSEA")
    outputIdGSEA <- paste("DatatableGSEA_", ontology, sep="")
    
    fluidRow(box(title = titleORA, status = "primary", solidHeader = TRUE, width = 12,
                        dataTableOutput(outputId = outputIdORA)) ,
            box(title = titleGSEA, status = "primary", solidHeader = TRUE, width = 12,
                        dataTableOutput(outputId = outputIdGSEA)) 
            )
    
  } else {
  title <- paste("Datatable after Analysis", analysisType)
  outputId <- paste("Datatable", analysisType, "_", ontology, sep="")
  
  fluidRow( box(title = title, status = "primary", solidHeader = TRUE, width = 12,
                    dataTableOutput(outputId = outputId)) 
      )
}}

# Fonction auxiliaire pour générer une box de Plot
generatePlotBox <- function(ontology, analysisType, plotType) {
  title <- paste(plotType)
  outputIdORA <- paste(plotType, "ORA_", ontology, sep="")
  outputIdGSEA <- paste(plotType, "GSEA_", ontology, sep="")
  
  if(analysisType == "both") {
    content <- fluidRow(
      column(6, plotOutput(outputId = outputIdORA)),
      column(6, plotOutput(outputId = outputIdGSEA))
    )
  } else {
    content <- fluidRow(plotOutput(outputId = if(analysisType == "ORA") outputIdORA else outputIdGSEA))
  }
  
  fluidRow(box(title = title, id = paste(plotType, "Box", ontology, sep="_"), width = 12,
      solidHeader = TRUE, status = 'primary', content))
}

# Fonction auxiliaire pour désactiver des éléments de l'UI
disableUiElements <- function() {
  list(shinyjs::disable("GoTermRunButton"), 
       shinyjs::disable("GoTermAnalysisType"), 
       shinyjs::disable("OntologiesSelect"))
}


generateGoTermAnalysisBoxes <- function(ontologiesSelected, analysisTypes) {
  GOallUiElements <- list()
  
  for(ontology in ontologiesSelected) {
    # Ajoutez un titre pour chaque ontologie
    GOallUiElements[[length(GOallUiElements) + 1]] <- tags$h3(ontology, style = "margin-top: 20px;")
    
    # Déterminez les box à ajouter en fonction des types d'analyse sélectionnés
    if("ORA" %in% analysisTypes && "GSEA" %in% analysisTypes) {
      analysisType <- "both"
    } else if("ORA" %in% analysisTypes) {
      analysisType <- "ORA"
    } else if("GSEA" %in% analysisTypes) {
      analysisType <- "GSEA"
    } else {
      next
    }
    
    plotTypes <- c("Barplot", "Cnetplot", "Emaplot", "Dotplot", if(analysisType == "both") "Goplot & Gseaplot" else if(analysisType == "ORA") "Goplot" else "Gseaplot")
    
    # Générer les boxes DataTable
    GOallUiElements[[length(GOallUiElements) + 1]] <- generateDataTableBox(ontology, analysisType)
    
    # Générer les boxes de plot pour chaque type de plot
    for(plotType in plotTypes) {
      GOallUiElements[[length(GOallUiElements) + 1]] <- generatePlotBox(ontology, analysisType, plotType)
    }
  }
   # Désactiver les éléments UI après exécution
  GOallUiElements <- c(GOallUiElements, disableUiElements())
  
  # Retournez tous les éléments UI en utilisant tagList pour les grouper sans ajouter de fluidRow supplémentaire
  do.call(tagList, GOallUiElements)
}





        #################################
####### FONCTIONS POUR PATHWAY ENRICHMENT ########
        #################################


######################################################################
##################### FONCTION KEGG ORA ANALYSIS #####################
######################################################################
kegg_ora_analysis <- function(filtered_data, pajust_method, organism = "mmu", key_type = "ENSEMBL") {
  
  gene_ORA <- as.character(filtered_data$ID)
  
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = gene_ORA,
                       column = "ENTREZID",
                       keytype = key_type,
                       multiVals = "first")
  
  
  tryCatch({
    kegg_results <- enrichKEGG(gene = entrez_ids, organism = organism, keyType = "ncbi-geneid",pAdjustMethod=pajust_method)
  }, error = function(e) {
    cat("Error in KEGG ORA:", e$message, "\n")
  })
  return(kegg_results)
}


######################################################################
#################### FONCTION KEGG GSEA ANALYSIS #####################
######################################################################
kegg_gsea_analysis <- function(data, pajust_method, organism = "mmu", key_type = "ENSEMBL", eps = 1e-300, pAdjustMethod = "BH") {
  
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
    gse_kegg_result <- gseKEGG(geneList = entrez_geneList_sorted, organism = organism, keyType = "ncbi-geneid", eps = eps, pAdjustMethod = pajust_method)
  }, error = function(e) {
    cat("Error in KEGG GSEA:", e$message, "\n")
  })
}

######################################################################
###################### FONCTION REACTOME GSEA ANALYSIS ###############  
######################################################################


Reactome_GSEA_Analysis <- function(data, pajust_method, pvalueCutoff = 0.05, organism = "mouse") {
  
  # Lire, nettoyer, transformer et trier les données en une seule séquence d'opérations
  data_sorted <- data%>%
    filter(!duplicated(ID)) %>% #enleve les duplicats si il y en a
    filter(!is.na(log2FC) & !is.na(ID) & !is.na(padj)) %>% #enleve les lignes avec des NaN
    #mutate(ranking = -log10(padj) * log2FC) %>% # Ajoute une colonne ranking
    arrange(desc(log2FC)) #tri par ordré décroissant par rapport au log2FC
  
  # Extraction des identifiants Ensembl
  ensembl_ids_GSEA <- data_sorted$ID
  
  #Conversion en EntrezID
  entrez_ids_GSEA <- mapIds(
    org.Mm.eg.db, 
    keys = ensembl_ids_GSEA, 
    keytype = "ENSEMBL", 
    column = "ENTREZID")
  
  # Ajout de la colonne EntrezID aux données triées
  data_sorted <- data_sorted %>%
    mutate(EntrezID = entrez_ids_GSEA)
  
  # Préparation de la liste de gènes pour GSEA
  gene_list <- setNames(data_sorted$log2FC, data_sorted$EntrezID)
  
  #Analyse d'enrichissement de voies géniques par GSEA
  reactome_results_GSEA <- gsePathway(geneList = gene_list, organism = organism, 
                                      exponent = 1, 
                                      nPerm = 1000,
                                      pvalueCutoff = pvalueCutoff, 
                                      pAdjustMethod = pajust_method, 
                                      by = "fgsea")
  
  return(reactome_results_GSEA)
}


######################################################################
###################### FONCTION REACTOME ORA ANALYSIS ###############  
######################################################################

Reactome_ORA_Analysis <- function(filtered_data, pajust_method, pvalueCutoff = 0.05, organism = "mouse") {
  # Extraction des identifiants Ensembl
  ensembl_ids_ORA <- filtered_data$ID
  
  #Conversion en EntrezID
  entrez_ids_ORA <- mapIds(
    org.Mm.eg.db, 
    keys = ensembl_ids_ORA, 
    keytype = "ENSEMBL", 
    column = "ENTREZID")
  
  # Ajout de la colonne EntrezID aux données filtrées
  filtered_data <- mutate(filtered_data, EntrezID = entrez_ids_ORA)
  
  # Extraction des noms de gènes
  gene_ORA <- as.character(filtered_data$EntrezID)
  
  # Analyse d'enrichissement de voies (ORA)
  reactome_results_ORA <- enrichPathway(gene = gene_ORA, 
                                        organism = organism,
                                        pvalueCutoff = pvalueCutoff,
                                        pAdjustMethod = pajust_method,
                                        qvalueCutoff = 0.05,
                                        universe = NULL, 
                                        minGSSize = 10, 
                                        maxGSSize = 500,
                                        readable = FALSE)
  
  return(reactome_results_ORA)
}





######################################################################
################# FONCTION PATHWAY APPARITION DE BOX #################
######################################################################

generatePathwayAnalysisBoxes <- function(database_selected, PathanalysisTypes) {
  # Initialisez une liste vide pour contenir tous les éléments UI
  PathAllUiElements <- list()
  
  # Parcourez chaque ontologie sélectionnée
  for(database in database_selected) {
    
    # Ajoutez un titre pour chaque ontologie
    PathAllUiElements[[length(PathAllUiElements) + 1]] <- tags$h3(database, style = "margin-top: 20px;")
    
    # Déterminez les box à ajouter en fonction des types d'analyse sélectionnés
    if("ORA" %in% PathanalysisTypes && "GSEA" %in% PathanalysisTypes) {
      # Ajoutez les box pour ORA et GSEA sans les envelopper dans un fluidRow ici
      PathAllUiElements[[length(PathAllUiElements) + 1]] <- fluidRow(box(title = paste("BarPlot"), 
                                                                         id = paste("PathwayOraGseaBarplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         fluidRow(
                                                                           column(6, plotOutput("barplot1")), # Pour le premier plot
                                                                           column(6, plotOutput("barplot2"))  # Pour le deuxième plot
                                                                         )),
                                                                     box(title = paste("CnetPlot"), 
                                                                         id = paste("PathwayOraGseaCnetplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         fluidRow(
                                                                           column(6, plotOutput("plot1")), # Pour le premier plot
                                                                           column(6, plotOutput("plot2"))  # Pour le deuxième plot
                                                                         )),
                                                                     box(title = paste("EmaPlot"), 
                                                                         id = paste("PathwayOraGseaEmaplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         fluidRow(
                                                                           column(6, plotOutput("plot3")), # Pour le premier plot
                                                                           column(6, plotOutput("plot4"))  # Pour le deuxième plot
                                                                         )),
                                                                     box(title = paste("DotPlot"), 
                                                                         id = paste("PathwayOraGseaDotplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         fluidRow(
                                                                           column(6, plotOutput("plot5")), # Pour le premier plot
                                                                           column(6, plotOutput("plot6"))  # Pour le deuxième plot
                                                                         )),
                                                                     box(title = paste("GoPlot & GseaPlot"), 
                                                                         id = paste("PathwayGoGseaplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         fluidRow(
                                                                           column(6, plotOutput("plot7")), # Pour le premier plot
                                                                           column(6, plotOutput("plot8"))  # Pour le deuxième plot
                                                                         ))
      )
      
    } else if("ORA" %in% PathanalysisTypes) {
      # Ajoutez les box pour ORA seulement
      PathAllUiElements[[length(PathAllUiElements) + 1]] <- fluidRow(box(title = paste("BarPlot"), 
                                                                         id = paste("PathwayOraBarplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         plotOutput("barplot2")),
                                                                     box(title = paste("CnetPlot"), 
                                                                         id = paste("PathwayOraCnetplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         plotOutput("plot2")),
                                                                     box(title = paste("EmaPlot"), 
                                                                         id = paste("PathwayOraEmaplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         plotOutput("plot4")),
                                                                     box(title = paste("DotPlot"), 
                                                                         id = paste("PathwayOraDotplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         plotOutput("plot6")),
                                                                     box(title = paste("GoPlot "), 
                                                                         id = paste("PathwayGoplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         plotOutput("plot8"))
      )
    } else if("GSEA" %in% PathanalysisTypes) {
      # Ajoutez les box pour GSEA seulement
      PathAllUiElements[[length(PathAllUiElements) + 1]] <- fluidRow(box(title = paste("BarPlot"), 
                                                                         id = paste("PathwayGseaBarplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         plotOutput("barplot1")),
                                                                     box(title = paste("CnetPlot"), 
                                                                         id = paste("PathwayGseaCnetplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         plotOutput("plot1")),
                                                                     box(title = paste("EmaPlot"), 
                                                                         id = paste("PathwayGseaEmaplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         plotOutput("plot3")),
                                                                     box(title = paste("DotPlot"), 
                                                                         id = paste("PathwayGseaDotplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         plotOutput("plot5")),
                                                                     box(title = paste("GseaPlot"), 
                                                                         id = paste("PathwayGseaplotBox", database, sep="_"), 
                                                                         width = 12, solidHeader = TRUE, status = 'primary',
                                                                         plotOutput("plot7"))
      )
    }
    
    # Notez que chaque box est ajoutée individuellement à GOallUiElements sans être enveloppée dans un fluidRow ici
  }
  
  # Désactivez les boutons et sélections après l'exécution
  PathAllUiElements <- c(PathAllUiElements, 
                         shinyjs::disable("PathwayRunButton"), 
                         shinyjs::disable("PathwayDatabase"), 
                         shinyjs::disable("PathwayAnalysisType"))
  
  # Retournez tous les éléments UI en utilisant tagList pour les grouper sans ajouter de fluidRow supplémentaire
  do.call(tagList, PathAllUiElements)
}


#################################################################################################################################################################################################################"
#################################################################################################################################################################################################################"
#################################################################################################################################################################################################################"






shinyServer(function(input, output, session) {
  raw_data <- reactiveVal()
  #############################################################################################################################################################################################################################################################"
  ######################################################## PARTIE WHOLE DATA INSPECTION: CHECK DE L'INPUT + VOLCANOPLOT + DOWNLOAD + ACTIVATION DES BOUTONS RUN DE GO ET PATHWAY ##############################################################################"
  #############################################################################################################################################################################################################################################################"

  
  
  
  
  observeEvent(input$link_whole_genome, {
    updateTabsetPanel(session, "WholeDataInspect")
  })
  
  
  observeEvent(input$link_go_term, {
    updateTabItems(session, "GoTermEnrichment")
  })
  
  observeEvent(input$link_pathway_enrichment, {
    updateTabItems(session, "PathwayEnrichment")
  }) 
  
  
  
  
  
    
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
      Ontologies_Selected <- input$OntologiesSelect
      
      GO_seuil_pval <- input$GoTermSliderP_val
      GO_seuil_fc <- input$GoTermSliderFoldChange
      GO_ora_option <- input$GoTermDegOptions
      GO_ajust_method <- methodes_ajustement[[input$GoTermAdjustmentMethod]]

      
      
      
      if ("ORA" %in% GO_choix_analyse) {
        print("c'est good")
        data_test <- raw_data() 
        filtered_data_ora <- subset(data_test,padj <= GO_seuil_pval & abs(log2FC) >= GO_seuil_fc)
        if (GO_ora_option == "Over expressed DEG") {
          # Gardez uniquement les lignes avec log2FC > 0 pour les gènes sur-exprimés
          filtered_data_ora <- filtered_data_ora[filtered_data_ora$log2FC > 0, ]
        } else if (GO_ora_option == "Under expressed DEG") {
          # Gardez uniquement les lignes avec log2FC < 0 pour les gènes sous-exprimés
          filtered_data_ora <- filtered_data_ora[filtered_data_ora$log2FC < 0, ]
        }
        res_GO_ora_globaux = GO_ora_analysis(filtered_data_ora, Ontologies_Selected, pajust_method=GO_ajust_method ,org_db = org.Mm.eg.db,simplifyCutoff = 0.7)
        print(res_GO_ora_globaux)
        
        ##### AJOUT DU 19 MARS #########
        
        for (ontology in Ontologies_Selected) {
          local({
            myOntology <- ontology
            datatableIdORA <- paste("DatatableORA_", myOntology, sep="")
            plotIdBarplotORA <- paste("BarplotORA_", myOntology, sep="")
            plotIdCnetplotORA <- paste("CnetplotORA_", myOntology, sep="")
            plotIdEmaplotORA <- paste("EmaplotORA_", myOntology, sep="")
            plotIdDotplotORA <- paste("DotplotORA_", myOntology, sep="")
            
            print(paste("Génération de Tableau:", datatableIdORA))
            print(paste("Création du plot:", plotIdBarplotORA)) 
            print(paste("Création du plot:", plotIdCnetplotORA))
            print(paste("Création du plot:", plotIdEmaplotORA))
            print(paste("Création du plot:", plotIdDotplotORA))
            
            
            output[[datatableIdORA]] <- renderDataTable({
              res_GO_ora_globaux[[myOntology]]@result
            })
            output[[plotIdBarplotORA]] <- renderPlot({
              barplot(res_GO_ora_globaux[[myOntology]], showCategory=20)
            })
            output[[plotIdCnetplotORA]] <- renderPlot({
              #ora_results_readable <- setReadable(res_GO_ora_globaux[[myOntology]], org.Mm.eg.db)
              cnetplot(res_GO_ora_globaux[[myOntology]])
            })
            output[[plotIdEmaplotORA]] <- renderPlot({
              ora_results_sim <- pairwise_termsim(res_GO_ora_globaux[[myOntology]])
              emapplot(ora_results_sim)
            })
            output[[plotIdDotplotORA]] <- renderPlot({
              dotplot(res_GO_ora_globaux[[myOntology]], showCategory=30) #good
            })
            
            
          })
        }
        ###############################
      }
      
      if ("GSEA" %in% GO_choix_analyse) {
        
        data_test <- raw_data() 
        res_GO_gsea_globaux = GO_gsea_analysis(data_test, Ontologies_Selected, pajust_method=GO_ajust_method , org_db = org.Mm.eg.db,simplifyCutoff = 0.7)
        print("c'est good pour la GSEA")
        print(res_GO_gsea_globaux)
           
        
        ##### AJOUT DU 19 MARS #########
        
        for (ontology in Ontologies_Selected) {
          local({
            myOntology <- ontology
            datatableIdGSEA <- paste("DatatableGSEA_", myOntology, sep="")
            plotIdBarplotGSEA <- paste("BarplotGSEA_", myOntology, sep="")
            plotIdCnetplotGSEA <- paste("CnetplotGSEA_", myOntology, sep="")
            plotIdEmaplotGSEA <- paste("EmaplotGSEA_", myOntology, sep="")
            plotIdDotplotGSEA <- paste("DotplotGSEA_", myOntology, sep="")
            plotIdGseaplot <- paste("Gseaplot_", myOntology, sep="")
            
            print(paste("Génération de Tableau:", datatableIdGSEA))
            print(paste("Création du plot:", plotIdBarplotGSEA)) 
            print(paste("Création du plot:", plotIdCnetplotGSEA))
            print(paste("Création du plot:", plotIdEmaplotGSEA))
            print(paste("Création du plot:", plotIdDotplotGSEA))
            print(paste("Création du plot:", plotIdGseaplot))
            
            
            output[[datatableIdGSEA]] <- renderDataTable({
              
              res_GO_gsea_globaux[[myOntology]]@result %>% select(-core_enrichment)
            })
            output[[plotIdBarplotGSEA]] <- renderPlot({
              barplot(res_GO_gsea_globaux[[myOntology]], showCategory=20)
            })
            output[[plotIdCnetplotGSEA]] <- renderPlot({
              #ora_results_readable <- setReadable(res_GO_gsea_globaux[[myOntology]], org.Mm.eg.db)
              cnetplot(res_GO_gsea_globaux[[myOntology]])
            })
            output[[plotIdEmaplotGSEA]] <- renderPlot({
              gsea_results_sim <- pairwise_termsim(res_GO_gsea_globaux[[myOntology]])
              emapplot(gsea_results_sim)
            })
            output[[plotIdDotplotGSEA]] <- renderPlot({
              dotplot(res_GO_gsea_globaux[[myOntology]], showCategory=30) #good
            })
            
            output[[plotIdGseaplot]] <- renderPlot({
                gseaplot(res_GO_gsea_globaux[[myOntology]], geneSetID = 2)
            })
            
            
          })
        }
        ###############################
      }
      
      generateGoTermAnalysisBoxes(Ontologies_Selected, GO_choix_analyse)
      
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
      
      
      Pathways_choix_analyse <- input$PathwayAnalysisType
      Databases_Select <- input$PathwayDatabase
      
      Pathways_seuil_pval <- input$PathwaySliderP_val
      Pathways_seuil_fc <- input$PathwaySliderFoldChange
      Pathways_ora_option <- input$PathwayDegOptions
      Pathways_ajust_method <- methodes_ajustement[[input$PathwayAdjustmentMethod]]
      
      
      if ("ORA" %in% Pathways_choix_analyse) {
        data_test <- raw_data() 
        filtered_data_ora_pathways <- subset(data_test,padj <= Pathways_seuil_pval & abs(log2FC) >= Pathways_seuil_fc)
        if (Pathways_ora_option == "Over expressed DEG") {
          # Gardez uniquement les lignes avec log2FC > 0 pour les gènes sur-exprimés
          filtered_data_ora_pathways <- filtered_data_ora_pathways[filtered_data_ora_pathways$log2FC > 0, ]
        } else if (Pathways_ora_option == "Under expressed DEG") {
          # Gardez uniquement les lignes avec log2FC < 0 pour les gènes sous-exprimés
          filtered_data_ora_pathways <- filtered_data_ora_pathways[filtered_data_ora_pathways$log2FC < 0, ]
        }
        if ("KEGG" %in% Databases_Select){
          print("ora kegg")
          res_kegg_ora = kegg_ora_analysis(filtered_data_ora_pathways, pajust_method = Pathways_ajust_method, organism = "mmu", key_type = "ENSEMBL")
          print(res_kegg_ora)}
        
        if ("Reactome" %in% Databases_Select){
          print("ora reactome")
          res_reactome_ora = Reactome_ORA_Analysis(filtered_data_ora_pathways, pajust_method = Pathways_ajust_method, pvalueCutoff = 0.05, organism = "mouse")
          print(res_reactome_ora)}
        
      }
      
      if ("GSEA" %in% Pathways_choix_analyse) {
        
        data_test <- raw_data() 
        if ("KEGG" %in% Databases_Select){
          print("gsea kegg")
          res_kegg_gsea = kegg_gsea_analysis(data_test, pajust_method = Pathways_ajust_method  ,organism = "mmu", key_type = "ENSEMBL", eps = 1e-300, pAdjustMethod = "BH")}
        print(res_kegg_gsea)
        if ("Reactome" %in% Databases_Select){
          print("gsea reactome")
          res_reactome_gsea =  Reactome_GSEA_Analysis(data_test, pajust_method = Pathways_ajust_method, pvalueCutoff = 0.05, organism = "mouse")
          print(res_reactome_gsea)}
      }
      
      generatePathwayAnalysisBoxes(Databases_Select, Pathways_choix_analyse)
      
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