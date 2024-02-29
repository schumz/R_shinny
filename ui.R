# Application Shiny v6
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
                    menuItem("About", tabName = "About", icon = icon("th", lib = "glyphicon")),
                    menuItem("Documentation", tabName = "Documentation", icon = icon("book", lib = "font-awesome"))
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
                            div(class = "custom-description",
                                HTML('
<img src="bims_logo.png" style="height:auto; width:10%;" />
<img src="logo_univ.png" style="height:auto; width:30%;" />

<h1>General Presentation</h1>
<p>Welcome to <strong>NOM A METTRE</strong>, your comprehensive solution for whole genome data analysis, Gene Ontology (GO) term enrichment, and pathway enrichment exploration. Designed with researchers and biologists in mind, our application facilitates the intricate process of analyzing genomic data, identifying significant GO terms, and uncovering enriched pathways to drive forward scientific discoveries and insights into gene functions, biological processes, and molecular interactions.</p>

                ')
                            ),
                            # Modification pour séparer la section "Core Features" en trois colonnes
                            fluidRow(
                              tags$h3("Core Features:", style = "text-align: left; margin-left: 15px;"),
                              div(style = "text-align: left; margin-left: 15px; margin-bottom: 20px;", 
                                  tags$p("Find an overview above. Click on the titles to navigate to the corresponding pages.")),
                              column(4,
                                     tags$h4(actionLink("link_whole_genome", "Whole Genome Data Inspection", style = "font-weight: bold; text-decoration: underline; text-align: center; cursor: pointer;")),
                                     div(style = "text-align: left;", 
                                         tags$p("Our application specializes in the visual representation of gene expression data through detailed volcano plots, a powerful tool for quickly identifying genes that are significantly upregulated or downregulated. This feature is pivotal for researchers aiming to pinpoint genes that show differential expression under various conditions or treatments."),
                                         HTML('<img src="volcano_plot.png" style="height:auto; width:90%;" />'),
                                         
                                     )
                              ),
                              column(4,
                                     tags$h4(actionLink("link_go_term", "GO Term Enrichment", 
                                                        style = "font-weight: bold; text-decoration: underline; text-align: center; cursor: pointer;")),
                                     div(style = "text-align: left;", 
                                         tags$p("With an integrated Gene Ontology database, our platform allows for the efficient identification of enriched GO terms associated with your gene sets. This feature aids in understanding the functional implications of your genomic data, highlighting biological processes, cellular components, and molecular functions tied to your genes of interest."),
                                         HTML('<img src="ora_gsea_plot.png" style="height:auto; width:90%; margin-bottom: 20px;" />'), # Ajout de margin-bottom ici
                                         HTML('<img src="process.png" style="height:auto; width:90%;" />')
                                     )
                                     
                              ),
                              column(4,
                                     tags$h4(actionLink("link_pathway_enrichment", "Pathway Enrichment Analysis", style = "font-weight: bold; text-decoration: underline; text-align: center; cursor: pointer;")),
                                     div(style = "text-align: left;", 
                                         tags$p("Dive deeper into the biological pathways with our pathway enrichment feature. By linking genomic data to known pathways, users can uncover significant associations, predict molecular interactions, and gain insights into the underlying mechanisms of diseases or phenotypes."),
                                         HTML('<img src="ora_gsea_plot.png" style="height:auto; width:90%; margin-bottom: 20px;" />'), # Ajout de margin-bottom ici
                                         HTML('<img src="kegg_reactome.png" style="height:auto; width:90%;" />')
                                         
                                     )
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
                            
                            uiOutput("GoTermAnalysisBox") %>% withSpinner(color="#0dc5c1",  type = 6)
                            
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
                            
                            uiOutput("PathwayAnalysisBox")%>% withSpinner(color="blue",  type = 6)
                            
                    ),
                    tabItem(tabName = "About",
                            tags$h2(style = "text-align: center; text-decoration: underline; font-size: 50px;", "About"),
                            br(),
                            br(),
                            div(class = "about_description",
                                HTML('
                                <h2>Project Information and Contributors</h2>
    <p>This project is publicly available for users, researchers, and developers alike on GitHub, providing a comprehensive resource for everything from getting started to diving deep into the application\'s functionalities. The repository includes detailed documentation on the Home section, as well as examples to help users familiarize themselves with the application and its features.</p>

    <h3>Access the Project on GitHub:</h3>
    <p>To explore the project, access code repositories, and view detailed documentation, please visit our GitHub page at the following link: https://github.com/schumz/R_shinny .Here, you will find all the necessary information and resources to get started, understand the application\'s capabilities, and learn how to utilize its features effectively.</p>

    <h3>Contributors:</h3>
    <p>This application was developed by three students from the Master\'s program in Bioinformatics and Modeling Systems (BIMS) at the University of Rouen:</p>
    <ul>
      <li><strong>Mathieu Zallio :</strong> mathieu.zallio@univ-rouen.fr</li>
      <li><strong>Adam Schumasher :</strong> adam.schumasher@univ-rouen.fr</li>
      <li><strong>Baptiste Herlemont :</strong> baptiste.herlemont@univ-rouen.fr</li>
    </ul>
    <p>Their collective effort and expertise have been instrumental in bringing this project to life, showcasing a blend of technical skills, bioinformatics knowledge, and a passion for advancing genomic research.</p>

    <h3>Technical Details:</h3>
    <p>The development of this application was carried out using R version 4.3.2, ensuring compatibility with the latest features and packages available in the R ecosystem. This choice of software version reflects our commitment to utilizing cutting-edge tools and technologies to provide a robust, efficient, and user-friendly application.</p>

    <h3>Acknowledgements:</h3>
    <p>We extend our heartfelt thanks to the University of Rouen, the BIMS program, and all who have supported this project. Your guidance and encouragement have been invaluable in the realization of this endeavor.</p>

    <p>For feedback, contributions, or inquiries, please do not hesitate to reach out to us via our GitHub page. We are continually looking to improve the application and welcome any suggestions or contributions from the community.</p>
  ')
                            )
                    ),
                    tabItem(tabName = "Documentation",
                            tags$h2(style = "text-align: center; text-decoration: underline; font-size: 50px;", "User Manual"),
                            br(),
                            br(),
                            div(class = "about_description",
                                HTML('
                   <h1>Usage</h1>
<p>Getting Started with <strong>NOM A METTRE</strong>. Our platform is designed to work seamlessly with data from <em>[analyse à mettre à jour]</em>, provided in a CSV format. This section guides you through preparing your data file and navigating through the initial steps to analyze your data.</p>

<h3>Preparing Your CSV File: <img src="select_file.png" style="height:auto; width:10%;" /> </h3>
<ul>
    <li><strong>Source Your Data:</strong> Begin by exporting your genomic data from <em>[analyse à mettre à jour]</em> into a CSV file. This step is crucial for ensuring that the data fed into <em>[Nom de l\'Application]</em> is in the correct format for analysis.</li>
  <li><strong>Format Your CSV File:</strong> For a smooth analysis process, your CSV file must include the following essential columns:
  <ul>
  <li>GeneName: The name of the gene.</li>
  <li>ID: A unique identifier for the gene.</li>
  <li>baseMean: The average expression of the gene across samples.</li>
  <li>log2FC (Log2 Fold Change): The logarithmic scale of fold change in gene expression, indicating upregulation or downregulation.</li>
  <li>pval (P-value): The statistical significance of the observed change in gene expression.</li>
  <li>padj (Adjusted P-value): The P-value adjusted for multiple testing corrections.</li>
  </ul>
  </li>
  </ul>
  <p>Ensuring these columns are present and correctly labeled in your CSV file is essential for the accurate analysis and visualization of your genomic data.</p>
  
  <h3>Upload Your Data:</h3>
  <p>With your CSV file prepared, navigate to the "Select a CSV file :" section. Here, you can upload your file directly into the application. Our platform includes verification mechanisms to check the integrity and format of your data, ensuring compatibility with our analysis tools.</p>
  
  <h2>Whole Data Inspection:
<img src="whole_data_inspection.png" style="height:auto; width:10%;" /></h2>
  <p>This section is instrumental for graphically filtering your data and offers the capability to display it in a table. It provides the flexibility to adjust two thresholds through sliders located in the option box:
  <ul>
  <li>P-value cutoff</li>
  <li>Log2 Fold Change cutoff</li>
  </ul>
  </p>
  <p>It\'s important to note that changes to these sliders affect the plot and the table in real time. The plot uses three colors to indicate genes that are above or below the Log2FC threshold in absolute value, with points falling between these values shown in gray. Additionally, there is an option to download the plot for further analysis or presentation purposes.</p>
<p>If you are interested in a specific gene, the table includes a search function that allows you to quickly locate it. This feature enhances the usability of the application by enabling targeted inspection of genes of interest.</p>

<h2>Go Term Enrichment:<img src="go_term_enrichment.png" style="height:auto; width:10%;" /></h2>
<p>In the GO Term Enrichment section of our application, you have several parameters at your disposal to tailor the analysis to your needs. First, you can choose the type of analysis:</p>
<ul>
    <li>ORA (Over-Representation Analysis): This method helps identify if a set of genes of interest (e.g., upregulated or downregulated genes) is significantly enriched for specific Gene Ontology (GO) terms compared to what would be expected by chance. This involves comparing the frequency of GO terms within your gene set against a reference frequency, usually the entire genome.</li>
    <li>GSEA (Gene Set Enrichment Analysis): Unlike ORA, which analyzes pre-selected gene sets, GSEA assesses whether members of a gene set (e.g., associated with a specific GO term) systematically rank higher or lower in a list of genes sorted by their differential expression. This allows for a more nuanced interpretation of GO term enrichment, considering the entire expression profile.</li>
</ul>
<p>You also have the option to select among different ontologies for a more focused analysis:</p>
<ul>
    <li>Biological Process: Describes the biological processes a gene is involved in, such as DNA replication or immune response.</li>
    <li>Molecular Function: Refers to the molecular activities performed by the gene products, like enzyme activities or binding functions.</li>
    <li>Cellular Component: Indicates the parts of a cell or its extracellular environment where the gene product is active, such as the nucleus, cytoplasm, or cell membrane.</li>
</ul>
<p>Once you\'ve made your selections among the different types of analysis and ontologies, note that for ORA, you can choose a p-value threshold and a log2fold change threshold to refine your gene set. After setting your parameters, click on "Run Analysis" to generate the corresponding plots showcasing the enrichment results.</p>
  <p>To initiate another type of analysis, simply click on "Clear Results" and make your new selection, allowing for seamless transition between different analytical approaches.</p>
  
  <h2>Pathway Enrichment:<img src="pathway_enrichment.png" style="height:auto; width:10%;" /></h2>
<p>In the Pathway Enrichment section of our application, you have several parameters at your disposal to tailor the analysis to your needs. First, you can choose the type of analysis:</p>
<ul>
    <li>ORA (Over-Representation Analysis): This method helps identify if a set of genes of interest (e.g., upregulated or downregulated genes) is significantly enriched for specific Gene Ontology (GO) terms compared to what would be expected by chance. This involves comparing the frequency of GO terms within your gene set against a reference frequency, usually the entire genome.</li>
    <li>GSEA (Gene Set Enrichment Analysis): Unlike ORA, which analyzes pre-selected gene sets, GSEA assesses whether members of a gene set (e.g., associated with a specific GO term) systematically rank higher or lower in a list of genes sorted by their differential expression. This allows for a more nuanced interpretation of GO term enrichment, considering the entire expression profile.</li>
</ul>  
<p>You also have the option to select patway databses for a more focused analysis:</p>
<ul>
<li>KEGG (Kyoto Encyclopedia of Genes and Genomes): A widely-used database that provides a wealth of information on genomes, biological pathways, diseases, drugs, and chemical substances. KEGG is invaluable for understanding high-level functions and utilities of the biological system.</li>
<li>Reactome: An open-source, curated database of biological pathways. Reactome helps researchers to visualize complex pathways and interpret their data in the context of human biology. It offers detailed descriptions of pathway components and their interactions, making it a powerful tool for pathway-based analysis.</li>
</ul>
<p> Choose between ORA or GSEA based on your research needs and select either KEGG or Reactome as your pathway database. For ORA analyses, adjust the p-value and log2 fold change cutoffs to focus on the most relevant pathways for your study. Once your parameters are set, click on "Run Analysis" to generate the enrichment plots and identify significant pathway associations.

To explore different analyses or adjust your parameters, simply click on "Clear Results" and reconfigure your selections as needed. This flexibility allows for iterative exploration and deep dives into the data, facilitating a thorough understanding of the biological implications of your study.</p>
                   
<h1>Script Structure</h1>    

<br>

<p>This application was developed on R version 4.2.2 and utilizes the following libraries:</p>
    <ul>
        <li>shiny</li>
        <li>shinydashboard</li>
        <li>plotly</li>
        <li>shinyalert</li>
        <li>htmlwidgets</li>
        <li>ggplot2</li>
        <li>shinyjs</li>
    </ul>
    <p>It is structured in the form of two scripts, <code>UI.R</code> and <code>Server.R</code>. A <code>www</code> directory must be a subdirectory of the working environment to properly use images.</p>
    <p>Regarding the script itself, we have utilized functions as much as possible to maximize reproducibility.</p>
<h1>Toy Examples</h1>
Please refer to the GitHub repository
https://github.com/schumz/R_shinny
to use the example.csv file in the application
<br>
')
                  )
                )
                  )
                )
  )
)
