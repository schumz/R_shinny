# **General Presentation**

Welcome to **NOM A METTRE,** your comprehensive solution for whole genome data analysis, Gene Ontology (GO) term enrichment, and pathway enrichment exploration. Designed with researchers and biologists in mind, our application facilitates the intricate process of analyzing genomic data, identifying significant GO terms, and uncovering enriched pathways to drive forward scientific discoveries and insights into gene functions, biological processes, and molecular interactions.


## **Core Features:**

- **Whole Genome Data Inspection**: Our application specializes in the visual representation of gene expression data through detailed volcano plots, a powerful tool for quickly identifying genes that are significantly upregulated or downregulated. This feature is pivotal for researchers aiming to pinpoint genes that show differential expression under various conditions or treatments.
- **GO Term Enrichment**: With an integrated Gene Ontology database, our platform allows for the efficient identification of enriched GO terms associated with your gene sets. This feature aids in understanding the functional implications of your genomic data, highlighting biological processes, cellular components, and molecular functions tied to your genes of interest.
- **Pathway Enrichment Analysis**: Dive deeper into the biological pathways with our pathway enrichment feature. By linking genomic data to known pathways, users can uncover significant associations, predict molecular interactions, and gain insights into the underlying mechanisms of diseases or phenotypes.

# **How to Download and Use**

### Step 1: Install Docker

Before you begin, you need to have Docker installed on your machine. If you don't have Docker installed, please visit [Docker's official website](https://www.docker.com/) and follow the instructions to install Docker on your operating system.

### Step 2: Pull the Image

Once Docker is installed, open a terminal and pull the Docker image of the application by running the following command:

```bash
sudo docker pull nutui/myshinyapp:latest
```

### Step 3: Launch the Application
To start the application, run the following command in your terminal:

```bash
sudo docker run --rm -p 3838:3838 nutui/myshinyapp:latest
```

### Step 4: Access the Application
After starting the application, open your web browser and go to the following URL to access the application:

```bash
http://localhost:3838
```

# **Usage**

Getting Started with **NOM A METTRE**. Our platform is designed to work seamlessly with data from **\[analyse à mettre à jour\],** provided in a CSV format. This section guides you through preparing your data file and navigating through the initial steps to analyze your data.

## **Preparing Your CSV File:**

·       **Source Your Data:** Begin by exporting your genomic data from **\[analyse à mettre à jour\]** into a CSV file. This step is crucial for ensuring that the data fed into **\[Nom de l'Application\]** is in the correct format for analysis.

·       **Format Your CSV File**: For a smooth analysis process, your CSV file must include the following essential columns:

\-GeneName: The name of the gene.

\-ID: A unique identifier for the gene.

\-baseMean: The average expression of the gene across samples.

\-log2FC (Log2 Fold Change): The logarithmic scale of fold change in gene expression, indicating upregulation or downregulation.

\-pval (P-value): The statistical significance of the observed change in gene expression.

\-padj (Adjusted P-value): The P-value adjusted for multiple testing corrections.

&nbsp;

Ensuring these columns are present and correctly labeled in your CSV file is essential for the accurate analysis and visualization of your genomic data.

&nbsp;

## **Upload Your Data:** 
With your CSV file prepared, navigate to the "Select a CSV file :" section. Here, you can upload your file directly into the application. Our platform includes verification mechanisms to check the integrity and format of your data, ensuring compatibility with our analysis tools.

&nbsp;

- **Whole Data Inspection:** This section is instrumental for graphically filtering your data and offers the capability to display it in a table. It provides the flexibility to adjust two thresholds through sliders located in the option box:

**\-P-value cutoff**

**\-Log2 Fold Change cutoff**

It's important to note that changes to these sliders affect the plot and the table in real time. The plot uses three colors to indicate genes that are above or below the Log2FC threshold in absolute value, with points falling between these values shown in gray. Additionally, there is an option to download the plot for further analysis or presentation purposes.

If you are interested in a specific gene, the table includes a search function that allows you to quickly locate it. This feature enhances the usability of the application by enabling targeted inspection of genes of interest.

- **Go Term Enrichment :**

In the GO Term Enrichment section of our application, you have several parameters at your disposal to tailor the analysis to your needs. First, you can choose the type of analysis:

·        **ORA (Over-Representation Analysis**): This method helps identify if a set of genes of interest (e.g., upregulated or downregulated genes) is significantly enriched for specific Gene Ontology (GO) terms compared to what would be expected by chance. This involves comparing the frequency of GO terms within your gene set against a reference frequency, usually the entire genome.

·        **GSEA (Gene Set Enrichment Analysis):** Unlike ORA, which analyzes pre-selected gene sets, GSEA assesses whether members of a gene set (e.g., associated with a specific GO term) systematically rank higher or lower in a list of genes sorted by their differential expression. This allows for a more nuanced interpretation of GO term enrichment, considering the entire expression profile.

&nbsp;

You also have the option to select among different ontologies for a more focused analysis:

·        **Biological Process**: Describes the biological processes a gene is involved in, such as DNA replication or immune response.

·        **Molecular Function**: Refers to the molecular activities performed by the gene products, like enzyme activities or binding functions.

·        **Cellular Component**: Indicates the parts of a cell or its extracellular environment where the gene product is active, such as the nucleus, cytoplasm, or cell membrane.

Once you've made your selections among the different types of analysis and ontologies, note that for ORA, you can choose a p-value threshold and a log2fold change threshold to refine your gene set. After setting your parameters, click on "Run Analysis" to generate the corresponding plots showcasing the enrichment results.

&nbsp;

To initiate another type of analysis, simply click on "Clear Results" and make your new selection, allowing for seamless transition between different analytical approaches.

- **Pathway Enrichment :**

In the Pathway Enrichment section of our application, you have several parameters at your disposal to tailor the analysis to your needs. First, you can choose the type of analysis.

·        **ORA** (Over-Representation Analysis): This method helps identify if a set of genes of interest (e.g., upregulated or downregulated genes) is significantly enriched for specific Gene Ontology (GO) terms compared to what would be expected by chance. This involves comparing the frequency of GO terms within your gene set against a reference frequency, usually the entire genome.

·        **GSEA** (Gene Set Enrichment Analysis): Unlike ORA, which analyzes pre-selected gene sets, GSEA assesses whether members of a gene set (e.g., associated with a specific GO term) systematically rank higher or lower in a list of genes sorted by their differential expression. This allows for a more nuanced interpretation of GO term enrichment, considering the entire expression profile.

You also have the option to select patway databses for a more focused analysis:

·        **KEGG** (Kyoto Encyclopedia of Genes and Genomes): A widely-used database that provides a wealth of information on genomes, biological pathways, diseases, drugs, and chemical substances. KEGG is invaluable for understanding high-level functions and utilities of the biological system.

·        **Reactome**: An open-source, curated database of biological pathways. Reactome helps researchers to visualize complex pathways and interpret their data in the context of human biology. It offers detailed descriptions of pathway components and their interactions, making it a powerful tool for pathway-based analysis.

&nbsp;

Choose between ORA or GSEA based on your research needs and select either KEGG or Reactome as your pathway database. For ORA analyses, adjust the p-value and log2 fold change cutoffs to focus on the most relevant pathways for your study. Once your parameters are set, click on "Run Analysis" to generate the enrichment plots and identify significant pathway associations.

&nbsp;

To explore different analyses or adjust your parameters, simply click on "Clear Results" and reconfigure your selections as needed. This flexibility allows for iterative exploration and deep dives into the data, facilitating a thorough understanding of the biological implications of your study.
