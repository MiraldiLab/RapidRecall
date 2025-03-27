library(shiny)
library(Seurat)
library(ggplot2)
library(InteractiveComplexHeatmap)
library(shinyjs)
library(profvis)
profvis({
  

options(shiny.maxRequestSize = 1024 * 1024^2)
source("PlotGene2.R")
source("plotHeatmap.R")

ui <- fluidPage(
  titlePanel("Interactive Gene Plot"),
  tabsetPanel(
    tabPanel("Gene Plot",
             sidebarLayout(
               sidebarPanel(
                 fileInput("counts_file", "Upload Counts File", accept = ".tsv"),
                 fileInput("meta_file", "Upload Meta File", accept = ".tsv"),
                 textInput("gene", "Enter Gene Name"),
                 actionButton("plot", "Plot Gene")
               ),
               mainPanel(
                 plotOutput("genePlot")
               )
             )
    ),
    tabPanel("UMAP Plot",
             sidebarLayout(
               sidebarPanel(
                 fileInput("seurat_file", "Upload Seurat RDS File", accept = ".rds"),
                 uiOutput("reduction_selector"),
                 textInput("feature_gene", "Enter Gene for FeaturePlot")
               ),
               mainPanel(
                 plotOutput("umapPlot"),
                 plotOutput("featurePlot")
               )
             )
    ),
    tabPanel("Heatmap",
             sidebarLayout(
               sidebarPanel(
                 fileInput("file1", "Upload First TSV File (Counts)", accept = ".txt"),
                 fileInput("file2", "Upload Second TSV File (Meta)", accept = ".txt"),
                 checkboxGroupInput("celltypes", "Select Cell Types",
                                    choices = c("Naive", "TCM", "Th1", "Th2", "Th17", "MHCII", "CTL", "Treg", "TCM/TEM"),
                                    selected = c("Naive", "TCM", "Th1", "Th2", "Th17", "MHCII", "CTL", "Treg", "TCM/TEM")
                 ),
                 textInput("genes", "Enter Genes (comma-separated)", 
                           value = "IFNG,IL4,IL5,IL13,IL17A,IL17F,IL21"),
                 actionButton("plot_heatmap", "Plot Heatmap")
               ),
               mainPanel(
                 #plotOutput("heatmapPlot")
                 InteractiveComplexHeatmapOutput(heatmap_id = "heatmap2",height1 = 250, width1 = 800)
               )
             )
    )
  )
)

server <- function(input, output, session) {
  counts <- reactive({
    req(input$counts_file)
    read.table(input$counts_file$datapath, row.names = 1, header = TRUE, sep = "\t")
  })
  
  meta <- reactive({
    req(input$meta_file)
    read.table(input$meta_file$datapath, header = TRUE, sep = "\t")
  })
  
  output$genePlot <- renderPlot({
    req(input$plot)
    counts_data <- counts()
    meta_data <- meta()
    gene <- input$gene
    
    if (!gene %in% rownames(counts_data)) {
      return()
    }
    
    celltypes <- c("Naive","TCM","Th1","Th2","Th17","TCM/TEM","MHCII","CTL","Treg")
    cols <- c("#E41A1C","#4A72A6","#48A462","#7E6E85","#D16948","#E1C62F","#B75F49","#EC83BA","#999999")
    outdir <- tempdir()
    
    plotGene(counts_data, meta_data, gene, celltypes, cols, outdir)
    
    # Directly generate and display the plot
    plotGene(counts_data, meta_data, gene, celltypes, cols, NULL)
  })
  
  seurat_obj <- reactive({
    req(input$seurat_file)
    obj <- readRDS(input$seurat_file$datapath)
    print("Seurat object loaded")
    print(names(obj@reductions))
    print(head(obj@meta.data))
    return(obj)
  })
  
  output$reduction_selector <- renderUI({
    req(seurat_obj())
    reductions <- names(seurat_obj()@reductions)
    selectInput("selected_reduction", "Select Dimension Reduction", choices = reductions, selected = "wnn.umap")
  })
  
  output$umapPlot <- renderPlot({
    req(seurat_obj(), input$selected_reduction)
    
    # Check if the selected reduction and metadata column exist
    if (!(input$selected_reduction %in% names(seurat_obj()@reductions))) {
      stop("Selected reduction not found in the Seurat object.")
    }
    
    if (!("celltype2" %in% colnames(seurat_obj()@meta.data))) {
      stop("Metadata column 'celltype2' not found in the Seurat object.")
    }
    
    # Plot the dimension reduction
    DimPlot(seurat_obj(), reduction = input$selected_reduction, group.by = "celltype2")
  })
  
  output$featurePlot <- renderPlot({
    req(seurat_obj(), input$feature_gene)
    
    # Check if the gene exists in the Seurat object
    if (input$feature_gene %in% rownames(seurat_obj())) {
      FeaturePlot(seurat_obj(), features = input$feature_gene, slot = "counts", max.cutoff = 10)
    } else {
      plot.new()
      text(0.5, 0.5, "Gene not found in the dataset", cex = 1.5)
    }
  })
  
  # Reactive expressions for the third tab
  counts_file1 <- reactive({
    req(input$file1)
    read.table(input$file1$datapath, row.names = 1, header = TRUE, sep = "\t")
  })
  
  meta_file2 <- reactive({
    req(input$file2)
    read.table(input$file2$datapath, header = TRUE, sep = "\t")
  })
  
  observeEvent(input$plot_heatmap,{
    #req(input$plot_heatmap)
    shinyjs::disable("plot_heatmap")
    counts_data <- counts_file1()
    meta_data <- meta_file2()
    selected_celltypes <- input$celltypes
    genes <- unlist(strsplit(input$genes, ",\\s*"))
    
    # Call the plotHeatmap function
    ht <- plotHeatmap(counts = counts_data, meta = meta_data, genes = genes, celltypes = selected_celltypes)
    HT <- draw(ht)
    makeInteractiveComplexHeatmap(input, output, session, HT, heatmap_id = "heatmap2")
    shinyjs::enable("plot_heatmap")
  })
}

#shinyApp(ui, server)
})
shinyApp(ui, server)
