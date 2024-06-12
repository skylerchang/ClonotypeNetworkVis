
## build ui.R -----------------------------------
## 1. header -------------------------------
header <- 
  dashboardHeader( title = HTML("Clonotype Network Visualization Analysis"), 
                   disable = FALSE, 
                   titleWidth  = 400
  )



## 2. siderbar ------------------------------
siderbar <- 
  dashboardSidebar( 
    width = 250,
    sidebarMenu(
      id = 'sidebar',
      style = "position: relative; overflow: visible;",
      menuItem("Data Input and Processing", icon = icon("database"),
               menuSubItem("Input Data", tabName = "tab1"),
               menuSubItem("Process Data", tabName = "tab2"),
               menuSubItem("Integrate Current Data", tabName = "tab3"),
               menuSubItem("Run Anchor Clustering", tabName = "tab4"),
               menuSubItem("Read Anchor Clustering Results", tabName = "tab5")
      ),
      menuItem("Statistical Analysis", icon = icon("chart-bar"),
               menuSubItem("Clonotype Distribution", tabName = "tab6"),
               menuSubItem("Clonotype Overlap", tabName = "tab7"),
               menuSubItem("Circos Plot", tabName = "tab8")
               
      ),
      menuItem("Network Analysis", icon = icon("code-branch"),
               menuSubItem("Scatter Plot of Network Parameters", tabName = "tab9"),
               menuSubItem("Grouped Clusters", tabName = "tab10"),
               menuSubItem("Interactive Visualization", tabName = "tab11")
      ),
      menuItem("Phylogenetic Analysis", icon = icon("network-wired"),
               menuSubItem("Phylogenetic plot", tabName = "tab12")
      )
    )
  )
## 3. body ------------------------------
body <- dashboardBody(
  tabItems(
    ##### UI: Input Data ########
    tabItem(tabName = "tab1",
            tags$style(button_color_css),
            titlePanel("Select Data and Pick Information as meta_info"),
            sidebarLayout(
              sidebarPanel(
                box(title = "", width = NULL, solidHeader = TRUE, status = "primary",
                fileInput("fileInput", "Choose CSV or TSV File",
                          accept = c(".csv", ".tsv")
                )
                ),
                hr(),
                box(title = "Selecting content for meta_info column:", width = NULL, solidHeader = TRUE, status = "info",
                uiOutput("columnsList"),
                actionButton("reset", "Reset Selections", icon = icon("refresh")))
              ),
              mainPanel(
                hr(),
                box(title = "Original Data Summary and View", width = NULL, solidHeader = TRUE, status = "primary",
                dataTableOutput("dataTable")),
                hr(),
                br(),
                box(title = "Information of meta_info View", width = NULL, solidHeader = TRUE, status = "info",
                tableOutput("selectedDataTable"))
              )
            )
            
            
    ),
    ##### UI: Process Data######## 
    tabItem(tabName = "tab2",
            tags$style(button_color_css),
            titlePanel("Filter and Prepare Data Format for Anchor Clustering:"),
            sidebarLayout(
              sidebarPanel(
                box(title = "Select columns to display:", width = NULL, solidHeader = TRUE, status = "primary",
                actionButton("reset2", "Reset Selections", icon = icon("refresh"))),
                hr(),
                box(title = "Filter by unique values where applicable:", width = NULL, solidHeader = TRUE, status = "primary",
                uiOutput("filterColumnsList"),
                uiOutput("columnSelectUI")),
                hr(),
                box(title = "Save the table for Anchor Clustering:", width = NULL, solidHeader = TRUE, status = "success",
                downloadButton("downloadData", "Download .tsv"))
              ),
              mainPanel(
                box(title = "Filtered Data View", width = NULL, solidHeader = TRUE, status = "primary",
                    dataTableOutput("filteredDataTable")),
                hr(),
                box(title = "Data Format in Anchor Clustering: sequence_id,  v_call,  j_call,  junction_length,  junction,  meta_info", width = NULL, solidHeader = TRUE, status = "success",
                    dataTableOutput("AnchorClusteringDataTable"))
              )
            )
    ),
    
    ##### UI: Integrate Current Data ########
    tabItem(tabName = "tab3",
            tags$style(button_color_css),
            titlePanel("Select Data to Integrate with the Input Data:"),
            sidebarLayout(
              sidebarPanel(
                box(title = "Select Data to integrate:", width = NULL, solidHeader = TRUE, status = "primary",
                actionButton("reset3", "Reset Selections", icon = icon("refresh")),
                uiOutput("filterColumnDefaultList"),
                uiOutput("columnUIdefault")),
                hr(),
                box(title = "Select columns to add in meta_info:", width = NULL, solidHeader = TRUE, status = "info",
                uiOutput("meta_infoList"),
                actionButton("resetBtn", "Reset Selections"),
                actionButton("processDF", "Process for meta_info")),
                hr(),
                box(title = "Merge selected data with your input file", width = NULL, solidHeader = TRUE, status = "danger",
                actionButton("mergeBtn", "Merge Datasets"),
                textOutput("message"),
                downloadButton("downloadMergeData", "Download Merged Data .tsv")),
                hr()
              ),
              mainPanel(
                box(title = "Filtered Data View", width = NULL, solidHeader = TRUE, status = "primary",
                dataTableOutput("dataTabledefault")),
                hr(),
                box(title = "Data Format in Anchor Clustering: sequence_id,  v_call,  j_call,  junction_length,  junction,  meta_info", width = NULL, solidHeader = TRUE, status = "info",
                dataTableOutput("DefaultACdata"))
              )
            )
     
    ),
    ##### UI: Run Anchor Clustering ########
    tabItem(tabName = "tab4",
            tags$style(button_color_css),
            titlePanel("Select Parameters and Run Anchor Clustering"),
           
            fluidRow(
              column(
                4,
                wellPanel(
                  box(
                    title = "Note the default parameters for applying Anchor Clustering on CDR3 amino acid sequences:", width = NULL, background = "green",
                  hr(),
                  h5(HTML("Population size: 1000, <br>
                           Random Material Rate: 50, <br>
                           Minimum Distance Ratio: 0.6, <br>
                           BIRCH Radius: 0.5, <br>
                           Fraction Data: 1, <br>
                           Size Threshold: 1000, <br>
                           Clustering threshold: 1, <br>
                           Distance Type: Hamming distance, <br>
                           Linkage Method: Single linkage, <br>
                           VJ Grouping: VJ <br> ")),
                  hr(),
                  h5(HTML("More specific usage methods, please check https://github.com/skylerchang/Anchor_Clustering_Nt"))),
                  
                  box(
                    title = "Select Anchor Clustering Parameters", width = NULL, solidHeader = TRUE, status = "primary",
                  numericInput("p", label = "Population Size", value = 1000),
                  numericInput("r", label = "Random Material Rate", value = 50, min = 25, max = 100, step = 25),
                  numericInput("m", label = "Minimum Distance Ratio", value = 0.6, min = 0.5, max = 0.9, step = 0.1),
                  numericInput("b", label = "BIRCH Radius", value = 0.5),
                  numericInput("f", label = "Fraction Data", value = 1, max = 1),
                  numericInput("z", label = "Size Threshold", value = 1000, max = 5000),
                  numericInput("t", label = "Clustering Threshold", value = 1),
                  radioButtons("s", label = "VJ Grouping", choices = list("V", "VJ", "None")),
                  radioButtons("d", label = "Distance Type", choices = list(
                    "Hamming distance" = "hd",
                    "Normalized Hamming distance" = "norm_hd"
                  )),
                  selectInput("l", "Linkage Method", choices = list(
                    "Single" = "single",
                    "Avergae" = "average",
                    "Complete" = "complete",
                    "Ward" = "ward",
                    "Centroid" = "centroid",
                    "Median" = "median",
                    "Weighted" = "weighted"
                  ))
                  )  
                )
              ),
              hr(),
              br(),
              mainPanel(
                box(
                  title = "Step1: Python Path Detection", width = NULL, solidHeader = TRUE, status= "primary",
                selectInput("selectedPython", "Choose Python Version:", choices = NULL),
                verbatimTextOutput("selectedVersion")),
                hr(),
                hr(),
                box(
                  title = "Step2: Choose File for Anchor Clustering", width = NULL, solidHeader = TRUE, status= "primary",
                fileInput("Anchorfile", "Upload your file for Anchor Clustering:",
                          accept = c("text/tsv", "text/tab-separated-values,text/plain", ".tsv")
                )),
                hr(),
                hr(),
                box(
                  title = "Step3: Run Anchor Clustering", width = NULL, solidHeader = TRUE, status= "primary",
                helpText("Run Default Anchor Clustering Parameters"),
                actionButton("rundefault", "Run Default Parameters!"),
                hr(),
                helpText("Or Run Selected Anchor Clustering Parameters"),
                actionButton("runselect", "Run Selected Parameters!")),
                hr(),
                hr(),
                box(
                  title = "Output of Anchor Clustering:", width = NULL, solidHeader = TRUE, status= "success",
                uiOutput("scriptOutput"))
              )
            )    
     
    ),
    ##### UI: Read Anchor Clustering Results ########
    tabItem(tabName = "tab5",
            tags$style(button_color_css),
            titlePanel("Split the Meta-Information for Anchor Clustering Results"),
            sidebarLayout(
              sidebarPanel(
                box(title = "Split Meta-info Column:", width = NULL, solidHeader = TRUE, status = "primary",
                fileInput("AnchorResultFile", "Upload your file:"),
                h5("Preview your uploaded data"),
                actionButton("showData", "Show First 10 Rows"),
                uiOutput("column_name_inputs"),
              
                h5("Give the column names after splitting meta-info:"),
                actionButton("split", "Split Columns"),
             
                h5("Preview your data after splitting meta-info by clicking Show Button Again:")),
                
                box(title = "Generate Network Cluster Features:", width = NULL, solidHeader = TRUE, status = "info",
                actionButton("analysis", "Generate Network Parameters"))
              ),
              mainPanel(
                box(title = "Splitting meta-info in Anchor Clustering results:", width = NULL, solidHeader = TRUE, status = "primary",
                DT::dataTableOutput("dataHead")),
                box(title = "Key Parameters of Network Cluster Features:", width = NULL, solidHeader = TRUE, status = "info",
                DT::dataTableOutput("processedTable"))
              )
            )
    ),
    ##### UI: Clonotype Distribution ########
    tabItem(tabName = "tab6",
            tags$style(button_color_css),
            titlePanel("Clonotype Distribution"),
            column(
              3,
              wellPanel(
                box(title = "Select Plot Features:", width = NULL, solidHeader = TRUE, status = "success",
                selectInput(
                  inputId = "bar_type",
                  label = "Select Type of Bar Plot",
                  choices = c("dodge", "stack"),
                  selected = "stack"
                ), 
                selectInput("selectedfactor1",
                            "Select factor to display (x axis) :",
                            choices = NULL
                ),
                conditionalPanel(
                  condition = "input.selectedfactor1 == 'junction_length'",
                  numericInput("x_min", "Minimum junction length ", value = 0),
                  numericInput("x_max", "Maximum  junction length", value = 30)
                ),
                selectInput("selectedfactor2",
                            "Select another factor to display (fill):",
                            choices = NULL
                ),
                selectInput("colorvar", "Color Palette",
                            color_vars,
                            selected = "Set1"
                ),
                selectInput("sortOrder", "Sort Order",
                            choices = c("Descending" = "desc", "Ascending" = "asc", "None" = "none"),
                            selected = "none"
                ),
                
                textInput(inputId = "x_label",
                          label = "Enter text for X-axis label"),
                
                selectInput("yAxisType", "Y Axis Type",
                            choices = c("Count" = "count", "Percentage" = "percentage"),
                            selected = "count"
                )),
                
                box(title = "Adjust your plot size:", width = NULL, solidHeader = TRUE, status = "warning",
                sliderInput(
                  inputId = "axis_size",
                  label = "Size of All Text ",
                  min = 3,
                  max = 40,
                  value = 20
                ),
                sliderInput(
                  inputId = "width",
                  label = "Width of the output plot ",
                  min = 400,
                  max = 1200,
                  value = 700
                ),
                sliderInput(
                  inputId = "height",
                  label = "Height of the output plot ",
                  min = 400,
                  max = 1200,
                  value = 600
                )
                )

              )
            ),
            fluidRow(
              mainPanel(
                box(status = "success", width = 36, height = "700px",
                tabsetPanel(
                  tabPanel("Bar Plot", plotOutput(outputId = "plot_lengthdis")),
                  tabPanel(
                    "Download Plot",
                    fluidRow(
                      column(
                        3,
                        numericInput("user_width", label = "Enter Width in cm", value = 20)
                      ),
                      column(
                        3,
                        numericInput("user_height", label = "Enter height in cm", value = 20)
                      ),
                      column(
                        3,
                        radioButtons(
                          "extension",
                          "Save As:",
                          choices = c("png", "pdf", "jpeg"),
                          inline = TRUE
                        )
                      ),
                      column(
                        3,
                        downloadButton("download", label = "Download plot", class = "butt")
                      )
                    ),
                    br(),
                    br(),
                    hr()
                  )
                )
              )
              )
            )
    ),
    
    ##### UI: Clonotype Overlap (Upset) ########
    tabItem(tabName = "tab7",
            tags$style(button_color_css),
            titlePanel("Compare Clonotype Overlap:"),
            column(
              3,
              wellPanel(
                box(title = "Select Plot Features:", width = NULL, solidHeader = TRUE, status = "success",
                selectInput("selectedColumnOverlap",
                            "Choose a Column:",
                            choices = NULL
                ),
                actionButton("applyButtonOverlap", "Apply Selection"), # Choices will be updated server-side
                hr(),
                uiOutput("sets"),
                selectInput(
                  "intersection_assignment_type",
                  "Intersection assignment type",
                  choices = c(
                    `Highest-order (UpSet)` = "upset",
                    `All associated intersections` = "all"
                  ),
                  selected = "upset"
                ),
                uiOutput("nsets"),
                sliderInput(
                  "nintersects",
                  label = "Number of intersections",
                  min = 2,
                  max = 40,
                  step = 1,
                  value = 20
                )),
                box(title = "Show Plot Details:", width = NULL, solidHeader = TRUE, status = "warning",
                checkboxInput("set_sort", "Sort sets by size?", value = TRUE),
                checkboxInput("show_set_size", "Show set sizes?", value = TRUE),
                checkboxInput("bar_numbers",
                              "Show bar numbers?",
                              value = TRUE
                ),
                checkboxInput(
                  "show_empty_intersections",
                  label = "Show empty intersections?",
                  value = TRUE
                )
                )
              )
            ),
            
            fluidRow(
              mainPanel(
                box(status = "success",  width = 36, height = "700px",
                tabsetPanel(
                  tabPanel("Upset Plot", plotlyOutput(outputId = "plotly_upset",width = "800px", height = "400px")),
                  tabPanel(
                    "Download Plot",
                    fluidRow(
                      column(
                        3,
                        numericInput("user_width_up", label = "Enter Width in px", value = 1000)
                      ),
                      column(
                        3,
                        numericInput("user_height_up", label = "Enter height in px", value = 800)
                      ),
                      column(
                        3,
                        radioButtons(
                          "extension_up",
                          "Save As:",
                          choices = c("png", "pdf", "jpeg", "tiff"),
                          inline = TRUE
                        )
                      ),
                      column(
                        3,
                        downloadButton("download_up", label = "Download plot")
                      )
                    ),
                    br(),
                    br(),
                    hr()
                  )
                )
              )
              )
        
            )
    ),
    
    ##### UI: Circos Plot #########
    tabItem(tabName = "tab8",
            tags$style(button_color_css),
            titlePanel("Circos Plot"),
            column(
              3,
              wellPanel(
                box(title = "Select Plot Features:", width = NULL, solidHeader = TRUE, status = "primary",
                    selectInput("selectedfactor_cp",
                                "Select group of count V and J:",
                                choices = NULL
                    ),
                    br(),
                    uiOutput("dataset"),
                    uiOutput("reads_ui"),
                    br()),
                
                box(title = "Select Plot Features:", width = NULL, solidHeader = TRUE, status = "success",
                    numericInput("top_cp", "Number of the maximum most V/J Pairs", value =10 , min = 1, max = 20),
                    selectInput("colorvar_Vcol", "V Gene Color Palette",
                                color_vars,
                                selected = "Set1"
                    ),
                    
                    selectInput("colorvar_Jcol", "J Gene Color Palette",
                                color_vars,
                                selected = "Set3"
                    )),
                
                box(title = "Adjust your plot size:", width = NULL, solidHeader = TRUE, status = "warning",
                    sliderInput(
                      inputId = "cex_cp",
                      label = "Size of the gene names",
                      min = 0.3,
                      max = 1.5,
                      value = 1
                    ),
                    
                    sliderInput(
                      inputId = "width_cp",
                      label = "Width of the output plot ",
                      min = 400,
                      max = 1200,
                      value = 700
                    ),
                    
                    sliderInput(
                      inputId = "height_cp",
                      label = "Height of the output plot ",
                      min = 400,
                      max = 1200,
                      value = 600
                    )
                )
                
              )
            ),
            fluidRow(
              mainPanel(
                box(status = "success", width = 36, height = "700px",
                    tabsetPanel(
                      tabPanel("Circos Plot", plotOutput(outputId = "circos_plot")),
                      tabPanel(
                        "Download Plot",
                        fluidRow(
                          column(
                            3,
                            numericInput("fheight_cp", "Height (cm)", min=2, max=15, step=1, value = 10)
                          ),
                          column(
                            3,
                            numericInput("fwidth_cp", "Width (cm)", min=2, max=15, step=1, value = 10)
                          ),
                          column(
                            3,
                            selectInput("fformat_cp", "File type", choices=c("png","tiff","jpeg"), 
                                        selected = "png", multiple = FALSE, selectize = TRUE)
                          ),
                          column(
                            3,
                            downloadButton('bn_download_cp', 'Download Plot')
                          ),
                          br(),
                          br(),
                          hr()
                        )
                      )
                    )
                )
              )
            )
    ),
    
    
    ##### UI: Scatter Plot  ########
    tabItem(tabName = "tab9",
            tags$style(button_color_css),
            titlePanel("Scatter Plot and Static Cluster Network"),
            column(3,
                   wellPanel(
                     box(title = "Select Plot Features:", width = NULL, solidHeader = TRUE, status = "success",
                     selectInput("xvar", "X-axis variable", axis_vars, selected = "mean_degree" ),
                     selectInput("yvar", "Y-axis variable", axis_vars, selected = "cluster_density" ),
                     selectInput("sizevar", "Size variable", axis_vars,  selected = "seqs" ),
                     selectInput("colorscatter", "Color Palette", color_vars, selected = "Set1"),
                     
                     tags$h4(
                       "Note: The input variables were shown when you hang over the nodes"
                     )),
                     br(),
                     br(),
                     
                   ),
                   
                   wellPanel(
                     box(title = "Select Network Plot Features:", width = NULL, solidHeader = TRUE, status = "info",
                     selectInput("group",
                                 "Choose group for the node labels:",
                                 choices = NULL
                     ),
                     selectInput("colors2", "Color Palette", color_vars, selected = "Paired"),
                     selectInput("layouts", "Network Layout", layouts, selected = "fr"),
                     sliderInput(
                       inputId = "scalesize",
                       label = "Size of node ",
                       min = 0.1,
                       max = 10,
                       value = 1.5
                     ),
                     numericInput("seed", "Layout seed", value = 123, min = 1, max = 20000),
                     sliderInput(
                       inputId = "legendfontsize",
                       label = "Size of Network legend ",
                       min = 0.1,
                       max = 10,
                       value = 2
                     )),
                     
                     box(title = "Adjust Network and Seqlogo Plot Size:", width = NULL, solidHeader = TRUE, status = "warning",
                     sliderInput("height1", "Height of the Network plot (px):", 
                                 min = 100, max = 1200, value = 500),
                     
                     # Slider for width
                     sliderInput("width1", "Width of the Network plot (px):", 
                                 min = 100, max = 1200, value = 400),
                     
                     sliderInput("height2", "Height of the Seqlogo plot (px):", 
                                 min = 100, max = 1200, value = 150),
                     
                     # Slider for width
                     sliderInput("width2", "Width of the Seqlogo plot (px):", 
                                 min = 100, max = 1200, value = 450),
                     
                     sliderInput(
                       inputId = "seqlogofontsize",
                       label = "Size of Seqlogo title ",
                       min = 1,
                       max = 20,
                       value = 12
                     ),
                     
                     br(),
                     br(),
                     actionButton("resetnet", "Reset Plots"),
                     actionButton("undo", "Undo")
                   ))
            ),
            
            fluidRow(
              mainPanel(
                tabPanel(
                  "Tab 1",
                ),
                column(width = 12,
                       h4("Network Parameter Scatter Plot"),
                       box(status = "success",  width = 36, height = "500px",
                       plotlyOutput("plot1", width = "750px", height = "500px")),
                       br(),
                       box(title = "Click the cluster in scatter plot to show network connections:",width = 36, background = "aqua"),
                       br(),
                       box(status = "info",  width = NULL,
                       uiOutput("additionalPlots"))
                )
              )
              )
    ),
    
    ##### UI: Static Grouped Visualization  ########
    tabItem(tabName = "tab10",
            tags$style(button_color_css),
            titlePanel("Show Grouped Clusters"),
            column(3,
                   wellPanel(
                     box(title = "Select Plot Features:", width = NULL, solidHeader = TRUE, status = "success",
                     numericInput("spread_lower", "Spread Lower Limit", value = 0, min = 0, max = 100),
                     numericInput("spread_upper", "Spread Upper Limit", value = 0, min = 0, max = 100),
                     numericInput("size_lower", "Size Lower Limit (log10)", value = 0, min = 0, max = 5),
                     numericInput("size_upper", "Size Upper Limit (log10)", value = 0, min = 0, max = 5),
                     br(),
                     br(),
                     selectInput("colors3", "Color Palette", color_vars, selected = "Set3"),
                     selectInput("layouts3", "Network Layout", layouts, selected = "fr"),
                     sliderInput(
                       inputId = "scalesize3",
                       label = "Size of node ",
                       min = 0.1,
                       max = 10,
                       value = 1.5
                     ),
                     
                     numericInput("seed3", "Layout seed", value = 123, min = 1, max = 20000),
                     
                     br(),
                     br(),
                     actionButton("goButton", "Generate Network")),
                     
                     box(title = "Adjust Plot Size:", width = NULL, solidHeader = TRUE, status = "warning",
                     sliderInput("height3", "Height of the Network plot (px):", 
                                 min = 100, max = 1200, value = 600),
                     
                     # Slider for width
                     sliderInput("width3", "Width of the Network plot (px):", 
                                 min = 100, max = 1200, value = 600),
                     br(),
                     br(),
                   )
                   )
            ),
            fluidRow(
              mainPanel(
                box(status = "success", width = 36, height = "800px",
                tabsetPanel(
                  tabPanel("Spread Distribution of all the clusters",
                    fluidRow(
                    column(
                      width = 12,
                      plotlyOutput("groupSpread", width = "750px", height = "500px")),
                    )
                    ),
                  tabPanel("Grouped Clusters", plotOutput(outputId ="groupVis")),
                  tabPanel(
                    "Download Plot",
                    fluidRow(
                      column(
                        3,
                        numericInput("fheight", "Height (cm)", min=2, max=15, step=1, value = 10)
                      ),
                      column(
                        3,
                        numericInput("fwidth", "Width (cm)", min=2, max=15, step=1, value = 10)
                      ),
                      column(
                        3,
                        selectInput("fformat", "File type", choices=c("png","tiff","jpeg","pdf"), 
                                    selected = "png", multiple = FALSE, selectize = TRUE)
                      ),
                      column(
                        3,
                        downloadButton('bn_download', 'Download Plot')
                    ),
                    br(),
                    br(),
                    hr()
                  )
                )
              )
            )
              )
            )
    ),
    
    #### UI: Interactive Visualization  ########
    tabItem(tabName = "tab11",
            tags$style(button_color_css),
            titlePanel("Interactive Visualization by Cluster"),
            column(3,
                   wellPanel(
                     box(title = "Select Plot Features:", width = NULL, solidHeader = TRUE, status = "success",
                     textInput("vis_cluster", "Cluster ID:"),
                     selectInput("vis_group", "Node Group",choices = NULL),
                     selectInput("vis_colors", "Color Palette", color_vars, selected = "Paired"),
                     selectInput("vis_layout", "Network Layout", layout_vars, selected = "layout_with_fr"),
                     tags$small(
                       paste0(
                         "Note: The Interactive Plot I were shown based on your specified layout, <br>
                          The Interactive Plot II were completely fluid."
                       )
                     ),
                     br(),
                     br(),
                     downloadButton('download_trend', 'Download cluster CSV data', class = "down")
                   )
                   )
            ),
            
            fluidRow(
              mainPanel(
                box(status = "success", width = 28, height = "700px",
                tabsetPanel(
                  tabPanel("Interactive Plot I", visNetworkOutput("networkVis2")),
                  tabPanel("Interactive Plot II", visNetworkOutput("networkVis1")),
                  tabPanel("Cluster table", DT::dataTableOutput('tbl'))
                )
                )
              )
            )
    ),
    ##### UI: Phylogenetic #######
    tabItem(tabName = "tab12",
            tags$style(button_color_css),
            titlePanel("Phylogenetic Plot"),
            column(3,
                   wellPanel(
                     box(title = "Select Phylogenetic Tree Features:", width = NULL, solidHeader = TRUE, status = "success",
                     textInput("phy_cluster", "Cluster ID:"),
                     selectInput("phy_group", "Color Group",choices = NULL),
                     selectInput("phy_color", "Color Palette", color_vars, selected = "Set1"),
                     selectInput("phy_layout", "Phylogenetic Tree Layout", ggtree_layouts, selected="circular"),
                     conditionalPanel(
                       condition = "input.phy_layout == 'fan'",
                       sliderInput("open_angle", "Open Angle", min = 0, max = 360, value = 120)
                     ),
                     conditionalPanel(
                       condition = "input.phy_layout == 'ellipse' || input.phy_layout == 'fan' || input.phy_layout == 'circular'",
                       selectInput("branch_length", "Branch Length", choices = c("none", "branch.length"), selected = "none")
                     )),
                     br(),
                     br(),
                     box(title = "Adjust Phylogenetic Tree Plot:", width = NULL, solidHeader = TRUE, status = "warning",
                     sliderInput(
                       inputId = "phy_size",
                       label = "Line Width of Branch",
                       min = 1,
                       max = 3,
                       value = 1.5
                     ),
                     
                     sliderInput(
                       inputId = "phy_node_size",
                       label = "Size of Node Label ",
                       min = 1,
                       max = 10,
                       value = 2
                     ),
                     
                     checkboxInput("include_msa", "Include MSA Plot", value = FALSE),
                     conditionalPanel(
                       condition = "input.include_msa == true",
                     sliderInput(
                       inputId = "msa_off",
                       label = "Distance between Tree and Multiple Sequence Alignment (MSA)",
                       min = 0.1,
                       max = 10,
                       value = 2
                     ),
                     
                     sliderInput(
                       inputId = "msa_width",
                       label = "Width of MSA plot",
                       min = 1,
                       max = 15,
                       value = 5
                     )
                     ),
                                      
                     
                     sliderInput(
                       inputId = "width_phy",
                       label = "Width of the Tree Plot ",
                       min = 400,
                       max = 1200,
                       value = 600
                     ),
                     
                     sliderInput(
                       inputId = "height_phy",
                       label = "Height of the Tree Plot ",
                       min = 400,
                       max = 1200,
                       value = 600
                     )  
                   )
                   )
            ),
            
            fluidRow(
              mainPanel(
                box(status = "success", width = 36, height = "700px",
                tabsetPanel(
                  tabPanel("Phylogenetic Plot", plotOutput(outputId = "phy_plot")),
                  tabPanel(
                    "Download Plot",
                    fluidRow(
                      column(
                        3,
                        numericInput("user_width_phy", label = "Enter Width in cm", value = 20)
                      ),
                      column(
                        3,
                        numericInput("user_height_phy", label = "Enter height in cm", value = 20)
                      ),
                      column(
                        3,
                        radioButtons(
                          "extension_phy",
                          "Save As:",
                          choices = c("png", "pdf", "jpeg"),
                          inline = TRUE
                        )
                      ),
                      column(
                        3,
                        downloadButton("download_phy", label = "Download plot", class = "butt")
                      )
                    ),
                    br(),
                    br(),
                    hr()
                  )
                )
                )
              )
            )
    )
  )
)







## UI ------------------------------
## 

ui <- tags$body(
  #tags$img(src = "www/plot2_logo.png", width = '60px'),
  #tags$span("KellerLab", style = "margin-left:25px;"),
  dashboardPage(header, siderbar, body,skin = "purple"))