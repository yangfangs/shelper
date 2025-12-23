# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# ============================================================
# 性能优化版本
# ============================================================

library(shiny)
library(DT)
library(shinyFiles)
library(shinyjs)

shinyUI(fluidPage(
  useShinyjs(),
  tags$head(
    # 预加载关键CSS
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    # DNS预解析和预连接（如有外部资源）
    tags$link(rel = "dns-prefetch", href = "//cdn.datatables.net"),
    tags$style(HTML("
      /* 加载提示容器 */
      .progress-container {
        display: none;
        position: fixed;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        z-index: 9999;
        background: white;
        padding: 30px;
        border-radius: 10px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        min-width: 300px;
        /* GPU加速 */
        will-change: opacity;
      }
      .progress-overlay {
        display: none;
        position: fixed;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background: rgba(0, 0, 0, 0.5);
        z-index: 9998;
        /* GPU加速 */
        will-change: opacity;
      }
      .loading-spinner {
        border: 4px solid #f3f3f3;
        border-top: 4px solid #3498db;
        border-radius: 50%;
        width: 40px;
        height: 40px;
        animation: spin 1s linear infinite;
        margin: 0 auto 15px;
        /* GPU加速 */
        will-change: transform;
      }
      @keyframes spin {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
      }
      /* 优化表格性能 */
      .dataTables_wrapper {
        contain: layout style;
      }
      /* 优化Tab切换 */
      .tab-pane {
        contain: content;
      }
    "))
  ),
  div(class = "progress-overlay", id = "progress-overlay"),
  div(class = "progress-container", id = "progress-container",
    div(class = "loading-spinner"),
    h4(id = "progress-text", style = "text-align: center; color: #333;", "Processing...")
  ),
 
  # Application title
  titlePanel("Shelper Tool"),
  
  # Main tabs or panels for Primer and Variant
  tabsetPanel(
    tabPanel("Primer-Design",
             sidebarLayout(
               sidebarPanel(
                 # 添加单选按钮
                 radioButtons("genomeVersion2", "1. Select Genome Version:",
                              choices = c("hg19 (GRCh37)" = "hg19", "hg38 (GRCh38)" = "hg38", "Custom Sequence" = "Custom"),
                              selected = "hg19"),
                 
                 conditionalPanel(
                    condition = "input.genomeVersion2 == 'Custom'",
                    textAreaInput("custom_seq_primer", "Paste Custom Sequence:", 
                                  placeholder = "Paste your reference sequence here...", 
                                  rows = 5)
                 ),
                 
                 # Input field for location
                 textInput("location_primer", "2. Input Genome Location:"),
                 bsTooltip("location_primer", "For Custom: Use 'Name-Pos-Ref-Alt' format. Pos is the 1-based index in your pasted sequence.", "right", options = list(container = "body")),
                 
                 # # File upload control 
                 # fileInput("file_primer", "File upload:"),
                 # Submit button
                 
                 checkboxInput('example1', tags$span(style="color: red", 'Load Example Data'), FALSE),
                 actionButton("submit_primer", "Submit")
               ),
               mainPanel(
                 # Nested tabs for Plot and Results
                 tabsetPanel(
                   id = "mainTabset2",  # 确保这个 ID 是唯一的
                   tabPanel("Help", includeHTML("primerHelp.html")),
                   # tabPanel("primer", tableOutput("primer"))  # 确保 tabPanel 的标题是 "primer"
                   tabPanel("primer",  uiOutput("tables"))  # 确保 tabPanel 的标题是 "primer"
                  
                 )
               )
             )
    ),
    tabPanel("Variant-Validation",
             sidebarLayout(
               sidebarPanel(
                 # 添加单选按钮
                 radioButtons("genomeVersion", "1. Select Genome Version:",
                              choices = c("hg19 (GRCh37)" = "hg19", "hg38 (GRCh38)" = "hg38", "Custom Sequence" = "Custom"),
                              selected = "hg19"),
                              
                 conditionalPanel(
                    condition = "input.genomeVersion == 'Custom'",
                    textAreaInput("custom_seq_variant", "Paste Custom Sequence:", 
                                  placeholder = "Paste your reference sequence here...", 
                                  rows = 5)
                 ),
                 
                 # 文本数据输入框
                 textAreaInput("textData", "2. Input vilidate location and file name (split by comma ',' )", "", rows = 10),
                 bsTooltip("textData", "For Custom: Use 'Name-Pos-Ref-Alt' format in your CSV/Text input.", "right", options = list(container = "body")),
                 
                 # 文件上传按钮
                 fileInput("file1", "3. Upload Chromatogram File (.ab1 or .scf)", 
                           multiple = TRUE,
                           accept = c(".ab1")),
                 checkboxInput('example', tags$span(style="color: red", 'Load Example Data'), FALSE),
                 actionButton("submit", "Submit"),
                 # div(class = "loader", id = "loading", style = "visibility: hidden;")
                 # div(class = "loader", id = "loading", "Now start validate，please wait a moment...")

               ),
               mainPanel(
                 
                 tabsetPanel(
                   id = "mainTabset",  # 添加此属性
                   # tabPanel("Help", tags$iframe(srcdoc = readLines("/Users/fangy/Desktop/sanger/SangerPriVar/Validation.html", warn = FALSE))),
                   tabPanel("Help", includeHTML("Validation2.html")),
                   
                   tabPanel("Results", 
                            uiOutput("dynamicTable"),  # 动态生成的表格
                            tags$br(style="clear:both"),
                            plotOutput("fig_text"),
                            tags$br(style="clear:both"),
                            uiOutput("dynamicUI")  # 动态生成的UI元素
                            )
                  
                 )
                 
              
               )

             )
    )
  )
))
