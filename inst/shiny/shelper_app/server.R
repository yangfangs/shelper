#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# ============================================================
# 性能优化版本
# ============================================================

library(shiny)
library(knitr)
library(sangerseqR)
library(future)
library(promises)
library(future.apply)
library(Shelper) # Ensure the package is loaded

# 根据CPU核心数设置并行workers数量，保留1-2个核心给系统
n_workers <- max(1, parallel::detectCores() - 2)
# Check if plan is already set to avoid warning?
# plan(multisession, workers = n_workers) 
# It is better to let the user or global setting handle this, 
# but for standalone app behavior, we can set it.
tryCatch(plan(multisession, workers = n_workers), error = function(e) warning(e))

shinyServer(function(input, output, session) {
  
  # ============================================================
  # 初始化 reactive values 用于存储计算结果，避免重复计算
  # ============================================================
  
  # 存储按钮点击的观察器，避免重复创建
  button_observers <- reactiveValues(created = FALSE)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~primer~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # primer_reactive <- eventReactive(input$submit_primer, {
  #   primerRes(input$location_primer,input$genomeVersion2)
  # })
  
  
  
  
  observe({
    req(input$example1)  # 确保在使用之前输入已经被提供
    if(input$example1) {
      updateTextInput(session, "location_primer", value = "13-32914859-G-GA")
      # print(input$location_primer)
    }
  })
  
  
  
  # Define the reactive value that captures the result of future computation
  primer_reactive <- eventReactive(input$submit_primer, {
    location_primer <- input$location_primer  # Assuming input UI element exists
    genomeVersion2 <- input$genomeVersion2    # Assuming input UI element exists
    
    future({
      # Simulating a function that fetches primer data
      # Load package in worker if needed
      library(Shelper)
      primerRes(location_primer, genomeVersion2)
    }) %...!% {
      # Catching exceptions from the future computation
      shinyjs::hide("progress-container")
      shinyjs::hide("progress-overlay")
      warning("Error in future: ", .)
      NULL
    }
  })
  
  # Observing the primer submission button
  observeEvent(input$submit_primer, {
    # 判断 input$location_primer 是否为空
    if (is.null(input$location_primer) || input$location_primer == "") {
      showModal(modalDialog(
        title = "Warning:",
        "Please input genome location",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return()
    }
    
    # 显示加载提示
    shinyjs::html("progress-text", "Designing primers, please wait...")
    shinyjs::show("progress-overlay")
    shinyjs::show("progress-container")
    
    updateTabsetPanel(session, "mainTabset2", selected = "primer")
    req(primer_reactive())  # Ensure primer_reactive() is not NULL or execution stops
    
    # Handle the results once the future resolves
    primer_reactive() %...>% {
      # 隐藏加载提示
      shinyjs::hide("progress-container")
      shinyjs::hide("progress-overlay")
      
      if (!is.data.frame(.)) {
        warning("result_primers is not a data frame")
        return(NULL)
      }
      
      split_primers <- split(., (seq(nrow(.)) - 1) %/% 2)
      
      output$tables <- renderUI({
        table_output_list <- lapply(seq_along(split_primers), function(i) {
          wellPanel(
            h4(paste("Primer pair", i)),
            DT::dataTableOutput(outputId = paste0("table_", i))
          )
        })
        do.call(tagList, table_output_list)
      })
      
      lapply(seq_along(split_primers), function(i) {
        output[[paste0("table_", i)]] <- DT::renderDataTable({
          split_primers[[i]]
        }, options = list(
          paging = FALSE, 
          searching = FALSE, 
          info = FALSE, 
          lengthChange = FALSE,
          # 优化: 减少DOM操作
          deferRender = TRUE
        ), server = FALSE)  # 小数据量使用客户端处理更快
      })
    }
  })
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~example for validate~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # 反应性值，存储textAreaInput的内容
    textContent <- reactiveVal("")
    
    observeEvent(input$example, {
      # 当example复选框被选中时
      if (input$example) {
        # 定义示例数据
        exampleText <- "17-41245245-CT-C,example1.ab1\n13-32893229-G-A,example2.ab1\n13-32914859-G-GA,example3.ab1\n13-20763530-CACACGTTCTTGCAGCC-C,example4.ab1\n13-20763210-CGTT-CGTTCGTT,example5.ab1"
        # 更新反应性值
        textContent(exampleText)
      }
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~variant upload~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 用于存储数据的反应值
    values <- reactiveValues(data = NULL, showTable = FALSE)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # 存储已解析的数据 - 使用 reactiveVal 避免重复计算
    resolved_data <- reactiveVal(NULL)
    
    # 触发计算的标志
    compute_trigger <- reactiveVal(0)
    
    observeEvent(input$submit, {
      
      # 显示加载提示
      shinyjs::html("progress-text", "Validating variants, please wait...")
      shinyjs::show("progress-overlay")
      shinyjs::show("progress-container")
      
      current_data <- NULL
      
      if(input$example) {
        updateTextAreaInput(session, "textData", value = textContent())
        values$showTable <- TRUE
        updateTabsetPanel(session, "mainTabset", selected = "Results")
        # 加载示例数据
        exampleData <- data.frame(GenomicCoordinate = c("17-41245245-CT-C", "13-32893229-G-A","13-32914859-G-GA",'13-20763530-CACACGTTCTTGCAGCC-C',"13-20763210-CGTT-CGTTCGTT"),
                                  FileName = c("example1.ab1", "example2.ab1","example3.ab1","example4.ab1","example5.ab1"),
                                  FilePath = c("example1.ab1", "example2.ab1","example3.ab1","example4.ab1","example5.ab1"))
        values$data <- exampleData
        current_data <- exampleData
        
      } else{
        # 判断 input$location_primer 是否为空
        if (is.null(input$textData) || input$textData == "") {
          # 隐藏加载提示
          shinyjs::hide("progress-container")
          shinyjs::hide("progress-overlay")
          showModal(modalDialog(
            title = "Warning:",
            "please input vilidate location and file name",
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
          return()
        }
        
        
        # selectedGenomeVersion <- input$genomeVersion
        #切换标签
        updateTabsetPanel(session, "mainTabset", selected = "Results")
        #控制table显示
        values$showTable <- TRUE
        # 确保输入了文本数据
        if (nchar(input$textData) == 0) {
          values$data <- data.frame(Error = "未输入文本数据")
        } else {
          # 解析文本数据
          lines <- strsplit(input$textData, "\n")[[1]]
          parsedData <- data.frame(matrix(ncol = 2, nrow = length(lines)))
          colnames(parsedData) <- c("GenomicCoordinate", "FileName")
          
          for (i in 1:length(lines)) {
            parsedData[i,] <- strsplit(lines[i], ",")[[1]]
          }
          
          # 将数据与文件名结合（如果已上传文件）
          if (!is.null(input$file1)) {
            # 创建一个以上传文件名为键，以文件路径为值的映射
            fileMap <- setNames(input$file1$datapath, input$file1$name)
            
            # 将文件路径添加到解析数据中
            parsedData$FilePath <- fileMap[parsedData$FileName]
            
            values$data <- parsedData
            current_data <- parsedData
            
          } else {
            values$data <- data.frame(Error = "未上传文件")
            # 隐藏加载提示
            shinyjs::hide("progress-container")
            shinyjs::hide("progress-overlay")
            return()
          }
        }
      }
      
      # 更新表格输出
      output$table1 <- renderTable({
        values$data
      })
      
      # 如果有有效数据，执行异步计算
      if (!is.null(current_data) && nrow(current_data) > 0 && !"Error" %in% names(current_data)) {
        input_genomeVersion <- input$genomeVersion
        
        # 异步执行计算
        future({
          library(Shelper)
          pickchromatogram2(current_data, input_genomeVersion)
        }) %...>% {
          result <- .
          # 存储解析后的数据
          resolved_data(result)
          # 隐藏加载提示
          shinyjs::hide("progress-container")
          shinyjs::hide("progress-overlay")
          # 触发UI更新
          compute_trigger(compute_trigger() + 1)
        } %...!% {
          # 处理错误
          shinyjs::hide("progress-container")
          shinyjs::hide("progress-overlay")
          warning("Error in validation: ", .)
          resolved_data(NULL)
        }
      }
    })
    
  
    # 使用服务器端处理的 DataTable - 提升大数据量时的性能
    output$table <- DT::renderDataTable({
      # 依赖 compute_trigger 和 resolved_data
      compute_trigger()
      req(resolved_data())
      
      valid_data <- resolved_data()
      req(!is.null(valid_data) && !is.null(valid_data$res_table))
      
      data <- valid_data$res_table
      
      # 转换Figure列和Check列
      data$Figure <- sapply(data$Figure, function(i) {
        as.character(actionButton(inputId = paste0("button", i), label = "Show"))
      })
      data$Check <- sapply(data$Check, function(check) {
        if (check) {
          as.character(tags$span("\u2714", style = "color:green;"))
        } else {
          as.character(tags$span("\u2718", style = "color:red;"))
        }
      })
      
      datatable(data, escape = FALSE, options = list(
        preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
        drawCallback = JS('function() { Shiny.bindAll(this.api().table().node()); }'),
        columnDefs = list(list(targets = 6, visible = TRUE)),
        # 服务器端处理选项
        processing = TRUE
      ))
    }, server = TRUE)  # 启用服务器端处理
    
    
    showTags <- reactiveVal(FALSE)
    
    # 存储当前选中的按钮索引
    selected_button <- reactiveVal(NULL)
    
    # 监听所有按钮点击 - 使用单一的观察器而非动态创建多个
    observe({
      compute_trigger()
      valid_data <- resolved_data()
      req(!is.null(valid_data) && !is.null(valid_data$res_table))
      
      n_rows <- nrow(valid_data$res_table)
      if (n_rows > 0) {
        # 为每个按钮创建观察器（只创建一次）
        lapply(seq_len(n_rows), function(i) {
          observeEvent(input[[paste0("button", i)]], {
            selected_button(i)
          }, ignoreInit = TRUE, once = FALSE)
        })
      }
    })
    
    # 当按钮被点击时更新显示
    observeEvent(selected_button(), {
      req(selected_button())
      valid_data <- resolved_data()
      req(!is.null(valid_data))
      
      i <- selected_button()
      showTags(TRUE)
      
      output$fig_text <- renderPlot({
        sangerseqR::chromatogram(valid_data$chromatogram_pick[[i]], width = 25, height = 3, trim5 = 0, trim3 = 0, showcalls = "both")
      })
      output$refseq <- renderText(valid_data$refseq[[i]])
      output$altseq <- renderText(valid_data$altseq[[i]])
      output$alignment <- renderText(
        gsub(pattern="^.*#={39}(.+?)#-{39}.*$",
             replacement="\\1",
             x=valid_data$alignment[[i]])
      )
      output$header <- renderText(
        gsub(pattern="(^.+)#\\n#\\n#={39}.+$",
             replacement="\\1",
             x=valid_data$alignment[[i]])
      )
    }, ignoreInit = TRUE)
    
    
    
    # 根据 showTags 的值动态生成UI
    output$dynamicUI <- renderUI({
      if(showTags()) {  # 如果 showTags 为 TRUE
        tagList(
          tags$h4("Alignment"), 
          verbatimTextOutput('alignment'),
          tags$h4("Reference Sequence"), 
          verbatimTextOutput('refseq'),
          tags$h4("Alternate Allele"),
          verbatimTextOutput('altseq'),
          tags$h4("Alignment Header"), 
          verbatimTextOutput('header')
        )
      }
    })
    # 根据 showTable 的值动态生成 dataTableOutput
    output$dynamicTable <- renderUI({
      if(values$showTable) {  
        tagList(
          tags$br(style="clear:both"),
          div(
            id = "download-buttons",
            downloadButton("downloadTable", "Download Table"),
            downloadButton("downloadPlots", "Download Plots"),
            style = "text-align: right;"  # 添加样式
          ),
          tags$h4("Viliation Result"),
          tags$br(),
          dataTableOutput("table")
        )
      } else {
        tags$p("Check results will be shown here when a valid sequencing file has been uploaded.")
      }
    })
      
    
    output$downloadTable <- downloadHandler(
      filename = function() {
        paste("table-data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        # 获取已解析的表格数据
        result_data <- resolved_data()
        if (!is.null(result_data) && !is.null(result_data$res_table)) {
          write.csv(result_data$res_table, file, row.names = FALSE)
        }
      }
    )
    
    save_plot <- function(data, folder, filename) {
      # 确保文件夹存在
      if (!dir.exists(folder)) {
        dir.create(folder)
      }
      
      # 完整的文件路径
      full_filename <- file.path(folder, filename)
      
      png(full_filename, width = 25*72, height = 3*72) # 72 是每英寸的像素数
      sangerseqR::chromatogram(data, width = 25, height = 3, trim5 = 0, trim3 = 0, showcalls = "both")
      dev.off()
    }
    
    
    output$downloadPlots <- downloadHandler(
      filename = function() {
        paste("Chromatogram-", Sys.Date(), ".zip", sep="")
      },
      content = function(zip_filename) {
        # 获取已解析的数据
        result_data <- resolved_data()
        
        if (!is.null(result_data) && !is.null(result_data$chromatogram_pick)) {
          # 创建一个临时目录来存储图表文件
          temp_dir <- tempdir()
          chromatogram_folder <- file.path(temp_dir, "Chromatogram")
          
          # 确保文件夹存在
          if (!dir.exists(chromatogram_folder)) {
            dir.create(chromatogram_folder, recursive = TRUE)
          }
          
          # 清空文件夹中的旧文件
          existing_files <- list.files(chromatogram_folder, full.names = TRUE)
          if (length(existing_files) > 0) {
            file.remove(existing_files)
          }
          
          # 为每个图表生成一个文件
          for (i in seq_along(result_data$chromatogram_pick)) {
            plot_filename <- paste0("plot_", i, ".png")
            save_plot(result_data$chromatogram_pick[[i]], chromatogram_folder, plot_filename)
          }
          
          # 获取当前工作目录
          old_wd <- getwd()
          
          # 将文件打包成 ZIP 文件
          # 切换到临时目录，这样 ZIP 文件中的路径不包含完整的临时路径
          setwd(temp_dir)
          utils::zip(zipfile = zip_filename, files = "Chromatogram", flags = "-r9Xq")
          setwd(old_wd)
        }
      }
    )
    
    
}
)
