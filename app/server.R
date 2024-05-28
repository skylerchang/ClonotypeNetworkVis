## Server -----------------------------
server <- function(input, output, session) {
  #### ----------------------------####
  #### Data Input and Processing #####
  #### ----------------------------####
  ######################## Input Data ######################

  data_input <- reactive({
    req(input$fileInput)
    ext <- tools::file_ext(input$fileInput$name)
    switch(ext,
      csv = vroom::vroom(input$fileInput$datapath, delim = ",", show_col_types = FALSE),
      tsv = vroom::vroom(input$fileInput$datapath, delim = "\t", show_col_types = FALSE),
      validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })

  output$dataTable <- DT::renderDataTable({
    req(data_input())
    datatable(data_input(), options = list(pageLength = 6,scrollX = TRUE))
  })

  output$columnsList <- renderUI({
    req(data_input())
    column_names <- names(data_input())
    lapply(seq_along(column_names), function(i) {
      checkboxInput(inputId = paste("firstPage_col", i, sep = "_"), label = column_names[i], value = TRUE)
    })
  })

  selectedData <- reactive({
    req(data_input())
    column_names <- names(data_input())

    # Simplify the access to inputs
    selected_columns <- vapply(column_names, function(col_name) {
      input_id <- paste("firstPage_col", which(column_names == col_name), sep = "_")
      if (!is.null(input[[input_id]])) {
        input[[input_id]]
      } else {
        FALSE # Default to FALSE if input does not exist
      }
    }, logical(1))

    if (any(selected_columns, na.rm = TRUE)) {
      selected_df <- data_input()[, selected_columns, drop = FALSE]
      selected_df$meta_info <- do.call(paste, c(selected_df, sep = "<->"))
      return(selected_df)
    } else {
      return(data.frame())
    }
  })

  output$selectedDataTable <- renderTable({
    req(selectedData())
    head(selectedData(), 6)
  })

  observeEvent(input$reset, {
    updateCheckboxGroupInput(session, "selectedColumns", selected = NULL)
    column_names <- names(data_input())
    lapply(seq_along(column_names), function(i) {
      updateCheckboxInput(session, paste("firstPage_col", i, sep = "_"), value = FALSE)
    })
  })

  ######################## Process Data ######################

  output$columnSelectUI <- renderUI({
    df <- data_input() # Use the full data to list all columns
    req(df)
    column_names <- names(df)
    checkboxGroupInput("selectedColumns", "Columns to display:", choices = column_names, selected = column_names)
    actionButton("applySelection", "Apply Selections", icon = icon("check"))
  })

  output$filterColumnsList <- renderUI({
    df <- data_input()
    req(df)
    column_names <- names(df)

    ui_elements <- lapply(column_names, function(col_name) {
      unique_values <- unique(df[[col_name]])
      if (length(unique_values) <= 6) {
        list(
          h5(col_name),
          checkboxGroupInput(
            inputId = paste("filter", col_name, sep = "_"),
            label = NULL,
            choices = unique_values,
            selected = unique_values
          )
        )
      }
    })
    do.call(tagList, ui_elements)
  })

  observeEvent(input$applySelection, {
    output$filteredDataTable <- DT::renderDataTable({
      df <- data_input()
      req(df)

      # Append the 'meta_info' column from the selected data
      selected_df <- selectedData()
      if (!is.null(selected_df) && "meta_info" %in% names(selected_df)) {
        df$meta_info <- selected_df$meta_info
      }

      # Filter the dataframe by selected columns from the checkboxGroupInput
      selected_columns <- input$selectedColumns
      if (!is.null(selected_columns)) {
        df <- df[, selected_columns, drop = FALSE]
      }

      # Apply filters based on selected unique values for each column
      column_names <- names(df)
      for (col_name in column_names) {
        filter_input_id <- paste("filter", col_name, sep = "_")
        selected_values <- input[[filter_input_id]]
        if (!is.null(selected_values)) {
          df <- df[df[[col_name]] %in% selected_values, , drop = FALSE]
        }
      }
      datatable(df, options = list(pageLength = 6, scrollX = TRUE))
      #head(df, 10)
    })
  })

  observeEvent(input$reset2, {
    df <- data_input()
    req(df)
    column_names <- names(df)

    # Update each filter control to have no selections
    for (col_name in column_names) {
      # Check if the column has associated filter controls
      if (length(unique(df[[col_name]])) <= 6) {
        updateCheckboxGroupInput(session, paste("filter", col_name, sep = "_"), selected = character(0))
      }
    }

    # Optionally, you can force a refresh of the table to display the full data set with 'meta_info'
    output$filteredDataTable <- DT::renderDataTable({
      # Re-fetch selected data to ensure 'meta_info' column is up-to-date
      selected_df <- selectedData()
      if (!is.null(selected_df) && "meta_info" %in% names(selected_df)) {
        df$meta_info <- selected_df$meta_info
      }
      
    })
  })

  filteredData <- reactive({
    df <- data_input()
    req(df)
    df$meta_info <- selectedData()$meta_info

    # # Apply filters based on selected unique values
    column_names <- names(df)
    for (col_name in column_names) {
      filter_input_id <- paste("filter", col_name, sep = "_")
      if (!is.null(input[[filter_input_id]]) && length(input[[filter_input_id]]) > 0) {
        df <- df[df[[col_name]] %in% input[[filter_input_id]], , drop = FALSE]
      }
    }
    df
  })


  output$filteredDataTable <- DT::renderDataTable({
    df <- filteredData()
    datatable(df, options = list(pageLength = 6, scrollX = TRUE
    ))
  })
  

  acData <- reactiveVal()
  observeEvent(input$applySelection, {
    df <- filteredData()

    ac <- df %>%
      select(
        v_call = matches("^(v_gene|v|V|V_gene|v_call|V_call)$"),
        j_call = matches("^(j_gene|j|J|J_gene|j_call|J_call)$"),
        junction_length = matches("^(cdr3_length|CDR3_length|CDR3 length|cdr3 length|CDR3 Length|junction_length|junction length)$"),
        junction = matches("^(cdr3|CDR3|junction|Junction|cdr3_aa)$"),
        meta_info
      ) %>%
      distinct() %>%
      filter(!is.na(junction)) %>%
      mutate(sequence_id = row_number()) %>%
      select(sequence_id, v_call, j_call, junction_length, junction, meta_info)

    acData(ac)
  })


  output$AnchorClusteringDataTable <- DT::renderDataTable({
    ac <- acData()
    datatable(ac, options = list(pageLength = 6))
  })


  output$downloadData <- downloadHandler(
    filename = function() {
      test <- acData()
      req(test)
      chain <- unique(gsub("J.*[0-9]", "", test$j_call))
      paste0("Process_for_AC_", chain, "_AA.tsv")
    },
    content = function(file) {
      vroom::vroom_write(acData(), file, delim = "\t", quote = "none")
    }
  )

  ######################## Integrate Current Data ######################

  datadef <- reactive({
    # Adjust the path to your actual default data location
    datadf <- read.table(here("data", "all.tsv"), stringsAsFactors = FALSE, sep = "\t", header = TRUE, fill = TRUE)
  })

  output$columnUIdefault <- renderUI({
    df <- datadef() # Use the full data to list all columns
    req(df)
    column_names <- names(df)
    checkboxGroupInput("selectedColumnshere", "Columns to display:", choices = column_names, selected = column_names)
    actionButton("applySelectionhere", "Apply Selections", icon = icon("check"))
  })

  output$filterColumnDefaultList <- renderUI({
    df <- datadef()
    req(df)
    column_names <- names(df)
    ui_elements <- lapply(column_names, function(col_name) {
      unique_values <- na.omit(unique(df[[col_name]]))
      unique_values <- unique_values[unique_values != ""]

      if (length(unique_values) <= 7) {
        list(
          h5(col_name),
          checkboxGroupInput(
            inputId = paste("filter", col_name, sep = "_"),
            label = NULL,
            choices = unique_values,
            selected = unique_values
          )
        )
      }
    })
    do.call(tagList, ui_elements)
  })

  observeEvent(input$reset3, {
    df <- datadef()
    req(df)
    column_names <- names(df)
    # Update each filter control to have no selections
    for (col_name in column_names) {
      # Check if the column has associated filter controls
      if (length(unique(df[[col_name]])) <= 7) {
        updateCheckboxGroupInput(session, paste("filter", col_name, sep = "_"), selected = character(0))
      }
    }
  })

  df_reactive <- reactiveVal()

  # Apply selection button functionality
  observeEvent(input$applySelectionhere, {
    df <- datadef()
    req(df)

    selected_columns <- input$selectedColumnshere
    if (!is.null(selected_columns)) {
      df <- df[, selected_columns, drop = FALSE]
    }

    column_names <- names(df)
    for (col_name in column_names) {
      filter_input_id <- paste("filter", col_name, sep = "_")
      selected_values <- input[[filter_input_id]]
      if (!is.null(selected_values)) {
        df <- df[df[[col_name]] %in% selected_values, , drop = FALSE]
      }
    }
    df_reactive(df)
  })


  output$dataTabledefault <- renderDataTable(
    req(df_reactive()),
    options = list(pageLength = 6,scrollX = TRUE)
  )


  output$meta_infoList <- renderUI({
    df <- df_reactive()
    if (is.null(df)) {
      return(NULL)
    }
    column_names2 <- names(df)
    column_names2 <- column_names2[-(1:4)]
    lapply(seq_along(column_names2), function(i) {
      checkboxInput(inputId = paste("secondPage_col", i, sep = "_"), label = column_names2[i], value = TRUE)
    })
  })


  observeEvent(input$resetBtn, {
    # updateCheckboxInput(session, "selectedColumns", selected = NULL)
    df <- df_reactive()
    if (is.null(df)) {
      return(NULL)
    }
    column_names2 <- names(df)
    column_names2 <- column_names2[-(1:4)]
    lapply(seq_along(column_names2), function(i) {
      updateCheckboxInput(session, inputId = paste("secondPage_col", i, sep = "_"), value = FALSE)
    })
  })

  selectedmeta_infoData <- reactiveVal()

  observeEvent(input$processDF, {
    df <- isolate(df_reactive())
    req(df)
    selected_columns2 <- sapply(1:length(names(df)[-(1:4)]), function(i) {
      isolate(input[[paste("secondPage_col", i, sep = "_")]])
    })
    names(selected_columns2) <- names(df)[-(1:4)]
    selected_columns2 <- names(selected_columns2[selected_columns2])


    if (length(selected_columns2) > 0) {
      selected_df2 <- df[, selected_columns2, drop = FALSE]
      selected_df2$meta_info <- do.call(paste, c(selected_df2, sep = "<->"))
      selectedmeta_infoData(selected_df2)
    } else {
      selectedmeta_infoData(data.frame())
    }
  })


  acdefault <- eventReactive(input$processDF, {
    df <- df_reactive()
    req(df)
    df_with_meta_info <- selectedmeta_infoData()
    req(df_with_meta_info)
    if ("meta_info" %in% names(df_with_meta_info)) {
      df$meta_info <- df_with_meta_info$meta_info
    } else {
      df$meta_info <- NA
    }
    df %>%
      select(v_call, j_call, junction, junction_length, meta_info, locus) %>%
      distinct() %>%
      mutate(sequence_id = row_number()) %>%
      select(sequence_id, v_call, j_call, junction_length, junction, meta_info)
  })



  output$DefaultACdata <- DT::renderDataTable({
    ac2 <- acdefault()
    datatable(ac2, options = list(pageLength = 6, scrollX = TRUE))

  })


  mergedData <- reactiveVal()

  observeEvent(input$mergeBtn, {
    df1 <- acData()
    df2 <- acdefault()
    req(df2)
    if (!is.null(df1) && nrow(df1) > 0) {
      merged <- rbind(df1, df2)
      merged <- merged %>%
        select(-sequence_id) %>%
        dplyr::mutate(sequence_id = row_number()) %>%
        select(sequence_id, v_call, j_call, junction_length, junction, meta_info)
      mergedData(merged)
      output$message <- renderText("Data has been merged successfully.")
    } else {
      merged2 <- df2
      mergedData(merged2)
      output$message <- renderText("No input data to provide, database data has been selected.")
    }
  })


  output$downloadMergeData <- downloadHandler(
    filename = function() {
      paste("Merged_dataset_for_AnchorClustering", Sys.Date(), ".tsv")
    },
    content = function(file) {
      vroom::vroom_write(mergedData(), file, delim = "\t", quote = "none")
    }
  )


  #################### Run Anchor Clustering ######################

  # List Python executables in the determined directory
  observe({
    python_path <- system("which python3", intern = TRUE)
    use_python(python_path)
    python_dir <- dirname(python_path)
    python_files <- list.files(path = python_dir, pattern = "^python", full.names = TRUE)
    updateSelectInput(session, "selectedPython", choices = python_files)
  })

  # Display the selected Python version
  output$selectedVersion <- renderText({
    paste("You have selected:", input$selectedPython)
  })


  observeEvent(input$rundefault, {
    req(input$selectedPython, input$Anchorfile)
    python_path <- input$selectedPython

    # Build the command
    command <- sprintf(
      '%s "%s" -F "%s" -p 1000 -r 50 -m 0.6 -b 0.5 -f 1 -z 1000 -t 1 -s VJ -d hd -l single',
      python_path, "Anchor_Clustering_PreVJ.py", input$Anchorfile
    )

    # Execute the command and capture output
    output$scriptOutput <- renderUI({
      tryCatch(
        {
          result <- system(command, intern = TRUE)
          paste("Command executed successfully:")
          HTML(paste(gsub("\n", "<br>", result), collapse = "<br>"))
        },
        error = function(e) {
          HTML(paste("Error during execution:", e$message))
        }
      )
    })
  })

  observeEvent(input$runselect, {
    # Ensure all inputs are available
    req(input$selectedPython, input$Anchorfile, input$p, input$r, input$m, input$b, input$f, input$z, input$t, input$s, input$d, input$l)

    # Construct the command string
    python_path <- input$selectedPython

    p <- input$p
    r <- input$r
    m <- input$m
    b <- input$b
    f <- input$f
    z <- input$z
    t <- input$t
    s <- input$s
    d <- input$d
    l <- input$l

    # Build the command
    command <- sprintf(
      '%s "%s" -F "%s" -p %d -r %d -m %.1f -b %.1f -f %.1f -z %d -t %d -s %s -d %s -l %s',
      python_path, "Anchor_Clustering_PreVJ.py", input$Anchorfile, input$p, input$r,
      input$m, input$b, input$f, input$z, input$t, input$s, input$d, input$l
    )

    # Execute the command and capture output
    output$scriptOutput <- renderUI({
      tryCatch(
        {
          result <- system(command, intern = TRUE)
          paste("Command executed successfully:")
          HTML(paste(gsub("\n", "<br>", result), collapse = "<br>"))
        },
        error = function(e) {
          HTML(paste("Error during execution:", e$message))
        }
      )
    })
  })



  ######################## Read Anchor Clustering Results ######################

  acresultdata <- reactiveVal()

  observeEvent(input$AnchorResultFile, {
    req(input$AnchorResultFile)
    df <- read.table(input$AnchorResultFile$datapath, sep = "\t", header = TRUE)
    acresultdata(df)
  })

  
  
  
  

  output$column_name_inputs <- renderUI({
    df <- acresultdata()
    req(df)
    if ("meta_info" %in% names(df)) {
      # Find the maximum number of dashes
      max_dashes <- max(sapply(df$meta_info, function(x) length(strsplit(as.character(x), "<->")[[1]]) - 1))
      # Create inputs for each potential new column
      lapply(1:(max_dashes + 1), function(i) { # +1 because n-dashes splits into n+1 pieces
        textInput(paste0("col_name", i), sprintf("Column %d name:", i))
      })
    }
  })


  # Observe button click to split columns
  observeEvent(input$split, {
    req(acresultdata())
    df <- acresultdata()
    # Prepare to split the 'meta_info' column and standardize the lengths
    max_dashes <- max(sapply(df$meta_info, function(x) length(strsplit(as.character(x), "<->")[[1]]) - 1))
    split_data <- t(sapply(df$meta_info, function(x) {
      pieces <- unlist(strsplit(as.character(x), "<->"))
      length(pieces) <- max_dashes + 1 # Ensure all are of the same length
      return(pieces)
    }))

    
    # Extract input names for new columns
    col_names <- sapply(1:ncol(split_data), function(i) input[[paste0("col_name", i)]])
    colnames(split_data) <- col_names
    df <- cbind(df, split_data)
    df$meta_info <- NULL
    rownames(df) <- NULL
    acresultdata(df) 
  })


  # Display the first 10 rows when the button is clicked
  observeEvent(input$showData, {
    req(acresultdata())
    output$dataHead <- DT::renderDataTable({
      req(acresultdata())
      datatable(acresultdata(),options = list(pageLength = 6,scrollX = TRUE))
    })
  })


  processed_reactive <- reactiveVal()

  # Apply selection button functionality
  observeEvent(input$analysis, {
    data <- acresultdata()
    req(data)
    data <- data %>% mutate(locus = gsub("V.*", "", v_call))
    prodata <- generate_network_parameters(data)
    processed_reactive(prodata)
  })


  output$processedTable <- DT::renderDataTable({
    req(processed_reactive()) 
    datatable(processed_reactive(), options = list(pageLength = 6,scrollX = TRUE ))
  })


  #### ----------------------------####
  #### Statistical Analysis     #####
  #### ----------------------------####
  ######################## Clonotype Distribution  ######################

  observe({
    
    if (is.null(acresultdata()) || nrow(acresultdata()) == 0) {
      # Initialize data_dis as an empty dataframe or handle the case accordingly
      data_dis <- data.frame()  # Create an empty dataframe with no columns
    } else {
      # Proceed with your existing code
      data_dis <- acresultdata()
      
      data_dis <- data_dis %>% 
        dplyr::mutate(v_family = gsub("-.*", "", v_call)) %>% 
        dplyr::mutate(j_family = gsub("-.*", "", j_call))
    }
    

    all_columns <- colnames(data_dis)

    # Specify columns to exclude
    exclude_columns <- c("sequence_id", "cluster", "junction", "reads", "read" )

    # Remove the specified columns from the list
    filtered_columns <- setdiff(all_columns, exclude_columns)

    # Update the selectInput with the filtered column names
    updateSelectInput(session, "selectedfactor1", choices = filtered_columns)
    updateSelectInput(session, "selectedfactor2", choices = filtered_columns)
  })

  plot_output <- reactive({
    data <- acresultdata() %>% dplyr::mutate(v_family = gsub("-.*", "", v_call)) %>% 
      dplyr::mutate(j_family = gsub("-.*", "", j_call))
    factor1 <- input$selectedfactor1 # Store reactive input in a local variable
    factor2 <- input$selectedfactor2

    data_new <- data %>%
      mutate(locus = gsub("V.*", "", v_call)) %>%
      group_by(locus, !!sym(factor1), !!sym(factor2)) %>% # Use !!sym() for dynamic grouping
      summarise(n = n_distinct(junction), .groups = "drop") # Ensure groups are dropped after summarise

    if (input$yAxisType == "percentage") {
      data_new <- data_new %>%
        group_by(locus, !!sym(factor1)) %>%
        mutate(percentage = n / sum(n) * 100) # Convert counts to percentages
    }


    if (input$sortOrder == "desc") {
      data_new <- data_new %>%
        group_by(locus, !!sym(factor1)) %>%
        mutate(total_n = sum(n)) %>%
        arrange(desc(factor1), -total_n)
      data_new[[factor1]] <- factor(data_new[[factor1]], unique(data_new[[factor1]]))
    } else if (input$sortOrder == "asc") {
      data_new <- data_new %>%
        group_by(locus, !!sym(factor1)) %>%
        mutate(total_n = sum(n)) %>%
        arrange(desc(factor1), total_n)
      data_new[[factor1]] <- factor(data_new[[factor1]], unique(data_new[[factor1]]))
    } else {
      data_new
    }


    palette <- brewer.pal(min(8, length(unique(data_new[[factor2]]))), input$colorvar)
    colors <- colorRampPalette(palette)(length(unique(data_new[[factor2]])))
    color_list <- setNames(colors, unique(data_new[[factor2]]))

    
    
    y_variable <- if(input$yAxisType == "percentage") "percentage" else "n"
    
    
    ggplot(data_new, aes(x = get(factor1), y = get(y_variable), fill = get(factor2))) +
      geom_bar(position = input$bar_type, stat = "identity") +
      scale_fill_manual(values = color_list) +
      facet_wrap(~locus) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank()
      ) +
      theme_bw() +
      labs(
        x = input$x_label, y = if (input$yAxisType == "percentage") "Percentage" else "Count",
        fill = factor2
      ) +
      theme(
        axis.title = element_text(angle = 0),
        text = element_text(size = input$axis_size, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold", size = input$axis_size)
      )


    # Return the plot object
  })

  output$plot_lengthdis <-
    renderPlot(
      plot_output(),
      width = function() {
        input$width
      },
      height = function() {
        input$height
      }
    )

  output$download <- downloadHandler(
    filename = function() {
      paste0(paste0("Clonotype_distribution.", input$extension))
    },
    content = function(file) {
      ggsave(
        file,
        plot = plot_output(),
        device = input$extension,
        width = input$user_width,
        height = input$user_height,
        units = "cm"
      )
    }
  )




  ######################## Shared Clonotypes with UpSetR  ######################

  # Render a control for users to pick from their supplied sets

  output$sets <- renderUI({
    valid_sets <- getValidSets()
    req(!is.null(valid_sets))

    selectInput(
      "sets",
      "Sets",
      choices = names(valid_sets),
      selectize = TRUE,
      multiple = TRUE,
      selected = names(valid_sets)
    )
  })

  # Render a control to decide how many sets to consider in the plot

  output$nsets <- renderUI({
    selected_sets <- getSelectedSets()
    req(!is.null(selected_sets))

    max_sets <-
      ifelse(length(selected_sets) > 20,
        20,
        length(selected_sets)
      )
    sliderInput(
      "nsets",
      label = "Number of sets to include in plot",
      min = 2,
      max = max_sets,
      step = 1,
      value = min(10, max_sets)
    )
  })

  # Accessor for user-selected sets

  getSelectedSetNames <- reactive({
    req(input$sets)
    input$sets
  })

  # Accessor for the nsets parameter

  getNsets <- reactive({
    req(!is.null(input$nsets))
    input$nsets
  })

  # Accessor for the nintersections parameter

  getNintersections <- reactive({
    validate(need(!is.null(input$nintersects), "Waiting for nintersects"))
    input$nintersects
  })

  getShowEmptyIntersections <- reactive({
    validate(need(
      !is.null(input$show_empty_intersections),
      "Waiting for empty intersections option"
    ))
    input$show_empty_intersections
  })

  # Accessor for the intersection assignment type

  getIntersectionAssignmentType <- reactive({
    validate(need(
      !is.null(input$intersection_assignment_type),
      "Waiting for group_by"
    ))
    input$intersection_assignment_type
  })

  # Set sorting

  getSetSort <- reactive({
    validate(need(!is.null(input$set_sort), "Waiting for set_sort"))
    input$set_sort
  })

  # Bar numbers

  getBarNumbers <- reactive({
    validate(need(!is.null(input$bar_numbers), "Waiting for bar numbers"))
    input$bar_numbers
  })


  data <- reactive({
    acresultdata()
  })

  # Update choices for selectInput based on column names of the data
  observe({
    list1 <- names(acresultdata())
    selected_indices1<- 7:length(list1)
    selected_colnames1 <- list1[selected_indices1]
    exclude_columns1 <- c( "reads", "read" )
    selected_colnames1 <- setdiff(selected_colnames1, exclude_columns1)
    updateSelectInput(session, "selectedColumnOverlap",
      choices = selected_colnames1
    )
  })

  

  processed_upset <- reactiveVal()

  observeEvent(input$applyButtonOverlap, {
    req(input$selectedColumnOverlap)
    group_sym <- rlang::sym(input$selectedColumnOverlap)

    df_temp <- data() %>%
      distinct(junction, !!group_sym) %>%
      group_by(junction, !!group_sym) %>%
      tally() %>%
      tidyr::spread(key = !!group_sym, value = n, fill = 0)

    processed_upset(df_temp)
  })


  # Read input data

  getValidSets <- reactive({
    withProgress(message = "Deriving input sets", value = 0, {
      setdata <- processed_upset()
      logical_cols <-
        colnames(setdata)[apply(setdata, 2, function(x) {
          all(x %in% c(0, 1))
        })]
      names(logical_cols) <- logical_cols
      lapply(logical_cols, function(x) {
        which(setdata[[x]] == 1)
      })
    })
  })


  # Subset sets to those selected

  getSelectedSets <- reactive({
    valid_sets <- getValidSets()
    chosen_sets <- getSelectedSetNames()
    sets <- valid_sets[chosen_sets]
    if (getSetSort()) {
      sets <- sets[order(unlist(lapply(sets, length)))]
    }
    sets
  })

  # Get the sets we're going to use based on nsets

  getSets <- reactive({
    selected_sets <- getSelectedSets()
    req(length(selected_sets) > 0)

    nsets <- getNsets()
    selected_sets[1:min(nsets, length(selected_sets))]
  })

  # Calculate intersections between sets

  calculateIntersections <- reactive({
    selected_sets <- getSets()

    withProgress(message = "Calculating set intersections", value = 0, {
      sets <- getSets()
      nsets <- length(sets)

      # Get all possible combinations of sets

      combinations <- function(items, pick) {
        x <- combn(items, pick)
        lapply(seq_len(ncol(x)), function(i) {
          x[, i]
        })
      }

      assignment_type <- getIntersectionAssignmentType()

      # No point starting at size 1 in a non-upset plot

      startsize <- ifelse(assignment_type == "upset", 1, 2)

      combos <- lapply(startsize:nsets, function(x) {
        combinations(1:length(selected_sets), x)
      })

      # Calculate the intersections of all these combinations

      withProgress(message = "Running intersect()", value = 0, {
        intersects <- lapply(combos, function(combonos) {
          lapply(combonos, function(combo) {
            Reduce(intersect, selected_sets[combo])
          })
        })
      })

      # For UpSet-ness, membership of higher-order intersections takes priority Otherwise just return the number of entries in each intersection

      intersects <- lapply(1:length(intersects), function(i) {
        intersectno <- intersects[[i]]
        members_in_higher_levels <-
          unlist(intersects[(i + 1):length(intersects)])
        lapply(intersectno, function(intersect) {
          if (assignment_type == "upset") {
            length(setdiff(intersect, members_in_higher_levels))
          } else {
            length(intersect)
          }
        })
      })

      combos <- unlist(combos, recursive = FALSE)
      intersects <- unlist(intersects)

      if (!getShowEmptyIntersections()) {
        combos <- combos[which(intersects > 0)]
        intersects <- intersects[which(intersects > 0)]
      }

      # Sort by intersect size

      combos <- combos[order(intersects, decreasing = TRUE)]
      intersects <-
        intersects[order(intersects, decreasing = TRUE)]
      list(combinations = combos, intersections = intersects)
    })
  })

  # Add some line returns to contrast names

  getSetNames <- reactive({
    selected_sets <- getSets()
    gsub("_", " ", names(selected_sets))
  })

  # Make the grid of points indicating set membership in intersections

  upsetGrid <- reactive({
    selected_sets <- getSets()
    ints <- calculateIntersections()

    intersects <- ints$intersections
    combos <- ints$combinations

    # Reduce the maximum number of intersections if we don't have that many

    nintersections <- getNintersections()
    nintersections <- min(nintersections, length(combos))

    # Fetch the number of sets

    nsets <- getNsets()
    setnames <- getSetNames()

    lines <-
      data.table::rbindlist(lapply(1:nintersections, function(combono) {
        data.frame(
          combo = combono,
          x = rep(combono, max(2, length(combos[[combono]]))),
          y = (nsets - combos[[combono]]) + 1,
          name = setnames[combos[[combono]]]
        )
      }))

    plot_ly(
      type = "scatter",
      mode = "markers",
      marker = list(color = "lightgrey", size = 8)
    ) %>%
      add_trace(
        type = "scatter",
        x = rep(
          1:nintersections,
          length(selected_sets)
        ),
        y = unlist(lapply(1:length(selected_sets), function(x) {
          rep(x - 0.5, nintersections)
        })),
        hoverinfo = "none"
      ) %>%
      add_trace(
        type = "scatter",
        data = group_by(lines, combo),
        mode = "lines+markers",
        x = lines$x,
        y = lines$y - 0.5,
        line = list(color = "black", width = 3),
        marker = list(
          color = "black",
          size = 10
        ),
        hoverinfo = "text",
        text = ~name
      ) %>%
      layout(
        xaxis = list(
          showticklabels = FALSE,
          showgrid = FALSE,
          zeroline = FALSE
        ),
        yaxis = list(
          showticklabels = FALSE,
          showgrid = TRUE,
          range = c(0, nsets),
          zeroline = FALSE,
          range = 1:nsets
        ),
        margin = list(t = 0, b = 30)
      )
  })


  # Make the bar chart illustrating set sizes or only show set names
  upsetSetSizeBarChart <- reactive({
    setnames <- getSetNames()
    selected_sets <- getSets()

    if (input$show_set_size) {
      # Show bar chart with set sizes
      plot_ly(
        x = unlist(lapply(selected_sets, length)),
        y = setnames,
        type = "bar",
        orientation = "h",
        marker = list(color = "black")
      ) %>% layout(
        bargap = 0.4,
        yaxis = list(
          categoryarray = rev(setnames),
          categoryorder = "array"
        )
      )
    } else {
      # Only show set names as clean text entries, with no additional plot elements
      plot_ly(
        type = "scatter",
        mode = "text",
        text = "",
        y = setnames,
        ### x = array(1.5, length(setnames)),  # Center text in the middle of the plot area
        orientation = "l"
        # textposition = "middle center",
        # textfont = list(color = "black", size = 12, font="bold")
      ) %>% layout(
        yaxis = list(
          categoryarray = rev(setnames),
          categoryorder = "array",
          showticklabels = TRUE,
          tickfont = list( # Custom font style for y-axis labels
            font = "bold", # Specify font family and style
            size = 10, # Specify font size
            color = "black"
          ),
          autorange = TRUE,
          showline = FALSE,
          showgrid = FALSE,
          zeroline = FALSE
        ),
        xaxis = list(
          showticklabels = FALSE, # No x-axis tick labels
          zeroline = FALSE,
          showline = FALSE,
          showgrid = FALSE, # No grid lines
          title = "" # No x-axis title
        ),
        margin = list(l = 30), # Adjust left margin if set names are cut off
        showlegend = FALSE # No legend
      )
    }
  })

  # Make the bar chart illustrating intersect size

  upsetIntersectSizeBarChart <- reactive({
    ints <- calculateIntersections()
    intersects <- ints$intersections
    combos <- ints$combinations
    nintersections <- getNintersections()

    p <-
      plot_ly(showlegend = FALSE) %>% add_trace(
        x = 1:nintersections,
        y = unlist(intersects[1:nintersections]),
        type = "bar",
        marker = list(
          color = "black",
          hoverinfo = "none"
        )
      )

    bar_numbers <- getBarNumbers()

    if (bar_numbers) {
      p <-
        p %>% add_trace(
          type = "scatter",
          mode = "text",
          x = 1:nintersections,
          y = unlist(intersects[1:nintersections]) + (max(intersects) * 0.05),
          text = unlist(intersects[1:nintersections]),
          textfont = list(color = "black")
        )
    }

    p
  })
  


  plot_upset <- reactive({
    grid <- upsetGrid()
    set_size_chart <- upsetSetSizeBarChart()
    intersect_size_chart <- upsetIntersectSizeBarChart()


    intersect_size_chart <-
      intersect_size_chart %>% layout(yaxis = list(title = "Intersections size"))


    s1 <-
      subplot(
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        plotly_empty(type = "scatter", mode = "markers"),
        set_size_chart,
        nrows = 2,
        widths = c(0.6, 0.4)
      )
    s2 <-
      subplot(intersect_size_chart,
              grid,
              nrows = 2,
              shareX = TRUE
      ) %>% layout(showlegend = FALSE)

    subplot(s1, s2, widths = c(0.3, 0.7))
  })

  
  output$plotly_upset <- renderPlotly({
    plot_upset()  # Use the reactive expression here
  })
  
  
  output$download_up <- downloadHandler(
    filename = function() {
      paste0("Upset_plot.",input$extension_up)
    },
    
    content = function(file) {

        tempFile <- tempfile(fileext = ".html")
        p <- plot_upset()
        htmlwidgets::saveWidget(p, tempFile, selfcontained = TRUE)
        
        # Use webshot to take a screenshot of the HTML file
        webshot(tempFile, file = file, delay = 0.2, vwidth = input$user_width_up, vheight = input$user_height_up, zoom=2)
        
        
      }
    )
    

  #### ----------------------------####
  #### Network Analysis          #####
  #### ----------------------------####
  ######################## Scatter Plot  ######################
  # A reactive expression with the ggvis plot
  # Update choices for selectInput based on column names of the data
  vis <- reactive({
    # Ensure these inputs are available
    req(input$xvar, input$yvar, input$sizevar, input$colorscatter)
    xvar_name <- names(axis_vars)[axis_vars == input$xvar]
    yvar_name <- names(axis_vars)[axis_vars == input$yvar]
    sizevar_name <- names(axis_vars)[axis_vars == input$sizevar]
    
    
    data <- processed_reactive() 
    data <- data %>% dplyr::mutate(v_family = gsub("-.*", "", v_call)) %>% 
      dplyr::mutate(j_family = gsub("-.*", "", j_call))
    
    data$tooltip_text <- paste("Cluster ID: ", data$cluster, "\n",
                               "V subgroup: ", data$v_family, "\n",
                               "J subgroup: ", data$j_family, "\n",
                               xvar_name, ": ", data[[input$xvar]], "\n", 
                               yvar_name, ": ", data[[input$yvar]], "\n", 
                               sizevar_name, ": ", data[[input$sizevar]], "\n",sep = "")
    
    if ("sample" %in% names(data)) {
      data <- data %>%
        group_by(cluster) %>%
        mutate(samples = n_distinct(sample))
      palette <- brewer.pal(8, input$colorscatter)
      colorCount <- length(unique(data$samples))
      colors <- colorRampPalette(palette)(colorCount)
      data$samples <- factor(data$samples)
      p <- ggplot(data, aes_string(x = input$xvar, y = input$yvar, size = input$sizevar, fill = "samples", text = "tooltip_text", key = "cluster")) +
          geom_point(pch = 21) +
          theme_bw() + scale_fill_manual(values = colors) +
          labs(x = xvar_name, y = yvar_name, fill = "Number of samples", size = "")
    } else {
      p <- ggplot(data, aes_string(x = input$xvar, y = input$yvar, size = input$sizevar, text = "tooltip_text", key = "cluster")) +
        geom_point(pch = 21,colour = "orange") +
        theme_bw() +
        labs(x = xvar_name, y = yvar_name, size = "")
    }
    ggplotly(p, tooltip = "text")
    
  })  
  
  # Output of scatter plot 
  output$plot1 <- renderPlotly({
    vis()
  })
  
  
  # Reactive value to store data for clicked clusters
  clicked_clusters <- reactiveVal(list())
  
    observe({
      
      list=names(acresultdata())
      # Create a vector of the desired column indices
      selected_indices <- c(2, 3, 7:length(list))

      # Extract the column names based on the selected indices
      selected_colnames <- list[selected_indices]
      exclude_columns2 <- c( "reads", "read" )
  
      selected_colnames <- setdiff(selected_colnames, exclude_columns2)

      # Update the select input with the selected column names
      updateSelectInput(session, "group",
                        choices = selected_colnames)
  })
  
  
  net2 <- reactive({
    data_p <- processed_reactive() %>% select(-sequence_id, -v_call, -j_call, -junction) %>% distinct()
    data_a <- acresultdata() %>% mutate(locus = gsub("V.*", "", v_call)) %>% dplyr::mutate(v_family = gsub("-.*", "", v_call)) %>% 
      dplyr::mutate(j_family = gsub("-.*", "", j_call))
    req(input$group)
    group_value <- rlang::sym(input$group)
    palette <- brewer.pal(8, input$colors2)
    colorCount <- length(unique(data_a[[group_value]]))
    colors <- colorRampPalette(palette)(colorCount)
    color_list<-setNames(colors, unique(data_a[[group_value]]))
    # Color mapping function
    get_color <- function(value) {
      return(color_list[as.character(value)])
    }
    
    # Apply the color mapping to the data dataframe
    data_a <- data_a %>% 
      mutate(group_col = sapply(get(group_value), get_color))
    return(list(data_p = data_p, data_a = data_a))
    
  }) 
  

  # Observe plotly click event
  observeEvent(event_data("plotly_click"), {
    click_info <- event_data("plotly_click")
    if (!is.null(click_info)) {
      net_data <- net2()
      cluster_data <- net_data$data_p[net_data$data_p$cluster == click_info$key,]
      # Append the new plot data to the list
      clicked_clusters(c(clicked_clusters(), list(cluster_data)))
    }
  })
  
  
  # Generate and display additional plots
  output$additionalPlots <- renderUI({
    plot_output_list <- lapply(seq_along(clicked_clusters()), function(i) {
      plotName <- paste("plot_cluster", i, sep = "_")
      seqlogoplotName <-  paste("seqlogoplot_cluster", i, sep = "_")
      titleName <- paste("title_cluster", i, sep = "_")
      list(
        textOutput(titleName),  # Title element
        plotOutput(plotName, height = input$height1, width = input$width1 ), # Plot element
        plotOutput(seqlogoplotName, height = input$height2, width = input$width2)# Seqlogo Plot element
      )
    })
    
    # Ensure plot_output_list is not empty
    if (length(plot_output_list) > 0) {
      rows <- lapply(seq(1, length(plot_output_list), by = 2), function(i) {
        plot_elements <- plot_output_list[i:min(i+1, length(plot_output_list))]
        fluidRow(column(6, plot_elements[[1]]),
                 if (length(plot_elements) > 1) column(6, plot_elements[[2]]) else NULL)
      })
      
      do.call(tagList, rows)
    } else {
      # Handle the case where there are no plots
      return(NULL)
    }
  
  })
  
  # Render each plot in the list
  observe({
    for (i in seq_along(clicked_clusters())) {
      local({
        plot_index <- i
        plot_data <- clicked_clusters()[[plot_index]]
        title_output_id <- paste("title_cluster", plot_index, sep = "_")
        plot_output_id <- paste("plot_cluster", plot_index, sep = "_")
        seqlogoplot_output_id <- paste("seqlogoplot_cluster", plot_index, sep = "_")
        
        # Set the title for each plot
        output[[title_output_id]] <- renderText({
          cluster_id <- unique(plot_data$cluster)
          paste("Cluster ", cluster_id, ":")
        })
        output[[plot_output_id]] <- renderPlot({
          if (!is.null(plot_data) && nrow(plot_data) > 0) {
            net_data <- net2()
            net_data_all<-net_data$data_a
            row.names(net_data_all) <- NULL
            plot_network(plot_data, 
                         net_data_all, 
                                group = input$group,
                                color_col = "group_col",
                                scale = input$scalesize, 
                                circle_color = "grey", 
                                layout = input$layouts, 
                                legend = "on",
                                seed=input$seed,
                                text_size=input$legendfontsize)
            
          }
        })
        output[[seqlogoplot_output_id]] <- renderPlot({
          if (!is.null(plot_data) && nrow(plot_data) > 0) {
            net_data <- net2()
            net_data_all<-net_data$data_a
            row.names(net_data_all) <- NULL
            plot_seqlogo(plot_data,net_data_all, input$seqlogofontsize)

          }
        })
      })
    }
  })
  
  
  # Observe event for reset button
  observeEvent(input$resetnet, {
    clicked_clusters(list())  
  })
  
  # Observe event for undo button
  observeEvent(input$undo, {
    current_list <- clicked_clusters()
    if (length(current_list) > 0) {
      clicked_clusters(head(current_list, -1))  # Remove the last cluster
    }
  }
  )
  
######################## Grouped Cluster  ######################
  
  spread <- reactive({
    small_set <- processed_reactive() %>% select(-sequence_id, -v_call, -j_call) %>% distinct()
    p<-ggplot(small_set, aes(log10(seqs), spread, text = paste("Cluster ID: ", small_set$cluster, "\n",
                                                               "Spread: ", small_set$spread , "\n",
                                                               "Cluster size (in log10): ", round(log10(small_set$seqs),2), "\n") , 
                             key = "cluster")) + 
      geom_point(color = "#008000", alpha=0.2)+
      labs(
        title = "",
        x = "Cluster size (log10) ",
        y = "Spread",
      ) +
      theme_bw()
  
    ggplotly(p, tooltip = "text") 
      
  })
  
  
  
  # Output of scatter plot 
  output$groupSpread <- renderPlotly({
    spread()
  })
  
  
  
  

  net3 <- reactive({
    data_p <- processed_reactive() %>% select(-sequence_id, -v_call, -j_call) %>% distinct()
    
    data_p <- data_p[data_p$spread >= input$spread_lower & data_p$spread <= input$spread_upper, ]
    data_p <- data_p[log10(data_p$seqs) >= input$size_lower & log10(data_p$seqs) <= input$size_upper, ]

    data_p$groupname <-"Group"
    
    colorCount <- length(unique(data_p$cluster))
    palette <- brewer.pal(8, input$colors3)
    colors <- colorRampPalette(palette)(colorCount)
    color_list<- setNames(colors, unique(data_p$cluster))

    
    data_p <- data_p %>%
      dplyr::mutate(cluster = as.character(cluster)) %>%
      left_join(data.frame(cluster = names(color_list), cluster_col = color_list), by = "cluster") %>%
      dplyr::mutate(cluster = as.double(cluster)) 
  }) 
  


    

  
  group_output <- eventReactive(input$goButton, {
    net_data3 <- net3()
    cluster_number <- length(unique(net_data3$cluster))
    showModal(modalDialog(
      title = paste( "Please Wait", cluster_number, "Clusters Detected"),
      "Generating the network plot. Please wait...",
      easyClose = FALSE,
      footer = NULL
    ))
    
    
     
    plot_network(net_data3, 
                           net_data3, 
                           group = "groupname",
                           color_col = "cluster_col",
                           scale = input$scalesize3, 
                           circle_color = "red", 
                           layout = input$layouts3, 
                           legend = "off",
                           seed = input$seed3)
    
    removeModal()
  })
  
  output$groupVis <- renderPlot(
    group_output(),
    width = function()
      input$width3,
    height = function()
      input$height3
    )
  
  

  
  # store some values
  store <- reactiveValues(dname="Grouped_Cluster")
  
  # create filename
  fn_downloadname <- reactive({
    
    if(input$fformat=="png") filename <- paste0(store$dname,".png",sep="")
    if(input$fformat=="tiff") filename <- paste0(store$dname,".tif",sep="")
    if(input$fformat=="jpeg") filename <- paste0(store$dname,".jpg",sep="")
    if(input$fformat=="pdf") filename <- paste0(store$dname,".pdf",sep="")
    return(filename)
  })  
  

  
  # download function
  fn_download <- function()
  {
    
    fheight <- input$fheight
    fwidth <- input$fwidth
    
    if(input$fformat=="pdf") fheight <- round(fheight*0.3937,2)
    if(input$fformat=="pdf") fwidth <- round(fwidth*0.3937,2)
    
    if(input$fformat=="png") png(fn_downloadname(), height=fheight, width=fwidth, res=330, units="cm")
    if(input$fformat=="tiff") tiff(fn_downloadname(), height=fheight, width=fwidth, res=330, units="cm",compression="lzw")
    if(input$fformat=="jpeg") jpeg(fn_downloadname(), height=fheight, width=fwidth, res=330, units="cm",quality=100)
    if(input$fformat=="pdf") pdf(fn_downloadname(), height=fheight, width=fwidth)
    net_data3 <- net3() 
    plot_network(net_data3, 
                 net_data3, 
                 group = "groupname",
                 color_col = "cluster_col",
                 scale = input$scalesize3, 
                 circle_color = "red", 
                 layout = input$layouts3, 
                 legend = "off",
                 seed = input$seed3)
    dev.off()
  } 
  

  # download handler
  output$bn_download <- downloadHandler(
    filename = fn_downloadname,
    content = function(file) {
      fn_download()
      file.copy(fn_downloadname(), file, overwrite=T)
    }
  )

######################## Interactive Visualization II  ######################
  rv <- reactiveValues()

  output$networkVis1 <- renderVisNetwork({
    net_data <- net2()
    plot_network_vis1(net_data$data_p[net_data$data_p$cluster == input$vis_cluster,],
                      net_data$data_a)

  })
  
  # Attribute the input value to the reactive variable
  observeEvent(input$selected_and_nearest_nodes, {
    rv$data <- input$selected_and_nearest_nodes
  })
  
  
  # filter based on reactive variable
  rt<-reactive({
    req(net2(), net2()$data_a, input$vis_cluster)
    net_data <- net2()
    edges = net_data$data_a[net_data$data_a$cluster == input$vis_cluster,]
    rownames(edges) <- NULL
    edges$group_col <- NULL
    edges
    
  })
  

  
  nearest_node_csv <- reactive({rt()})
  
  output$download_trend <- downloadHandler(filename =  'nearest_to_center_nodes_data.csv',
                                           content = function(file) {write.csv(nearest_node_csv(), file, row.names=FALSE)})
  

  output$tbl <- renderDT({
      rt()
  })
  
  
  output$tbl <- DT::renderDataTable({
    req(rt())
    datatable(rt(), options = list(pageLength = 6,scrollX = TRUE))
  })
  
  
  observe({
    
    list=names(acresultdata())
    # Create a vector of the desired column indices
    selected_indices <- c(2, 3, 7:length(list))
    
    # Extract the column names based on the selected indices
    selected_colnames <- list[selected_indices]
    exclude_columns2 <- c( "reads", "read" )
    
    selected_colnames <- setdiff(selected_colnames, exclude_columns2)
    
    # Update the select input with the selected column names
    updateSelectInput(session, "vis_group",
                      choices = selected_colnames)
  })
  
  
  
  
  net <- reactive({
    data_p <- processed_reactive() %>% select(-sequence_id, -v_call, -j_call, -junction) %>% distinct()
    data_a <- acresultdata() %>% mutate(locus = gsub("V.*", "", v_call))
    req(input$vis_group)
    group_value <- rlang::sym(input$vis_group)
    palette <- brewer.pal(8, input$vis_colors)
    colorCount <- length(unique(data_a[[group_value]]))
    colors <- colorRampPalette(palette)(colorCount)
    color_list<-setNames(colors, unique(data_a[[group_value]]))
    # Color mapping function
    get_color <- function(value) {
      return(color_list[as.character(value)])
    }
    
    # Apply the color mapping to the data dataframe
    data_a <- data_a %>% 
      mutate(group_col = sapply(get(group_value), get_color))
    return(list(data_p = data_p, data_a = data_a))
    
  }) 
  
  output$networkVis2 <- renderVisNetwork({
    net_data <- net()
    plot_network_vis2(net_data$data_p[net_data$data_p$cluster == input$vis_cluster,],
                      net_data$data_a,
                      layout = input$vis_layout,
                      group = input$vis_group,
                      color_col = "group_col")
  })

  
  #### ----------------------------####
  #### Phylogenetic Analysis      #####
  #### ----------------------------####
  ######################## Phylogenetic Plot ######################  

  observe({
    
    list=names(acresultdata())
    # Create a vector of the desired column indices
    selected_indices <- c(2, 3, 7:length(list))
    
    # Extract the column names based on the selected indices
    selected_colnames <- list[selected_indices]
    exclude_columns2 <- c( "reads", "read" )
    
    selected_colnames <- setdiff(selected_colnames, exclude_columns2)
    
    # Update the select input with the selected column names
    updateSelectInput(session, "phy_group",
                      choices = selected_colnames)
  })
  
  
  

  
  # filter based on reactive variable
  tree_output <- reactive({

    data_phy <- acresultdata() %>% mutate(locus = gsub("V.*", "", v_call)) %>% dplyr::mutate(v_family = gsub("-.*", "", v_call)) %>% 
      dplyr::mutate(j_family = gsub("-.*", "", j_call))
    req(input$phy_group)
    group_value <- rlang::sym(input$phy_group)
    palette <- brewer.pal(8, input$phy_color)
    colorCount <- length(unique(data_phy[[group_value]]))
    colors <- colorRampPalette(palette)(colorCount)
    color_list<-setNames(colors, unique(data_phy[[group_value]]))
    # Color mapping function
    get_color <- function(value) {
      return(color_list[as.character(value)])
    }
    
    # Apply the color mapping to the data dataframe
    net_phy_data <- data_phy %>% 
      mutate(group_col = sapply(get(group_value), get_color))
    
    
    phy_data <- net_phy_data[net_phy_data$cluster == input$phy_cluster,]
    
    pure_data <- phy_data %>% distinct(cluster, junction)
    
    seqs=AAStringSet(pure_data$junction)
    names(seqs)<-paste(1:(length(pure_data$junction)))
    info<- as_tibble(phy_data)
    aligned_seqs <- AlignSeqs(seqs)
    seqs_AAbin <- as.AAbin(aligned_seqs)
    
    id_mapping <- data.frame(node = as.numeric(names(aligned_seqs)),
                              junction = gsub("-","", as.character(aligned_seqs)))
    
    phy_group_value <- rlang::sym(input$phy_group)
    
    id_mapping2 <- id_mapping %>%
      merge(info, by = "junction") %>%
      group_by(junction) %>%
      dplyr::mutate(sum_histo = n_distinct(!!phy_group_value)) %>%
      dplyr::mutate(
        Label = ifelse(
          sum_histo == 1,
          (!!phy_group_value),
          paste(unique(!!phy_group_value), collapse = "="))) %>%
      dplyr::mutate(Label2 = ifelse(sum_histo != 1, "public", !!phy_group_value)) %>%
      ungroup() %>%
      dplyr::mutate(color2 = ifelse(Label2 == "public", "black", group_col)) %>%
      distinct(node, junction, Label2, Label, color2) %>%
      arrange(node)
    
    names(seqs_AAbin) <- id_mapping2$Label

    aligned_seqs_phyDat <- as.phyDat(as(aligned_seqs, "matrix"), type = "AA")


    distmat <- dist.hamming(aligned_seqs_phyDat)


    tree <- NJ(distmat)


    tidytree <- as.treedata(tree)

    tree$tip.label <- as.character(id_mapping2$Label2)


    groupInfo <- split(tree$tip.label,substr(tree$tip.label, start=1, stop=80))
    new_tree3 <- groupOTU(tree, groupInfo)

    color_map <- unique(id_mapping2[, c('Label2', 'color2')])
    color_vector <- setNames(color_map$color2, color_map$Label2)
    
    unique_ids <- paste(id_mapping2$junction, seq_along(id_mapping2$Label2), sep="_")

    names(seqs_AAbin) <- unique_ids

    new_tree3$tip.label <- unique_ids
    

    p <- ggtree(new_tree3, size=input$phy_size, aes(color=group), layout=input$phy_layout) %<+% id_mapping2 + 
      scale_colour_manual(values = color_vector) + 
      geom_tiplab(aes(label=Label),size=input$phy_node_size) + theme(legend.position = "none")

    msaplot(p, seqs_AAbin, offset=input$msa_off, width=input$msa_width)
    
    
  })

  output$phy_plot <-
    renderPlot(
      tree_output(),
      width = function()
         input$width_phy, 
       height = function()
         input$height_phy
    )
  
  
  output$download_phy <- downloadHandler(
    filename = function() {
      paste0("Tree Plot.", input$extension_phy)
    },
    
    content = function(file) {
      ggsave(
        file,
        plot =  tree_output(),
        device = input$extension_phy,
        width = input$user_width_phy,
        height = input$user_height_phy,
        units = "cm",
        dpi=330,
        bg="white"
      )
    }
    
  )
  

  
  
}