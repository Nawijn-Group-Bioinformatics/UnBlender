library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(MuSiC)
library(SingleCellExperiment)
library(Seurat)
library(tibble)

source("config.r")
source(paste0(basedir, "functions.r"))

set.seed(123)
# Set he maximum upload size
options(shiny.maxRequestSize = 60 * 1024^2)
# Silence warnings (I know what i'm doing)
options(shiny.devmode = FALSE)
# Read the down sampled dat structure
so_small <- readRDS(paste0(datadir, "so_downs_5000.rds"))
# Read the files with tissue types for the bulk RNA
tissue_types <- list.files(paste0(datadir, "pseudobulks/"))
# Read a simplefied version the cell annotation dat for building the menu
cell_annotations <- read.table(
  paste0(datadir, "cellannotations.txt"),
  h = T,
  sep = "\t"
)

# Create some default data sets with the associated clusters
mytissue <- "parenchyma"
# A set of immune cells
immune <- cell_annotations %>%
  filter(level1 == "Immune", arc == mytissue) %>%
  pull(level4)
# A set of epithelial cells
epithelial <- cell_annotations %>%
  filter(level1 == "Epithelial", arc == mytissue) %>%
  pull(level4)


myclusters_org <- list(
  "Immune" = immune,
  "Epithelial" = epithelial
)


##### SERVER #####

shinyServer(function(input, output, session) {
  user_data_real <- reactiveValues(
    collections = NULL, # Will be filled with cell types or the merges
    geneset = NULL,
    tissue_type = NULL,
    pseudobulk_number = NULL, # Will be filled with the number of pseudobulks used
    music_error = "",
    deconv_error = "",
    upload_error = "No data uploaded (yet).",
    music_results = NULL,
    eval_results = NULL
  )
  # The following can be used for demo purposes
  user_data_demo <- reactiveValues(
    collections = myclusters_org,
    geneset = NULL,
    tissue_type = mytissue,
    #  music_results = readRDS(paste0(datadir,"music_output_processed.rds")),
    music_error = "",
    #  eval_results = readRDS(paste0(datadir,"cached_result_ground_thruth_eval.rds"))
  )

  session$onSessionEnded(function() {
    message("Session ended; stopping app.")
    stopApp()
  })

  user_data <- user_data_real
  if (DEVELOPMENT == TRUE) {
    user_data <- user_data_demo
  }

  ##### START PAGE #########
  output$intro <- renderUI({
    includeMarkdown(file.path(basedir, "www/intro.md"))
  })

  observeEvent(input$goto_start, {
    updateTabItems(session, "tabs", "select_tissue")
  })

  ##### BUILDING COLLECTIONS #####

  #### Select tissue ####
  observeEvent(input$tissue_type, {
    # Check whether there is a an input for tissue tyoe
    req(input$tissue_type)
    user_data$tissue_type <- input$tissue_type
    user_data$tissue_name <- names(tissues)[tissues == user_data$tissue_type]
    # Reset all collections
    user_data$collections <- NULL
  })

  observeEvent(input$goto_select_celltypes, {
    updateTabItems(session, "tabs", "create_sets")
  })

  #### Select cells ####
  observeEvent(input$goto_select_tissuetypes, {
    updateTabItems(session, "tabs", "select_tissue")
  })

  observeEvent(input$goto_select_removedegs, {
    updateTabItems(session, "tabs", "remove_degs")
  })

  output$no_tissue_error <- renderUI({
    x <- ""
    if (is.null(user_data$tissue_type)) {
      x <- no_tissue_selected_error()
    }
    # Retrn the non-empty erromesasage is there is no tissue selected
    x
  })

  output$current_tissue <- renderUI({
    req(user_data$tissue_type)
    HTML(paste0("Current tissue : <b>", user_data$tissue_name))
  })

  output$tissue_confirm <- renderUI({
    req(user_data$tissue_type)
    HTML(paste0("Current tissue : <b>", user_data$tissue_name, "</b>"))
  })

  output$mytree <- renderTree({
    req(user_data$tissue_type)
    mytissue <- user_data$tissue_type
    # copy ifelse logic from slide Tessa
    if (mytissue == "parenchyma") {
      current_tree <- cell_annotations %>% filter(arc == mytissue)
    }
    if (mytissue == "bronchial_biopsy") {
      current_tree <- cell_annotations %>%
        filter(
          sample_type == "biopsy",
          arc %in% c('airway', 'Intermediate Bronchi', "Trachea")
        )
    }
    if (mytissue == "nasal_brush") {
      current_tree <- cell_annotations %>%
        filter(
          sample_type %in% c("brush", "scraping"),
          arc %in% c('nose', 'Inferior turbinate')
        )
    }
    if (mytissue == "bronchial_brush") {
      current_tree <- cell_annotations %>%
        filter(sample_type %in% c("brush"), arc %in% c("Distal Bronchi"))
    }

    dfToTree(
      current_tree,
      hierarchy = c("level1", "level2", "level3", "level4")
    )
  })

  # Method to ensure that the selected cells from the tree are added as individual collections
  observeEvent(input$mytree, {
    myids <- get_selected(input$mytree, format = "names") %>% unlist()
    myids <- myids[myids %in% (cell_annotations %>% pull(level4))] %>% unique()
    # Only proceed if any IDS are left
    user_data$collections <- list()
    for (myid in myids) {
      collection_name <- gsub(myid, pattern = "\\s+", replace = "_")
      user_data$collections[[collection_name]] <- myid
    }
  })

  observeEvent(input$delete_collections, {
    print("Deleting")
    req(input$current_collections_rows_selected)
    print("Deleting 2")

    mydat <- user_data$collection_table[
      input$current_collections_rows_selected,
    ] %>%
      pull(collection)
    print("Deleting")

    for (mycollection in mydat) {
      user_data$collections[[mycollection]] <- NULL
    }
  })

  observeEvent(input$merge_collections, {
    print("Merging")
    req(input$current_collections_rows_selected)
    req(input$current_collections_rows_selected)
    print("Merging 2")

    mycollections <- user_data$collection_table[
      input$current_collections_rows_selected,
    ] %>%
      pull(collection)
    print("Merging")

    user_data$collections[[input$new_collection_name]] <- mycollections

    for (mycollection in mycollections) {
      user_data$collections[[mycollection]] <- NULL
    }
  })

  output$current_collections <- DT::renderDataTable({
    mylist <- list()
    for (myname in names(user_data$collections)) {
      mylist[[myname]] <- tibble(
        "collection" = myname,
        "cells" = paste(user_data$collections[[myname]], collapse = ", ")
      )
    }
    mydf <- bind_rows(mylist)
    user_data$collection_table <- mydf
    DT::datatable(mydf, 
      rownames = F, 
      filter = "top",
      options = list(pageLength = 100))
  }
)

  #### Remove DEGs ####

  observeEvent(input$remove_degs, {
    #     output$unmatched_degs <- renderText({
    input_degs <- unlist(strsplit(input$inserted_degs, split = "\n"))
    input_degs <- input_degs %>% gsub(pattern = " ", replace = "")
    geneset <- row.names(so_small)

    matched <- input_degs[(input_degs %in% geneset)]
    user_data$matched_DEGS <- matched

    # kan efficienter, door de gevonden met de input te vergelijken en NOT
    unmatched <- input_degs[!(input_degs %in% geneset)]

    str_unmatched <- paste0(unmatched, collapse = "<br>")
    output$unmatched_degs <- renderUI(HTML(str_unmatched))
    output$matched_degs <- renderUI(HTML(paste0(matched, collapse = "<br>")))
  })

  # observeEvents(input$goto_select_validate, {
  #   updateTabItems(session, "tabs", "current_sets")
  # })

  observeEvent(input$goto_select_validate, {
    updateTabItems(session, "tabs", "current_sets")
  })
  ##### EVALUATING COLLECTIONS  ######

  error_check <- reactive({
    mystring <- ""
    if (is.null(user_data$tissue_type)) {
      mystring <- paste0(mystring, no_tissue_selected_error())
    }
    if (is.null(user_data$collections)) {
      mystring <- paste0(mystring, no_cells_selected_error())
    }
    HTML(mystring)
  })

  output$pre_validation_error <- renderUI({
    x <- error_check()
  })

  output$pre_validation_summary <- renderUI({
    HTML(
      paste0("<b>Tissue type:</b> ", user_data$tissue_type),
      '<br>',
      "<b>Cell collections: </b>",
      selected_genes_formatter(user_data$collections)
    )
  })

  output$correlation_warning_text <- renderUI({
    tp <- get_correlations_per_celltype(user_data$eval_results$corr_df)
    names <- list(tp$cluster_name[is.na(tp$mycor)])
    HTML(
      paste0( paste0(unlist(names), collapse = ", "), " are missing due to a standard deviation of 0" )
    )
  })

  output$evaluate_button <- renderUI({
    mydisabled <- FALSE
    if (is.null(user_data$tissue_type)) {
      mydisabled <- TRUE
    }
    if (is.null(user_data$collections)) {
      mydisabled <- TRUE
    }
    # actionButton('evaluate', "Evaluate", disabled = mydisabled)
    shiny::actionButton(
      "evaluate",
      label = "Evaluate cell collections",
      icon = icon("check"),
      class = "btn-primary",
      style = "color: white; background-color: #28a745; border-color: #28a745;",
      disabled = mydisabled
      )
  })

  run_music_for_gt <- function(so_small, tissue_type) {

    so_small_sub <- create_subset_so(so_small, tissue_type)

    myclusters <- user_data$collections %>%
      stack() %>%
      dplyr::rename(cluster_member = values, cluster_name = ind) %>%
      dplyr::mutate(across(everything(), as.character)) %>%
      distinct()
    #
    gt <- create_ground_truth(so = so_small_sub, user_clusters = myclusters)

    message("Done creating ground truth")
    message("Subsetting Single Cell based on cell selection")
    lung.sce <- as.SingleCellExperiment(so_small_sub)
    new_sce <- create_input_sce(start_sce = lung.sce, new_clusters = myclusters)
    message("Done creating Single Cell object")

    message("Reading Pseudobulk")
    df <- fread(paste0(
      pseudobulksdir,
      "/pseudobulks_",
      user_data$tissue_type,
      ".csv"
    ))

    names(df)[1] <- "V1"
    lmx_bulk <- data.matrix(as_tibble(df) %>% column_to_rownames("V1"))
    #
    message("Done reading pseudobulk")
    message("Running music")
    estimated_properties <- run_music_algorithm2(
      bulk_data = lmx_bulk,
      sce_object = new_sce,
      celltypes = c("other", myclusters %>% pull(cluster_name) %>% unique())
    )
    message("Done running music")
    list("ground_truth" = gt, "music_results" = estimated_properties)
  }

  observeEvent(input$evaluate, {
    withCallingHandlers(
      {
        shinyjs::html("evaluation_progress", "")
        message("Running evaluation")
        req(user_data$collections)
        req(user_data$tissue_type)

        # This runs the actual pipeline
        myresults <- run_music_for_gt(so_small, user_data$tissue_type)
        tp <- eval_ground_truth(
          music_result = myresults[["music_results"]],
          ground_truth = myresults[["ground_truth"]]
        )

        user_data$eval_results <- tp
        message("Done")
      },
      message = function(m) {
        shinyjs::html(
          id = "evaluation_progress",
          html = paste0(m$message, '...<br>'),
          add = TRUE
        )
      }
    )
    Sys.sleep(1)
    updateTabItems(session, "tabs", "results_evaluation")
  })

  ###### Visualisation #####

  output$no_gt_eval_error <- renderUI({
    x <- ""
    if (is.null(user_data$eval_results)) {
      x <- error_message(
        "No evaluation data available.<br> Please select tissue type and cells and run an evaluation."
      )
    }
    x
  })

  # Only create control elements if there are actual data
  output$flip_mape <- renderUI({
    req(user_data$eval_results)
    checkboxInput("flip_mape", "Flip chart")
  })

  output$flip_corplot <- renderUI({
    req(user_data$eval_results)
    checkboxInput("flip_corplot", "Flip chart")
  })

  output$show_se <- renderUI({
    req(user_data$eval_results)
    checkboxInput("show_se", "Show ribbon")
  })

  output$aggregate_cor_plot <- renderPlot({
    req(user_data$eval_results)
    p <- plot_decision_cor(
      user_data$eval_results[["corr_df"]],
      flip = input$flip_corplot
    )
    p
  })

  output$aggregate_mape_plot <- renderPlot({
    req(user_data$eval_results)
    p <- plot_decision_mape(
      user_data$eval_results[["mape"]],
      flip = input$flip_mape
    )
    p
  })
  
  
  # This function handles the downloading of the Mape results
  output$download_mape <- downloadHandler(
    filename = function() {
      "UnBlender_deconv_evaluation.txt"
    },
    content = function(file) {
      # Lmerge data
      mape <- user_data$eval_results["mape"][[1]]
      corr <- user_data$eval_results["corr_df"][[1]] %>% select(c("cluster_name", "mycor", "mycor_over_outliers"))
      uni_corr <- unique(corr) %>%
                    group_by(cluster_name) %>%
                    summarise(
                      correlation = mean(mycor, na.rm=T),
                      correlation_over_outliers = mean(mycor_over_outliers)
                    )
      tp <- mape %>% left_join(uni_corr, by="cluster_name")

      # Write mape to table
      write.table(tp, file, , sep="\t", quote=FALSE, row.names = FALSE)

      # Create Analysis information (is appended to table txt file)
      file_text = "# UnBlender analysis information: "
      file_text = paste0(file_text, "\n# ",selected_genes_formatter(user_data$collections))
      nmbr_unique_pseudobulks = length(user_data$eval_results["corr_df"][[1]] %>% ungroup() %>% distinct(sample_id) %>% pull(sample_id))
      file_text = paste0(file_text, "\n# ", "Selected sample type: ", 
                          user_data$tissue_type, " (", nmbr_unique_pseudobulks, " pseudobulk samples used)")
      cat(file_text, file=file, append=T)
    }
  )

  output$gt_stats <- DT::renderDataTable({
    req(user_data$eval_results)
    mape <- user_data$eval_results[["mape"]]
    mycor <- user_data$eval_results[["corr_df"]] %>%
      dplyr::select(cluster_name, mycor, mycor_over_outliers) %>%
      summarise(
        mycor = mean(mycor, na.rm = TRUE),
        mycor_over_outliers  = mean(mycor_over_outliers)
      ) 

    toshow <- mape %>% inner_join(mycor, by = "cluster_name")
    user_data$gt_stats <- toshow

    toshow <- toshow %>%
      dplyr::rename(
        "Collection" = cluster_name,
        "Correlation" = mycor,
        "Error" = mape,
        "Correlation over Outliers" = mycor_over_outliers
      )

    DT::datatable(
      toshow,
      rownames = F,
      selection = list(mode = "single", selected = 1)
    ) %>%
      formatSignif(columns = c("Correlation", "Error"))
  })

  output$prop_error <- DT::renderDataTable({
    toshow <- user_data$eval_results[["prop_error"]] %>%
      dplyr::select(
        sample_id,
        cluster_name,
        percentage_found,
        percentage_true,
        prop_error
      )
    toshow <- toshow %>%
      dplyr::rename(
        sample = sample_id,
        collection = cluster_name,
        "True fraction" = percentage_true,
        "Estimated fraction" = percentage_found,
        "Error" = prop_error
      )

    DT::datatable(toshow, rownames = FALSE) %>%
      formatSignif(columns = c("True fraction", "Estimated fraction", "Error"))
  })

  output$mape <- DT::renderDataTable({
    toshow <- user_data$eval_results[["mape"]] 
    toshow <- toshow %>%
      dplyr::rename(
        collection = cluster_name,
        "mean proportional error" = mape,
      )
    DT::datatable(toshow, rownames = FALSE) %>%
      formatSignif(columns = c("mean proportional error"))
  })

  output$correlation_table <- DT::renderDataTable({
    to_show <- user_data$eval_results[["corr_df"]] %>%
      dplyr::select(cluster_name, mycor) %>%
      distinct()
    to_show <- to_show %>% dplyr::rename(correlation = mycor)
    DT::datatable(
      to_show,
      rownames = FALSE,
      selection = list(mode = "single", selected = 1)
    )
  })

  output$correlation_plot <- renderPlot({
    req(input$gt_stats_rows_selected)
    tp2 <- user_data$gt_stats
    collection_to_show <- tp2[input$gt_stats_rows_selected, ] %>%
      pull(cluster_name)
    tp <- tp2 %>% filter(cluster_name == collection_to_show)
    corframe <- user_data$eval_results[["corr_df"]] %>%
      filter(cluster_name == collection_to_show)
    plot_corr_df(corframe, show_se = input$show_se, show_sample_ids = FALSE)
  })

  ##### DECONVOLUTION #####

  ###### Data upload ######
  output$no_deconvolution_results_error <- renderUI({
    x <- ""
    if (is.null(user_data$music_results)) {
      x <- no_deconvolution_results()
    }
    x
  })

  output$pre_deconvolute_validation <- renderUI({
    x <- error_check()
  })

  observeEvent(input$user_bulk_file, {
    print("Getting bulk")
    user_data$upload_error <- ""
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    req(!is.null(input$user_bulk_file))
    print("Getting bulk")
    if (!is.null(input$user_bulk_file)) {
      filename <- input$user_bulk_file$datapath
      print(filename)
      myres <- grepl(filename, pattern = "\\.csv")
      if (!myres) {
        err <- 1
        user_data$upload_error <- error_message("Please upload a csv file.")
        req(err == 0)
      }
    }

    tryCatch(
      {
        df <- fread(filename)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        print("Here is the safe error")
        err <- 1
        user_data$upload_error <- "No valid fileformat detected."
        req(err == 0)
        print(safeError(e))
        print("This was the safer error")
        stop(safeError(e))
      }
    )
    names(df)[1] <- "gene"
    user_data$bulk <- df
  })

  shiny::observeEvent(input$goto_celltype_deconvolution,{
    updateTabItems(session, "tabs", "upload_data")
  })

  observeEvent(input$use_example_file, {
    if (input$use_example_file == TRUE) {
      df <- fread(paste0(datadir, "/example_gse76225.csv"))
      user_data$upload_error <- ''
      user_data$bulk <- df
    }
    if (input$use_example_file == FALSE) {
      #NO example
      user_data$bulk <- NULL
    }
  })

  output$deconvolute_button <- renderUI({
    mydisabled <- FALSE
    if (is.null(user_data$bulk)) {
      mydisabled <- TRUE
    }
    if (is.null(user_data$tissue_type)) {
      mydisabled <- TRUE
    }
    if (is.null(user_data$collections)) {
      mydisabled <- TRUE
    }
    actionButton("run_music", 
      "Deconvolute your data",
      style = "color: white; background-color: #28a745; border-color: #28a745;",
      disabled = mydisabled)
  })

  output$music_error <- renderUI({
    mymessage <- ""
    if (user_data$music_error != "") {
      mymessage <- user_data$music_error
    }
    HTML(mymessage)
  })

  output$userdata_snap <- DT::renderDataTable({
    df <- user_data$bulk
    req(!is.null(df))

    tp <- df %>% as_tibble()

    DT::datatable(
      df[1:5, 1:3],
      rownames = F,
      options = list(dom = "", ordering = F)
    ) %>%
      formatSignif(columns = c(2:3))
  })

  output$example_data <- DT::renderDataTable({
    df <- fread(paste0(datadir, "example_sample_gse76225.csv"))
    DT::datatable(
      df[1:5, 1:3],
      rownames = F,
      options = list(dom = "", ordering = F)
    ) %>%
      formatSignif(columns = c(2:3))
  })

  output$upload_error <- renderUI({
    req(is.null(user_data$bulk))
    user_data$upload_error <- 'No data uploaded'
    error_message(user_data$upload_error)
  })

  output$deconv_error <- renderUI({
    user_data$deconv_error
  })

  observeEvent(input$run_music, {
    # shinyCatch(stop("Error in running music"),  blocking_level = "error")
    user_data$deconv_error <- ""
    withCallingHandlers(
      {
        shinyjs::html("deconvolution_progress", "")
        #  shinyjs::html("text", "")
        message("Creating subset")
        # so_small_sub <- subset(so_small, subset= anatomical_region_coarse == input$tissue_type)

        ##tryCatch({
        so_small_sub <- create_subset_so(so_small, user_data$tissue_type)
        # },
        #error=function(e){
        #user_data$music_error <- "Wynand"
        # })

        message("Creating single cell object")
        sce.data2 <- as.SingleCellExperiment(so_small_sub)
        sce.data2$ann_level_4 <- gsub(
          sce.data2$ann_level_4,
          pattern = "^\\d_",
          replace = ""
        )
        #
        message("Parsing user data")

        df <- user_data$bulk
        df <- df %>%
          dplyr::filter(!gene == "") %>%
          dplyr::filter(!grepl(gene, pattern = "[/#_]")) %>%
          dplyr::filter(!is.na(gene))

        df$duprow <- duplicated(df$gene)

        df2 <- df %>% dplyr::filter(duprow == FALSE, !is.na(gene))
        lmx_bulk <- as.matrix(
          df2 %>% dplyr::select(-duprow) %>% column_to_rownames("gene")
        )

        # lmx_bulk <- data.matrix(as_tibble(df) %>% column_to_rownames("V1"))

        message("Getting clusters")
        myclusters <- user_data$collections %>%
          stack() %>%
          dplyr::rename(cluster_member = values, cluster_name = ind) %>%
          mutate(cluster_name = as.character(cluster_name))

        message("Updating single cell object")
        new_sce <- create_input_sce(
          start_sce = sce.data2,
          new_clusters = myclusters
        )
        message("Running calculation")
        tryCatch(
          {
            estimated_properties <- run_music_algorithm2(
              bulk_data = lmx_bulk,
              sce_object = new_sce,
              celltypes = c(
                "other",
                myclusters %>% pull(cluster_name) %>% unique()
              )
            )
          },
          error = function(e) {
            err <- 1
            user_data$deconv_error <- deconv_error(safeError(e))
            req(err == 0)
            stop(safeError(e))
          }
        )
        message("Saving results")
        tp <- estimated_properties
        tp$r.squared.full <- NULL
        tp$Weight.gene <- NULL
        tp1 <- lapply(tp, function(x) {
          x <- as_tibble(x, rownames = "sample_id")
        })
        all_res <- bind_rows(tp1, .id = "resname")

        tp <- all_res %>%
          filter(resname == "Est.prop.weighted") %>%
          pivot_longer(
            cols = c(-resname, -sample_id),
            names_to = "cell_type",
            values_to = "fraction"
          )

        user_data$music_results <- tp
        message("done")
        Sys.sleep(1)
        updateTabItems(session, "tabs", "user_results")
      },

      message = function(m) {
        shinyjs::html(
          id = "deconvolution_progress",
          html = paste0(m$message, '...<br>'),
          add = TRUE
        )
      }
    )
  })

  output$download_music_results <- downloadHandler(
    filename = function() {
      "deconvolution_results.xlsx"
    },
    content = function(file) {
      tp <- user_data$music_results
      toshow <- rawdata_music(tp, pivot_it = input$pivot_music)
      write_xlsx(toshow, path = file)
    }
  )

  # create heatmap of cibersort results (fractions)
  output$music_results <- DT::renderDataTable({
    req(user_data$music_results)
    tp <- user_data$music_results
    toshow <- rawdata_music(tp, pivot_it = input$pivot_music)
    mycolnames <- toshow %>% dplyr::select(!any_of(c("cell_type", "sample_id")))
    DT::datatable(
      toshow,
      rownames = F,
      filter = "top",
      options = list(scrollX = TRUE)
    ) %>%
      formatRound(columns = names(mycolnames), digits = 3)
  })

  # output$dotplot <- renderPlot({
  #   req(user_data$music_results)
  #   tp <- user_data$music_results
  #   p <- dotplot_music(tp)
  #   p
  # })

  # output$heatmap <- renderPlot({
  #   req(user_data$music_results)
  #   tp <- user_data$music_results
  #   p <- heatmap_music(tp, flipit = input$flip_heatmap)
  #   p
  # })

  output$cibersort_stackedbar <- renderPlotly({
    # stacked barchart
    req(user_data$music_results)
    tp <- user_data$music_results
    p <- stacked_bar_music(tp, flipit = TRUE)
    ggplotly(p) %>% 
      plotly::layout(
        yaxis = list(showticklabels = FALSE, ticks = "", showline = FALSE)
      )
  })
})
