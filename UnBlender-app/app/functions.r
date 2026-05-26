library(Seurat)
library(tidyr)
library(dplyr)
library(SingleCellExperiment)
library(MuSiC)
library(plotly)

##### HTML functions ######

#theme_cyberlung


selected_genes_formatter <- function(collection){
  #' This function pastes all selected genes provided in the collection together. 
  #' It formats rows with multiple genes the following way: merge_name(gene1, gene2, gene3)
  gene_list = c()
  for (i in 1:length(collection)){
      if (length(collection[[i]]) > 1){
          merged_genes = paste0(collection[[i]], collapse = ", " )
          formatted_text = paste0(names(collection[i]), "(", merged_genes, ")")
          gene_list = c(gene_list, formatted_text)
      }
      else{
          gene_list = c(gene_list, collection[[i]])
      }
  }
  gene_string = paste(gene_list, collapse = "; ")
  return (gene_string)
}

no_tissue_selected_error <- function(x) {
  error_message("Please select a tissue for building a cell collection.")
}

no_cells_selected_error <- function(x) {
  error_message(
    "Please select cells to add to your cell collection before running the evaluation."
  )
}

deconv_error <- function(x) {
  error_message(paste0(
    "The deconvolution caused an error. The error message is shown below.<br>",
    x
  ))
}
no_deconvolution_results <- function(x) {
  error_message(
    "No deconvolution results yet. Please run a deconvolution first."
  )
}

error_message <- function(x) {
  shiny::HTML(paste0('<div style="color:#FF0000;">', x, "</div>"))
}

instruction_message <- function(x) {
  shiny::HTML(paste0('<div class="instruction_message">', x, "</div>"))
}

mybox <- function(
  ...,
  title = NULL,
  footer = NULL,
  status = 'primary',
  solidHeader = TRUE,
  background = NULL,
  width = 12,
  height = NULL,
  collapsible = FALSE,
  collapsed = FALSE,
  help = NULL
) {
  boxClass <- "box"
  if (solidHeader || !is.null(background)) {
    boxClass <- paste(boxClass, "box-solid")
  }
  if (!is.null(status)) {
    boxClass <- paste0(boxClass, " box-", status)
  }
  if (collapsible && collapsed) {
    boxClass <- paste(boxClass, "collapsed-box")
  }
  if (!is.null(background)) {
    shinydashboard::validateColor(background)
    boxClass <- paste0(boxClass, " bg-", background)
  }
  style <- NULL
  if (!is.null(height)) {
    style <- paste0("height: ", validateCssUnit(height))
  }
  titleTag <- NULL
  if (!is.null(title)) {
    titleTag <- h3(class = "box-title", title)
  }
  collapseTag <- NULL
  if (collapsible) {
    buttonStatus <- status
    collapseIcon <- if (collapsed) {
      "plus"
    } else {
      "minus"
    }
    collapseTag <- div(
      class = "box-tools pull-right",
      tags$button(
        class = paste0("btn btn-box-tool"),
        `data-widget` = "collapse",
        shiny::icon(collapseIcon)
      )
    )
  }

  helpTag <- NULL
  if (!is.null(help)) {
    helpTag <- shiny::HTML(paste0(
      '<div class="pull-right" style="display:block" onClick=toggle("',
      help,
      '")>',
      '<i class="fa fa-question-circle" style="font-size:20px;" aria-hidden="false"></i></div>'
    ))
  }

  headerTag <- NULL
  if (!is.null(titleTag) || !is.null(collapseTag)) {
    headerTag <- div(class = "box-header", titleTag, HTML(helpTag), collapseTag)
  }
  div(
    class = if (!is.null(width)) {
      paste0("col-m-", width)
    },
    div(
      class = boxClass,
      style = if (!is.null(style)) {
        style
      },
      headerTag,
      div(class = "box-body", ...),
      if (!is.null(footer)) {
        div(class = "box-footer", footer)
      }
    )
  )
}


##### MUSIC FUNCTIONS #######

music_prop2 <- function(
  bulk.mtx,
  sc.sce,
  markers = NULL,
  clusters,
  samples,
  select.ct = NULL,
  cell_size = NULL,
  ct.cov = FALSE,
  verbose = TRUE,
  iter.max = 1000,
  nu = 1e-04,
  eps = 0.01,
  centered = FALSE,
  normalize = FALSE,
  ...
) {
  bulk.gene = rownames(bulk.mtx)[rowMeans(bulk.mtx) != 0]
  bulk.mtx = bulk.mtx[bulk.gene, ]
  if (is.null(markers)) {
    sc.markers = bulk.gene
  } else {
    sc.markers = intersect(bulk.gene, unlist(markers))
  }
  sc.basis = music_basis2(
    sc.sce,
    non.zero = TRUE,
    markers = sc.markers,
    clusters = clusters,
    samples = samples,
    select.ct = select.ct,
    cell_size = cell_size,
    ct.cov = ct.cov,
    verbose = verbose
  )
  cm.gene = intersect(rownames(sc.basis$Disgn.mtx), bulk.gene)
  if (is.null(markers)) {
    if (length(cm.gene) < 0.2 * min(length(bulk.gene), nrow(sc.sce))) {
      stop("Too few common genes!")
    }
  } else {
    if (length(cm.gene) < 0.2 * length(unlist(markers))) {
      stop("Too few common genes!")
    }
  }
  if (verbose) {
    message(paste("Used", length(cm.gene), "common genes..."))
  }
  m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx))
  m.bulk = match(cm.gene, bulk.gene)
  D1 = sc.basis$Disgn.mtx[m.sc, ]
  M.S = colMeans(sc.basis$S, na.rm = T)
  if (!is.null(cell_size)) {
    if (!is.data.frame(cell_size)) {
      stop(
        "cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes"
      )
    } else if (sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))) {
      stop("Cell type names in cell_size must match clusters")
    } else if (any(is.na(as.numeric(cell_size[, 2])))) {
      stop("Cell sizes should all be numeric")
    }
    my_ms_names <- names(M.S)
    cell_size <- cell_size[my_ms_names %in% cell_size[, 1], ]
    M.S <- cell_size[match(my_ms_names, cell_size[, 1]), ]
    M.S <- M.S[, 2]
    names(M.S) <- my_ms_names
  }
  Yjg = relative.ab(bulk.mtx[m.bulk, ])
  N.bulk = ncol(bulk.mtx)
  if (ct.cov) {
    Sigma.ct = sc.basis$Sigma.ct[, m.sc]
    Est.prop.allgene = NULL
    Est.prop.weighted = NULL
    Weight.gene = NULL
    r.squared.full = NULL
    Var.prop = NULL
    for (i in 1:N.bulk) {
      if (sum(Yjg[, i] == 0) > 0) {
        D1.temp = D1[Yjg[, i] != 0, ]
        Yjg.temp = Yjg[Yjg[, i] != 0, i]
        Sigma.ct.temp = Sigma.ct[, Yjg[, i] != 0]
        if (verbose) {
          message(paste(
            colnames(Yjg)[i],
            "has common genes",
            sum(Yjg[, i] != 0),
            "..."
          ))
        }
      } else {
        D1.temp = D1
        Yjg.temp = Yjg[, i]
        Sigma.ct.temp = Sigma.ct
        if (verbose) {
          message(paste(
            colnames(Yjg)[i],
            "has common genes",
            sum(Yjg[, i] != 0),
            "..."
          ))
        }
      }
      lm.D1.weighted = music.iter.ct(
        Yjg.temp,
        D1.temp,
        M.S,
        Sigma.ct.temp,
        iter.max = iter.max,
        nu = nu,
        eps = eps,
        centered = centered,
        normalize = normalize
      )
      Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp = rep(NA, nrow(Yjg))
      weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
      Weight.gene = cbind(Weight.gene, weight.gene.temp)
      r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
    }
  } else {
    Sigma = sc.basis$Sigma[m.sc, ]
    valid.ct = (colSums(is.na(Sigma)) == 0) &
      (colSums(is.na(D1)) == 0) &
      (!is.na(M.S))
    if (sum(valid.ct) <= 1) {
      stop("Not enough valid cell type!")
    }
    if (verbose) {
      message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
    }
    D1 = D1[, valid.ct]
    M.S = M.S[valid.ct]
    Sigma = Sigma[, valid.ct]
    Est.prop.allgene = NULL
    Est.prop.weighted = NULL
    Weight.gene = NULL
    r.squared.full = NULL
    Var.prop = NULL
    for (i in 1:N.bulk) {
      if (sum(Yjg[, i] == 0) > 0) {
        D1.temp = D1[Yjg[, i] != 0, ]
        Yjg.temp = Yjg[Yjg[, i] != 0, i]
        Sigma.temp = Sigma[Yjg[, i] != 0, ]
        if (verbose) {
          message(paste(
            colnames(Yjg)[i],
            "has common genes",
            sum(Yjg[, i] != 0),
            "..."
          ))
        }
      } else {
        D1.temp = D1
        Yjg.temp = Yjg[, i]
        Sigma.temp = Sigma
        if (verbose) {
          message(paste(
            colnames(Yjg)[i],
            "has common genes",
            sum(Yjg[, i] != 0),
            "..."
          ))
        }
      }
      lm.D1.weighted = music.iter(
        Yjg.temp,
        D1.temp,
        M.S,
        Sigma.temp,
        iter.max = iter.max,
        nu = nu,
        eps = eps,
        centered = centered,
        normalize = normalize
      )
      Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp = rep(NA, nrow(Yjg))
      weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
      Weight.gene = cbind(Weight.gene, weight.gene.temp)
      r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
    }
  }
  colnames(Est.prop.weighted) = colnames(D1)
  rownames(Est.prop.weighted) = colnames(Yjg)
  colnames(Est.prop.allgene) = colnames(D1)
  rownames(Est.prop.allgene) = colnames(Yjg)
  names(r.squared.full) = colnames(Yjg)
  colnames(Weight.gene) = colnames(Yjg)
  rownames(Weight.gene) = cm.gene
  colnames(Var.prop) = colnames(D1)
  rownames(Var.prop) = colnames(Yjg)
  return(list(
    Est.prop.weighted = Est.prop.weighted,
    Est.prop.allgene = Est.prop.allgene,
    Weight.gene = Weight.gene,
    r.squared.full = r.squared.full,
    Var.prop = Var.prop
  ))
}

music_basis2 <- function(
  x,
  non.zero = TRUE,
  markers = NULL,
  clusters,
  samples,
  select.ct = NULL,
  cell_size = NULL,
  ct.cov = FALSE,
  verbose = TRUE
) {
  if (!is.null(select.ct)) {
    x = x[, x@colData[, clusters] %in% select.ct]
  }
  if (non.zero) {
    x <- x[rowSums(counts(x)) > 0, ]
  }
  clusters <- as.character(colData(x)[, clusters])
  samples <- as.character(colData(x)[, samples])
  M.theta <- sapply(unique(clusters), function(ct) {
    my.rowMeans(
      sapply(unique(samples), function(sid) {
        y = counts(x)[, clusters %in% ct & samples %in% sid]
        if (is.null(dim(y))) {
          return(y / sum(y))
        } else {
          return(rowSums(y) / sum(y))
        }
      }),
      na.rm = TRUE
    )
  })
  if (verbose) {
    message("Creating Relative Abundance Matrix...")
  }
  if (ct.cov) {
    nGenes = nrow(x)
    n.ct = length(unique(clusters))
    nSubs = length(unique(samples))
    Theta <- sapply(unique(clusters), function(ct) {
      sapply(unique(samples), function(sid) {
        y = counts(x)[,
          clusters %in%
            ct &
            samples %in%
              sid
        ]
        if (is.null(dim(y))) {
          return(y / sum(y))
        } else {
          return(rowSums(y) / sum(y))
        }
      })
    })
    if (!is.null(select.ct)) {
      m.ct = match(select.ct, colnames(Theta))
      Theta = Theta[, m.ct]
    }
    Sigma.ct = sapply(1:nGenes, function(g) {
      sigma.temp = Theta[nGenes * (0:(nSubs - 1)) + g, ]
      Cov.temp = cov(sigma.temp)
      Cov.temp1 = cov(sigma.temp[
        rowSums(is.na(Theta[
          nGenes *
            (0:(nSubs - 1)) +
            1,
        ])) ==
          0,
      ])
      Cov.temp[which(colSums(is.na(sigma.temp)) > 0), ] = Cov.temp1[
        which(colSums(is.na(sigma.temp)) > 0),
      ]
      Cov.temp[, which(colSums(is.na(sigma.temp)) > 0)] = Cov.temp1[,
        which(colSums(is.na(sigma.temp)) > 0)
      ]
      return(Cov.temp)
    })
    colnames(Sigma.ct) = rownames(x)
    if (!is.null(markers)) {
      ids <- intersect(unlist(markers), rownames(x))
      m.ids = match(ids, rownames(x))
      Sigma.ct <- Sigma.ct[, m.ids]
    }
    if (verbose) {
      message("Creating Covariance Matrix...")
    }
  } else {
    Sigma <- sapply(unique(clusters), function(ct) {
      apply(
        sapply(unique(samples), function(sid) {
          y = counts(x)[,
            clusters %in%
              ct &
              samples %in%
                sid
          ]
          if (is.null(dim(y))) {
            return(y / sum(y))
          } else {
            return(rowSums(y) / sum(y))
          }
        }),
        1,
        var,
        na.rm = TRUE
      )
    })
    if (!is.null(select.ct)) {
      m.ct = match(select.ct, colnames(Sigma))
      Sigma = Sigma[, m.ct]
    }
    if (!is.null(markers)) {
      ids <- intersect(unlist(markers), rownames(x))
      m.ids = match(ids, rownames(x))
      Sigma <- Sigma[m.ids, ]
    }
    if (verbose) {
      message("Creating Variance Matrix...")
    }
  }
  S <- sapply(unique(clusters), function(ct) {
    my.rowMeans(
      sapply(unique(samples), function(sid) {
        y = counts(x)[, clusters %in% ct & samples %in% sid]
        if (is.null(dim(y))) {
          return(sum(y))
        } else {
          return(sum(y) / ncol(y))
        }
      }),
      na.rm = TRUE
    )
  })
  if (verbose) {
    message("Creating Library Size Matrix...")
  }
  S[S == 0] = NA
  M.S = colMeans(S, na.rm = TRUE)
  if (!is.null(cell_size)) {
    if (!is.data.frame(cell_size)) {
      stop(
        "cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes"
      )
    } else if (sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))) {
      stop("Cell type names in cell_size must match clusters")
    } else if (any(is.na(as.numeric(cell_size[, 2])))) {
      stop("Cell sizes should all be numeric")
    }
    my_ms_names <- names(M.S)
    cell_size <- cell_size[my_ms_names %in% cell_size[, 1], ]
    M.S <- cell_size[match(my_ms_names, cell_size[, 1]), ]
    M.S <- M.S[, 2]
    names(M.S) <- my_ms_names
  }
  D <- t(t(M.theta) * M.S)
  if (!is.null(select.ct)) {
    m.ct = match(select.ct, colnames(D))
    D = D[, m.ct]
    S = S[, m.ct]
    M.S = M.S[m.ct]
    M.theta = M.theta[, m.ct]
  }
  if (!is.null(markers)) {
    ids <- intersect(unlist(markers), rownames(x))
    m.ids = match(ids, rownames(x))
    D <- D[m.ids, ]
    M.theta <- M.theta[m.ids, ]
  }
  if (ct.cov) {
    return(list(
      Disgn.mtx = D,
      S = S,
      M.S = M.S,
      M.theta = M.theta,
      Sigma.ct = Sigma.ct
    ))
  } else {
    return(list(
      Disgn.mtx = D,
      S = S,
      M.S = M.S,
      M.theta = M.theta,
      Sigma = Sigma
    ))
  }
}


run_music_algorithm2 <- function(bulk_data, sce_object, celltypes) {
  Est.prop.paren <- music_prop2(
    bulk.mtx = bulk_data,
    sc.sce = sce_object,
    clusters = "new_clusters",
    samples = 'sample',
    select.ct = celltypes,
    verbose = F
  )
}


create_input_sce <- function(start_sce, new_clusters) {
  # Input is a tibble with cluster names and cluster_mebers
  toreplace <- setNames(
    obj = new_clusters %>% pull(cluster_name),
    nm = new_clusters %>% pull(cluster_member)
  )
  tomatch <- tibble(level4 = start_sce$ann_level_4)
  tomatch <- tomatch %>%
    dplyr::mutate(level4 = gsub(level4, pattern = "^\\d_", replace = "")) %>%
    left_join(new_clusters, by = c("level4" = "cluster_member")) %>%
    dplyr::mutate(cluster_name = as.character(cluster_name)) %>%
    dplyr::mutate(
      cluster_name = ifelse(is.na(cluster_name), "other", cluster_name)
    )

  start_sce$new_clusters <- tomatch %>% pull(cluster_name)
  start_sce
}

####### SUBSETTING ############
create_subset_so <- function(so_small, mytissue) {
  # message("Subsetting")
  #  mytissue <- user_data$tissue_type

  print(paste0("Subsetting for ", mytissue))
  #  current_tree <- cell_annotations %>% filter(tissue == input$tissue_type)
  # copy ifelse logic from slide Tessa
  if (mytissue == "parenchyma") {
    so_small_sub <- subset(
      so_small,
      subset = anatomical_region_coarse == mytissue
    )
  }
  if (mytissue == "bronchial_biopsy") {
    so_small_sub <- subset(so_small, subset = sample_type == "biopsy") %>%
      subset(
        subset = anatomical_region_coarse %in%
          c('airway', 'Intermediate Bronchi', "Trachea")
      )
  }

  if (mytissue == "nasal_brush") {
    so_small_sub <- subset(
      so_small,
      subset = sample_type %in% c("brush", "scraping")
    ) %>%
      subset(
        subset = anatomical_region_coarse %in% c('nose', 'Inferior turbinate')
      )
  }
  if (mytissue == "bronchial_brush") {
    so_small_sub <- subset(so_small, subset = sample_type %in% c("brush")) %>%
      subset(subset = anatomical_region_coarse %in% c("Distal Bronchi"))
  }

  so_small_sub
}


####### GROUND TRUTH EVALUATION  #########

create_ground_truth <- function(so, tissue, user_clusters, sample_type) {
  gt <- so@meta.data %>%
    as_tibble(rownames = "barcode") %>%
    dplyr::select(barcode, sample, anatomical_region_coarse, ann_level_4) %>%
    dplyr::mutate(across(everything(), as.character)) %>%
    #  filter(anatomical_region_coarse == tissue) %>%
    dplyr::mutate(
      ann_level_4 = gsub(ann_level_4, pattern = "^\\d_", replace = "")
    ) %>%
    dplyr::rename(sample_id = sample)

  #Update the gt to reflect the new clusters

  gt <- gt %>%
    left_join(user_clusters, by = c("ann_level_4" = "cluster_member")) %>%
    dplyr::mutate(
      cluster_name = ifelse(is.na(cluster_name), "other", cluster_name)
    )

  # Now add the count per cluster
  gt <- gt %>%
    group_by(sample_id, cluster_name) %>%
    dplyr::mutate(cluster_count = n())

  gt <- gt %>%
    ungroup() %>%
    group_by(sample_id) %>%
    dplyr::mutate(total_cells = n())

  gt <- gt %>%
    dplyr::mutate(percentage_true = cluster_count / total_cells) %>%
    dplyr::select(cluster_name, sample_id, percentage_true) %>%
    distinct()

  gt
}

eval_ground_truth <- function(music_results, ground_truth) {
  music_results <- music_results$Est.prop.weighted
  samples_pseudobulk <- row.names(music_results)
  samples_gt <- ground_truth %>% pull(sample_id) %>% unique()

  only_in_gt <- setdiff(samples_gt, samples_pseudobulk)
  only_in_pseudobulk <- setdiff(samples_pseudobulk, samples_gt)

  #We now should see recalculate the portions of clusters per sample for the groudn truth

  gt_narrow <- ground_truth %>% filter(!sample_id %in% only_in_gt)

  music_data <- music_results %>%
    as_tibble(rownames = "sample_id") %>%
    pivot_longer(
      cols = -sample_id,
      names_to = "cluster_name",
      values_to = "percentage_found"
    )

  eval_results <- music_data %>%
    full_join(
      gt_narrow,
      by = c("sample_id" = "sample_id", "cluster_name" = "cluster_name")
    ) %>%
    dplyr::mutate(
      percentage_found = ifelse(
        is.na(percentage_found),
        0.01,
        percentage_found
      ),
      percentage_true = ifelse(is.na(percentage_true), 0.01, percentage_true)
    ) %>%
    dplyr::mutate(
      prop_error = (percentage_true - percentage_found) / percentage_true
    )

  # Get mean absolute prop error

  mape <- eval_results %>%
    # filter(cluster_name!="other") %>%
    group_by(cluster_name) %>%
    summarize(mape = mean(abs(prop_error)))

  # Get the correlations
  eval_results <- eval_results %>%
    group_by(cluster_name) %>%
    dplyr::mutate(
      iqr_high = quantile(percentage_true, probs = .75) %>% as.numeric(),
      iqr_low = quantile(percentage_true, probs = .25) %>% as.numeric(),
      H = 1.5 * IQR(percentage_true)
    )
 
  # Calculate the correlation. Both over and not over the outliers.
  corr_df <- eval_results %>%
    mutate( outlier = !
      (percentage_true < iqr_high + H) & 
      (percentage_true > iqr_low - H)
      ) %>% 
    group_by(cluster_name) %>%
      dplyr::mutate(
        mycor = case_when(
          !outlier ~ floor(100 * suppressWarnings(cor(percentage_found[!outlier], percentage_true[!outlier]))) / 100
        )
      ) %>% 
    group_by(cluster_name) %>%
      dplyr::mutate(
          mycor_over_outliers = floor(100 * suppressWarnings(cor(percentage_found, percentage_true))) / 100
      )
  
  # Add if the standard deviation = 0 per cell type 
  # This causes missing points in the plot later
  no_standard_deviation <- corr_df %>%
    group_by(cluster_name) %>%
    summarise(all_same = sd(percentage_found)) %>%
    mutate(no_standard_deviation = if_else(all_same==0 | (is.na(all_same)), TRUE, FALSE)) %>%
    select(cluster_name, no_standard_deviation)
  corr_df <- dplyr::full_join(corr_df, no_standard_deviation, by="cluster_name")

  return(list(prop_error = eval_results, mape = mape, corr_df = corr_df))
}


plot_corr_df <- function(
  correlation_df,
  show_sample_ids = FALSE,
  show_se = FALSE
) {
  # Add color column
  correlation_df <- correlation_df %>%
    dplyr::select(cluster_name, percentage_true, percentage_found, mycor, outlier) %>%
    mutate(color = ifelse(outlier, "Outlier" ,"Non Outlier"))

  p <- ggplot(correlation_df, aes(x = percentage_true, y = percentage_found))
  p <- p + geom_point(size = 5, alpha = 0.6, aes(color = color))
  if (show_sample_ids == TRUE) {
    p <- p + geom_text(aes(label = sample_id), size = 5)
  }
  p <- p + geom_smooth(method = 'lm', formula = y ~ x, se = show_se)
  # Add title take first of sorted to prevent two plots one with na one with cor.
  p <- p + facet_wrap(. ~ paste0(cluster_name, " : ", sort(mycor)[1]), scales = "free")
  p <- p + labs(title = "Correlation", y = "Estimated fraction", x = "True fraction")
  
  p <- p + scale_color_manual(values = c("Non Outlier" = "#0000ff", "Outlier" = outlier_color),
                              name = "Outlier") 

  p <- p + theme_bw()
  p <- p + theme(axis.text = element_text(size = 12))
  p <- p + theme(strip.text = element_text(size = 14))
  p <- p + theme(axis.title = element_text(size = 14)) 
  p
}


get_correlations_per_celltype <- function(correlation_df){
  tp <- correlation_df %>%
    dplyr::group_by(cluster_name) %>%
    dplyr::summarise(
      avg_cor = mean(mycor, na.rm = TRUE),
      avg_cor_over_outliers = mean(mycor_over_outliers, na.rm = TRUE)
    ) %>%
    dplyr::mutate(outlier = dplyr::if_else(is.na(avg_cor), TRUE, FALSE)) %>%
    dplyr::mutate(mycor = dplyr::if_else(is.na(avg_cor), avg_cor_over_outliers, avg_cor)) %>%
    dplyr::select(cluster_name, mycor, outlier)
  return (tp)
}

plot_decision_cor <- function(correlation_df, flip = FALSE) {
  tp <- get_correlations_per_celltype(correlation_df)

  p <- ggplot(
    tp, 
    aes(
      x = reorder(cluster_name, mycor),
      y = mycor,
      color = outlier)) +
    annotate(
      geom = 'rect',
      ymin = 0.7,
      ymax = Inf,
      xmin = -Inf,
      xmax = Inf,
      fill = 'darkolivegreen2',
      alpha = 1
    ) +
    annotate(
      geom = 'rect',
      ymin = -Inf,
      ymax = 0.7,
      xmin = -Inf,
      xmax = Inf,
      fill = 'sandybrown',
      alpha = 1
    ) + 
    geom_point(size = 5) + 
    labs(x = "", y = "Correlation", color= "Outlier") +
    scale_color_manual(values = c(`FALSE` = "black", `TRUE` = outlier_color)) +
    theme_bw() +
    theme(axis.text = element_text(size = 12)) + 
    theme(strip.text = element_text(size = 14)) +
    theme(axis.title = element_text(size = 14)) + 
    coord_flip()
  p
}

plot_decision_mape <- function(mape, flip = FALSE) {
  tp <- mape %>% 
    dplyr::select(cluster_name, mape) %>% 
    distinct()

  p <- ggplot(tp, 
    aes(
      x = reorder(cluster_name, mape), 
      y = mape)) +
    ylim(0, NA) +
    annotate(
      geom = 'rect',
      ymin = -Inf,
      ymax = 1,
      xmin = -Inf,
      xmax = Inf,
      fill = 'darkolivegreen2',
      alpha = 1
    ) +
    annotate(
      geom = 'rect',
      ymin = 1,
      ymax = Inf,
      xmin = -Inf,
      xmax = Inf,
      fill = 'sandybrown',
      alpha = 1
    ) +
    geom_point(size = 5) +
    labs(x = "", y = "MAPE") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 12)) +
    theme(strip.text = element_text(size = 14)) +
    theme(axis.title = element_text(size = 14)) +
    coord_flip()
  p
}

####### MUSIC EVALUATION VISUALISATION #######
rawdata_music <- function(x, pivot_it = F) {
  x <- x %>% 
    dplyr::select(-resname)
  if (pivot_it == T) {
    x <- x %>% 
      pivot_wider(
        names_from = "cell_type", 
        values_from = "fraction")}
  x
}

dotplot_music <- function(x) {
  p <- ggplot(x, aes(x = cell_type, y = fraction)) +
    geom_point() +
    facet_wrap(. ~ sample_id) + 
    coord_flip()
  p
}


heatmap_music <- function(x, show_fractions = FALSE, flipit = FALSE) {
  p <- ggplot(
    x,
    aes(
      x = cell_type, 
      y = sample_id, 
      fill = fraction, 
      label = fraction)
  ) + geom_tile()
  if (show_fractions == T) {
    p < p + geom_text(size = 2)
  }
  p <- p + theme_bw()
  p <- p + scale_fill_gradient(high = "blue", low = "white")
  p <- p + labs(x = "", y = "", title = "Frequency distribution")
  p <- p + theme(axis.text = element_text(size = 14))
  #  p <- p + facet_wrap(.~sample_id)
  if (flipit == TRUE) {
    p <- p + coord_flip()
  }
  p <- p + theme(axis.text.x = element_text(angle = -270))
  p
}

stacked_bar_music <- function(x, flipit = FALSE) {
  p <- ggplot(x, aes(x = sample_id, y = fraction, fill = cell_type))
  p <- p + geom_bar(stat = "identity")
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = +90, hjust = 0))
  p <- p + labs(x = "", y = "", title = "Frequency distribution")
  if (flipit == TRUE) {
    p <- p + coord_flip()
  }
  p
}
