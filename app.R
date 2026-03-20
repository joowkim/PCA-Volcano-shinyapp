library(shiny)
library(ggplot2)
library(ggrepel)
library(readr)

guess_col <- function(names, patterns) {
  matches <- grep(paste(patterns, collapse = "|"), names, ignore.case = TRUE, value = TRUE)
  if (length(matches) > 0) matches[1] else names[1]
}

ui <- fluidPage(
  titlePanel("PCA & Volcano Plot"),
  tabsetPanel(

    # ── Volcano tab ────────────────────────────────────────────────────────────
    tabPanel("Volcano Plot",
      sidebarLayout(
        sidebarPanel(
          actionButton("load_toy_volcano", "Load toy data", class = "btn-sm btn-default"),
          hr(),
          fileInput("file", "Upload CSV or TSV file",
                    accept = c(".csv", ".tsv", ".txt")),
          uiOutput("col_selectors"),
          numericInput("pval_cutoff", "P-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput("lfc_cutoff", "log2FC cutoff (±)", value = 1.0, min = 0, step = 0.1),
          numericInput("top_n_labels", "Label top N genes (0 = none)", value = 10, min = 0, step = 1),
          textInput("volcano_title", "Plot title", value = ""),
          downloadButton("download_volcano", "Download PNG"),
          downloadButton("download_volcano_svg", "Download SVG")
        ),
        mainPanel(
          plotOutput("volcano", height = "600px")
        )
      )
    ),

    # ── PCA tab ────────────────────────────────────────────────────────────────
    tabPanel("PCA Plot",
      sidebarLayout(
        sidebarPanel(
          actionButton("load_toy_pca", "Load toy data", class = "btn-sm btn-default"),
          hr(),
          fileInput("pca_expr_file", "Upload expression matrix CSV",
                    accept = ".csv"),
          helpText("Rows = features (genes), columns = samples. First column = feature IDs."),
          fileInput("pca_meta_file", "Upload metadata CSV",
                    accept = ".csv"),
          helpText("Must contain a column of sample IDs matching expression matrix column names."),
          uiOutput("pca_controls"),
          textInput("pca_title", "Plot title", value = ""),
          downloadButton("download_pca", "Download PNG"),
          downloadButton("download_pca_svg", "Download SVG")
        ),
        mainPanel(
          uiOutput("pca_error_ui"),
          plotOutput("pca_plot", height = "600px")
        )
      )
    )
  )
)

server <- function(input, output, session) {

  # ── Volcano ──────────────────────────────────────────────────────────────────

  data <- reactiveVal(NULL)

  # Load from file upload
  observeEvent(input$file, {
    ext <- tools::file_ext(input$file$name)
    if (ext == "tsv" || ext == "txt") {
      data(read_tsv(input$file$datapath, show_col_types = FALSE))
    } else {
      data(read_csv(input$file$datapath, show_col_types = FALSE))
    }
  })

  # Load toy DE results bundled with the app
  observeEvent(input$load_toy_volcano, {
    data(read_csv("toy_de_results.csv", show_col_types = FALSE))
  })

  output$col_selectors <- renderUI({
    req(data())
    nms <- names(data())
    lfc_default   <- guess_col(nms, c("log2FoldChange", "log2FC", "logFC", "lfc"))
    pval_default  <- guess_col(nms, c("padj", "FDR", "adj.*p", "pvalue", "pval", "p\\.value", "p_value"))
    label_default <- guess_col(nms, c("gene", "symbol", "name", "id", "feature"))
    tagList(
      selectInput("lfc_col",   "log2FC column",     choices = nms, selected = lfc_default),
      selectInput("pval_col",  "P-value column",    choices = nms, selected = pval_default),
      selectInput("label_col", "Gene name column",  choices = nms, selected = label_default)
    )
  })

  volcano_data <- reactive({
    req(data(), input$lfc_col, input$pval_col)
    df <- data()
    df$lfc  <- as.numeric(df[[input$lfc_col]])
    df$pval <- as.numeric(df[[input$pval_col]])
    df <- df[!is.na(df$lfc) & !is.na(df$pval) & df$pval > 0, ]
    df$neg_log10_p <- -log10(df$pval)
    df$category <- ifelse(
      df$lfc >  input$lfc_cutoff & df$pval < input$pval_cutoff, "Up",
      ifelse(df$lfc < -input$lfc_cutoff & df$pval < input$pval_cutoff, "Down", "NS")
    )
    df$category <- factor(df$category, levels = c("Up", "Down", "NS"))
    df
  })

  make_volcano <- reactive({
    req(volcano_data())
    df <- volcano_data()
    p <- ggplot(df, aes(x = lfc, y = neg_log10_p, color = category)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c(Up = "#d62728", Down = "#1f77b4", NS = "grey60"),
                         name = NULL) +
      geom_vline(xintercept = c(-input$lfc_cutoff, input$lfc_cutoff),
                 linetype = "dashed", color = "black", linewidth = 0.4) +
      geom_hline(yintercept = -log10(input$pval_cutoff),
                 linetype = "dashed", color = "black", linewidth = 0.4) +
      labs(x = input$lfc_col, y = paste0("-log\u2081\u2080(", input$pval_col, ")"),
           title = input$volcano_title) +
      theme_bw(base_size = 14)

    # Label top N genes by significance (smallest p-value), skipping if N = 0 or no label column
    n <- as.integer(input$top_n_labels)
    if (!is.null(input$label_col) && input$label_col %in% names(df) && n > 0) {
      top <- head(df[order(df$pval), ], n)
      p <- p + geom_text_repel(
        data = top,
        aes(label = .data[[input$label_col]]),
        size = 3, max.overlaps = Inf, show.legend = FALSE
      )
    }
    p
  })

  output$volcano <- renderPlot({ make_volcano() })

  output$download_volcano <- downloadHandler(
    filename = function() "volcano_plot.png",
    content  = function(file) ggsave(file, plot = make_volcano(), width = 8, height = 6, dpi = 150)
  )

  output$download_volcano_svg <- downloadHandler(
    filename = function() "volcano_plot.svg",
    content  = function(file) ggsave(file, plot = make_volcano(), width = 8, height = 6)
  )

  # ── PCA ──────────────────────────────────────────────────────────────────────

  pca_expr <- reactiveVal(NULL)
  pca_meta <- reactiveVal(NULL)

  # Load from file uploads
  observeEvent(input$pca_expr_file, {
    pca_expr(read_csv(input$pca_expr_file$datapath, show_col_types = FALSE))
  })
  observeEvent(input$pca_meta_file, {
    pca_meta(read_csv(input$pca_meta_file$datapath, show_col_types = FALSE))
  })

  # Load toy datasets bundled with the app
  observeEvent(input$load_toy_pca, {
    pca_expr(read_csv("toy_expression.csv", show_col_types = FALSE))
    pca_meta(read_csv("toy_metadata.csv",   show_col_types = FALSE))
  })

  # Validation: check that sample IDs match exactly
  pca_validation <- reactive({
    req(pca_expr(), pca_meta())
    expr     <- pca_expr()
    meta     <- pca_meta()
    # Expression matrix: first column is feature IDs, rest are sample columns
    sample_cols <- names(expr)[-1]

    if (is.null(input$pca_id_col) || !input$pca_id_col %in% names(meta)) {
      return(list(ok = FALSE, msg = "Select the sample ID column from the metadata."))
    }
    meta_ids <- as.character(meta[[input$pca_id_col]])

    missing_in_meta <- setdiff(sample_cols, meta_ids)
    missing_in_expr <- setdiff(meta_ids, sample_cols)

    if (length(missing_in_meta) > 0 || length(missing_in_expr) > 0) {
      msg <- character(0)
      if (length(missing_in_meta) > 0)
        msg <- c(msg, paste("In expression but not metadata:", paste(missing_in_meta, collapse = ", ")))
      if (length(missing_in_expr) > 0)
        msg <- c(msg, paste("In metadata but not expression:", paste(missing_in_expr, collapse = ", ")))
      return(list(ok = FALSE, msg = paste(msg, collapse = "\n")))
    }
    list(ok = TRUE, msg = NULL)
  })

  output$pca_error_ui <- renderUI({
    req(input$pca_id_col)
    v <- pca_validation()
    if (!v$ok) {
      div(style = "color: red; margin-bottom: 10px;",
          strong("Sample ID mismatch:"), br(),
          pre(v$msg))
    }
  })

  # Dynamic controls rendered after both files are uploaded
  output$pca_controls <- renderUI({
    req(pca_expr(), pca_meta())
    meta_cols <- names(pca_meta())
    n_samples  <- ncol(pca_expr()) - 1  # exclude feature ID column
    max_pc     <- min(n_samples, 10)
    pc_choices <- paste0("PC", seq_len(max_pc))

    tagList(
      selectInput("pca_id_col",    "Sample ID column (metadata)", choices = meta_cols),
      selectInput("pca_color_col", "Color points by",             choices = meta_cols),
      selectInput("pca_shape_col", "Shape points by",             choices = c("None", meta_cols)),
      selectInput("pca_x", "X axis", choices = pc_choices, selected = "PC1"),
      selectInput("pca_y", "Y axis", choices = pc_choices, selected = "PC2"),
      numericInput("pca_log2_pseudo", "log2 pseudocount (0 = no transform)", value = 0, min = 0, step = 0.01),
      checkboxInput("pca_scale",  "Scale features",    value = TRUE),
      checkboxInput("pca_labels", "Show sample labels", value = FALSE)
    )
  })

  pca_result <- reactive({
    req(pca_expr(), pca_meta(), input$pca_id_col, input$pca_color_col)
    v <- pca_validation()
    req(v$ok)

    expr <- pca_expr()
    meta <- pca_meta()

    # Matrix: drop feature ID column, transpose so rows = samples
    mat <- as.matrix(expr[, -1])
    mat <- t(mat)
    storage.mode(mat) <- "numeric"

    if (!is.null(input$pca_log2_pseudo) && input$pca_log2_pseudo > 0) {
      mat <- log2(mat + input$pca_log2_pseudo)
    }

    # Remove features with zero variance
    mat <- mat[, apply(mat, 2, var) > 0, drop = FALSE]

    pca <- prcomp(mat, center = TRUE, scale. = isTRUE(input$pca_scale))

    # Build scores data frame
    scores <- as.data.frame(pca$x)
    scores$sample_id <- rownames(scores)

    # Join metadata (reorder to match)
    meta$sample_id <- as.character(meta[[input$pca_id_col]])
    scores <- merge(scores, meta, by = "sample_id")

    # Variance explained
    var_exp <- summary(pca)$importance["Proportion of Variance", ] * 100

    list(scores = scores, var_exp = var_exp)
  })

  make_pca <- reactive({
    req(pca_result(), input$pca_x, input$pca_y, input$pca_color_col)
    res     <- pca_result()
    scores  <- res$scores
    var_exp <- res$var_exp

    x_pct <- round(var_exp[input$pca_x], 1)
    y_pct <- round(var_exp[input$pca_y], 1)

    p <- ggplot(scores, aes(x = .data[[input$pca_x]],
                            y = .data[[input$pca_y]],
                            color = .data[[input$pca_color_col]])) +
      geom_point(size = 3, alpha = 0.85) +
      labs(x = paste0(input$pca_x, " (", x_pct, "%)"),
           y = paste0(input$pca_y, " (", y_pct, "%)"),
           color = input$pca_color_col,
           title = input$pca_title) +
      theme_bw(base_size = 14)

    if (!is.null(input$pca_shape_col) && input$pca_shape_col != "None") {
      p <- p + aes(shape = .data[[input$pca_shape_col]]) +
               scale_shape_discrete(name = input$pca_shape_col)
    }

    if (isTRUE(input$pca_labels)) {
      p <- p + geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = Inf, show.legend = FALSE)
    }
    p
  })

  output$pca_plot <- renderPlot({
    req(pca_validation()$ok)
    make_pca()
  })

  output$download_pca <- downloadHandler(
    filename = function() "pca_plot.png",
    content  = function(file) ggsave(file, plot = make_pca(), width = 8, height = 6, dpi = 150)
  )

  output$download_pca_svg <- downloadHandler(
    filename = function() "pca_plot.svg",
    content  = function(file) ggsave(file, plot = make_pca(), width = 8, height = 6)
  )
}

shinyApp(ui, server)
