library(shiny)
library(readxl)
library(h2o)
library(ggplot2)
library(dplyr)
library(cluster)
library(kernlab)
library(MASS)
library(energy)
library(copula)
library(tibble)
library(caret)
library(scales)
library(matrixStats)

ui <- fluidPage(
  titlePanel("Campionamento AI-based con Autoencoder e Controllo Bias"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Carica file Excel", accept = ".xlsx"),
      uiOutput("varSelectUI"),
      numericInput("soglia_bias", "Soglia bias variabili chiave (%)", value = 5, min = 1, max = 20),
      actionButton("run", "Estrai campione"),
      actionButton("run2", "Forza variabili chiave (hard constraint)")
    ),
    mainPanel(
      plotOutput("biasPlot"),
      verbatimTextOutput("kOutput"),
      tableOutput("biasTable")
    )
  )
)

server <- function(input, output, session) {
  dati <- reactiveVal(NULL)
  variabili <- reactiveVal(NULL)
  report_iniziale <- reactiveVal(NULL)
  best_k_global <- reactiveVal(NULL)

  observeEvent(input$file, {
    df <- read_excel(input$file$datapath)
    df_num <- df[sapply(df, is.numeric)]
    df_num <- na.omit(df_num)
    df_num <- df_num[apply(df_num, 1, function(r) all(is.finite(r))), ]
    dati(list(df = df, df_num = df_num))

    updateSelectInput(session, "key_vars", choices = names(df_num))
  })

  output$varSelectUI <- renderUI({
    req(dati())
    selectInput("key_vars", "Scegli variabili chiave (opzionale)", 
                choices = names(dati()$df_num), multiple = TRUE)
  })

  estraiCampione <- function(df, df_num, key_vars = NULL, forza_chiave = FALSE, soglia = 5) {
    X <- scale(as.matrix(df_num))
    h2o.init()
    data_h2o <- as.h2o(X)

    latent_dims <- seq(from = floor(ncol(X) / 2), to = 10, by = -2)
    hidden_layers <- lapply(latent_dims, function(d) c(50, d, 50))
    errors <- numeric(length(latent_dims))
    features_list <- list()

    for (i in seq_along(hidden_layers)) {
      model <- h2o.deeplearning(
        x = colnames(data_h2o),
        training_frame = data_h2o,
        autoencoder = TRUE,
        hidden = hidden_layers[[i]],
        activation = "Tanh",
        epochs = 100,
        seed = 42
      )
      recon <- h2o.anomaly(model, data_h2o, per_feature = FALSE)
      errors[i] <- mean(as.data.frame(recon)$Reconstruction.MSE)
      latent_feats <- h2o.deepfeatures(model, data_h2o, layer = 2)
      features_list[[i]] <- as.data.frame(latent_feats)
    }

    get_elbow_index <- function(x, y) {
      point1 <- c(x[1], y[1])
      point2 <- c(x[length(x)], y[length(y)])
      line_vec <- point2 - point1
      line_vec <- line_vec / sqrt(sum(line_vec^2))
      vec_from_first <- cbind(x - point1[1], y - point1[2])
      proj_lengths <- vec_from_first %*% line_vec
      proj_points <- matrix(rep(point1, each = length(x)), ncol = 2, byrow = TRUE) + proj_lengths %*% t(line_vec)
      dists <- sqrt(rowSums((vec_from_first - (proj_points - matrix(rep(point1, each = length(x)), ncol = 2, byrow = TRUE)))^2))
      which.max(dists)
    }

    elbow_idx <- get_elbow_index(latent_dims, errors)
    best_latent_dim <- latent_dims[elbow_idx]
    latent <- features_list[[elbow_idx]]

    compute_mmd <- function(X, Y, sigma = 1) {
      Kxx <- kernelMatrix(rbfdot(sigma = sigma), X)
      Kyy <- kernelMatrix(rbfdot(sigma = sigma), Y)
      Kxy <- kernelMatrix(rbfdot(sigma = sigma), X, Y)
      mean(Kxx) + mean(Kyy) - 2 * mean(Kxy)
    }

    errore_on_copula <- function(X1, X2) {
      n <- min(nrow(X1), nrow(X2))
      u1 <- apply(X1[1:n, ], 2, rank) / (n + 1)
      u2 <- apply(X2[1:n, ], 2, rank) / (n + 1)
      c1 <- pobs(u1)
      c2 <- pobs(u2)
      mean((c1 - c2)^2)
    }

    results <- data.frame()
    best_error <- Inf
    best_sample <- NULL
    best_k <- NULL
    max_k <- floor(0.2 * nrow(X))

    for (k in 10:max_k) {
      kmeans_result <- kmeans(latent, centers = k, nstart = 10)
      centers <- kmeans_result$centers
      idx <- sapply(1:k, function(i) {
        which.min(rowSums((latent - matrix(centers[i, ], nrow(latent), ncol(latent), byrow=TRUE))^2))
      })
      sample_scaled <- X[idx, , drop = FALSE]
      sample_df <- df_num[idx, ]

      pop_means <- colMeans(df_num, na.rm = TRUE)
      sample_means <- colMeans(sample_df, na.rm = TRUE)
      abs_error <- abs(pop_means - sample_means)
      rel_error <- abs_error / ifelse(pop_means == 0, NA, abs(pop_means)) * 100

      if (forza_chiave && !is.null(key_vars)) {
        chiave_bias <- rel_error[names(rel_error) %in% key_vars]
        if (any(chiave_bias > soglia, na.rm = TRUE)) next
      }

      mean_diff <- colMeans(X) - colMeans(sample_scaled)
      cov_inv <- ginv(cov(X))
      mahal <- sqrt(t(mean_diff) %*% cov_inv %*% mean_diff)
      mmd <- compute_mmd(X, sample_scaled, sigma = 1)
      errore <- errore_on_copula(X, sample_scaled)
      mmd_norm <- min(1, mmd / 0.2)
      mahal_norm <- min(1, mahal / 5)
      errore_norm <- min(1, errore / 0.1)
      combined_error <- 0.4 * mmd_norm + 0.4 * mahal_norm + 0.2 * errore_norm
      penalty_k <- log(k) / log(nrow(X))
      combined_error <- combined_error + 0.3 * penalty_k^2

      if (combined_error < best_error) {
        best_error <- combined_error
        best_k <- k
        best_sample <- sample_df
      }
    }

    if (is.null(best_sample)) {
      showNotification("❌ Nessun campione soddisfa i vincoli richiesti sulle variabili chiave.", type = "error")
      return(NULL)
    }

    best_k_global(best_k)
    pop_means <- colMeans(df_num, na.rm = TRUE)
    sample_means <- colMeans(best_sample, na.rm = TRUE)
    abs_error <- abs(pop_means - sample_means)
    rel_error <- abs_error / ifelse(pop_means == 0, NA, abs(pop_means)) * 100
    sample_se <- colSds(as.matrix(best_sample), na.rm = TRUE) / sqrt(nrow(best_sample))

    report <- data.frame(
      Variabile = names(pop_means),
      Pop_Mean = round(pop_means, 2),
      Sample_Mean = round(sample_means, 2),
      Abs_Error = round(abs_error, 2),
      Rel_Error_Percent = round(rel_error, 2),
      Sample_SE = round(sample_se, 2)
    )
    report <- report[order(-report$Rel_Error_Percent), ]
    report_iniziale(report)
  }

  observeEvent(input$run, {
    req(dati())
    estraiCampione(dati()$df, dati()$df_num, input$key_vars, forza_chiave = FALSE, soglia = input$soglia_bias)
  })

  observeEvent(input$run2, {
    req(dati(), best_k_global())
    estraiCampione(dati()$df, dati()$df_num, input$key_vars, forza_chiave = TRUE, soglia = input$soglia_bias)
  })

  output$biasPlot <- renderPlot({
    req(report_iniziale())
    ggplot(report_iniziale()[report_iniziale()$Rel_Error_Percent > 5, ], 
           aes(x = reorder(Variabile, -Rel_Error_Percent), y = Rel_Error_Percent)) +
      geom_bar(stat = "identity", fill = "salmon") +
      coord_flip() +
      labs(title = "Bias % sopra soglia per variabili", x = "Variabile", y = "% Bias") +
      geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
      theme_minimal()
  })

  output$biasTable <- renderTable({
    req(report_iniziale())
    head(report_iniziale(), 10)
  })

  output$kOutput <- renderPrint({
    req(report_iniziale(), best_k_global())
    rel_error <- report_iniziale()$Rel_Error_Percent
    cat("K ottimale:", best_k_global(),
        "
Bias%:", round(quantile(rel_error, 0.5, na.rm = TRUE), 2), "%")
  })
}

shinyApp(ui = ui, server = server)
