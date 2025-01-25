custom_plot_binary_TSP <-function (Data, classifier, ref = NULL, prediction = NULL, platform = NULL, 
          classes = NULL, platforms_ord = NULL, top_anno = c("ref", 
                                                             "prediction", "platform")[1], title = "", binary_col = c("white", 
                                                                                                                      "black", "gray"), ref_col = NULL, pred_col = NULL, platform_col = NULL, 
          show_ref = TRUE, show_predictions = TRUE, show_platform = TRUE, 
          show_scores = TRUE, show_rule_name = TRUE, legend = TRUE, 
          cluster_cols = TRUE, cluster_rows = TRUE, anno_height = 0.03, 
          score_height = 0.03, margin = c(0, 5, 0, 5)) 
{
  if (class(classifier)[1] != "OnevsrestScheme_TSP") {
    stop("classifier should be OnevsrestScheme_TSP object from train_one_vs_rest_TSP function!")
  }
  else {
    C <- classifier
  }
  if (!class(Data)[1] %in% c("multiclassPairs_object", "ExpressionSet", 
                             "data.frame", "matrix")) {
    stop("Data should be class:\n  matrix/data.frame/ExpressionSet/multiclassPairs_object from ReadData function!")
  }
  if (is.data.frame(Data)) {
    D <- Data
  }
  if (is.matrix(Data)) {
    D <- as.data.frame(Data, stringsAsFactors = FALSE)
  }
  if (class(Data)[1] == "ExpressionSet") {
    if (!requireNamespace("Biobase", quietly = TRUE)) {
      message("ExpressionSet is used and 'Biobase' package from Bioconductor is needed!")
      stop("Visit their website or install Biobase package using:\n      if (!requireNamespace('BiocManager', quietly = TRUE)) {\n      install.packages('BiocManager')\n      }\n      BiocManager::install('Biobase')", 
           call. = FALSE)
    }
    else {
      requireNamespace("Biobase")
    }
    D <- as.data.frame(Biobase::exprs(Data), stringsAsFactors = FALSE)
  }
  if (class(Data)[1] == "multiclassPairs_object") {
    D <- Data$data$Data
  }
  if (is.null(rownames(D))) {
    stop("Provide feature/gene names as rownames in the Data matrix!")
  }
  if (!is.null(classes)) {
    if (any(!classes %in% names(classifier[["classifiers"]]))) {
      message("These classes are not found in the classifier object:")
      message(paste0(classes[!classes %in% names(classifier[["classifiers"]])], 
                     collapse = " "))
      stop("classes names in classes argument should be similar to the names of the classifiers in classifier object!")
    }
  }
  else {
    classes <- names(classifier[["classifiers"]])
  }
  if (any(!names(classifier[["classifiers"]]) %in% classes)) {
    message("Classes argument misses the following classes, they will be removed from the heatmap:")
    message(paste0(names(classifier[["classifiers"]])[!names(classifier[["classifiers"]]) %in% 
                                                        classes], collapse = " "))
  }
  if (class(Data)[1] == "multiclassPairs_object") {
    if (is.null(ref)) {
      L <- Data$data$Labels
    }
  }
  if ((is.character(ref) | is.factor(ref)) & class(Data)[1] != 
      "ExpressionSet") {
    L <- as.character(ref)
  }
  if (class(Data)[1] == "ExpressionSet") {
    if (is.character(ref) & length(ref) == 1) {
      if (ref %in% Biobase::varLabels(Data)) {
        L <- as.character(Biobase::pData(Data)[, ref])
      }
      else {
        message(capture.output(cat("Phenotype data has these variables:", 
                                   Biobase::varLabels(Data), fill = TRUE)))
        stop("Ref label variable is not found in the phenotype data of your ExpressionSet")
      }
    }
    if ((is.character(ref) | is.factor(ref)) & length(ref) != 
        1) {
      L <- as.character(ref)
    }
  }
  if (is.null(ref) & class(Data)[1] != "multiclassPairs_object") {
    L <- NULL
  }
  if (!is.null(ref)) {
    if (length(L) != ncol(D)) {
      message("Number of samples: ", ncol(D))
      message("Labels length: ", length(L))
      stop("Labels vector length are not equal to\n       samples in data")
    }
  }
  if (class(Data)[1] == "multiclassPairs_object") {
    if (is.null(platform)) {
      P <- Data$data$Platform
    }
  }
  if ((is.character(platform) | is.factor(platform)) & class(Data)[1] != 
      "ExpressionSet") {
    P <- as.character(platform)
  }
  if (class(Data)[1] == "ExpressionSet") {
    if (is.character(platform) & length(platform) == 1) {
      if (platform %in% Biobase::varLabels(Data)) {
        P <- as.character(Biobase::pData(Data)[, platform])
      }
      else {
        message(capture.output(cat("Phenotype data has these variables:", 
                                   Biobase::varLabels(Data), fill = TRUE)))
        stop("Platform/study label variable is not found in the phenotype data of your ExpressionSet")
      }
    }
    if ((is.character(platform) | is.factor(platform)) & 
        length(platform) != 1) {
      P <- as.character(platform)
    }
  }
  if (is.null(platform) & class(Data)[1] != "multiclassPairs_object") {
    P <- NULL
  }
  if (!is.null(platform)) {
    if (length(P) != ncol(D)) {
      message("Number of samples: ", ncol(D))
      message("Labels length: ", length(P))
      stop("Platform labels vector length are not equal to\n       samples in data")
    }
  }
  if (!is.null(platforms_ord) & !is.null(P)) {
    if (any(!platforms_ord %in% P)) {
      message("These platform/study in platforms_ord are not in found the platform labels:")
      message(platforms_ord[!platforms_ord %in% P])
      stop("platforms_ord argument should have similar names of the platforms/studies in the data!")
    }
  }
  else {
    platforms_ord <- unique(P)
  }
  if (!is.null(prediction) & any(class(prediction) %in% "OneVsRestTSP prediction")) {
    pred <- prediction
  }
  else {
    pred <- NULL
  }
  if (!is.null(prediction)) {
    if (nrow(pred) != ncol(D)) {
      message("Number of samples in the data: ", ncol(D))
      message("Number of samples in the prediction: ", 
              ncol(pred))
      stop("prediction should be for the same data!\n     Use predict_one_vs_rest_TSP to generate it for this data!")
    }
  }
  if (is.null(pred) & is.null(L) & is.null(P)) {
    stop("No available ref, prediction, or platform labels!\n     One of them atleast is needed.")
  }
  if (top_anno == "ref" & is.null(L)) {
    stop("top annotation (top_anno) is ref while there is no ref labels available!")
  }
  if (top_anno == "ref" & !show_ref) {
    message("show_ref was turned to TRUE because top_anno is 'ref'!")
    show_ref <- TRUE
  }
  if (top_anno == "prediction" & is.null(pred)) {
    stop("top annotation (top_anno) is prediction while there is no prediction dataframe available! Use predict_one_vs_rest_TSP function to generate it!")
  }
  if (top_anno == "prediction" & !show_predictions) {
    message("show_predictions was turned to TRUE because top_anno is 'prediction'!")
    show_predictions <- TRUE
  }
  if (top_anno == "platform" & is.null(P)) {
    stop("top annotation (top_anno) is platform while there is no platform labels available!")
  }
  if (top_anno == "platform" & !show_platform) {
    message("show_platform was turned to TRUE because top_anno is 'platform'!")
    show_platform <- TRUE
  }
  if (any(!top_anno %in% c("ref", "prediction", "platform")) | 
      !is.character(top_anno) | length(top_anno) != 1) {
    stop("Top annotation argument should be character with one of these options:\n     ref prediction platform")
  }
  if (!is.character(title) | length(title) != 1) {
    stop("Title argument should be character input!")
  }
  if (!is.character(binary_col) | length(binary_col) != 3) {
    stop("binary_col should be character input with length of 3!\n     Three colors are needed for rules with false, true and NAs.\n     By default it is c('white','black','gray')")
  }
  if (show_ref & !is.null(L) & !is.null(ref_col)) {
    if (!is.character(ref_col) | any(!classes %in% names(ref_col))) {
      stop("ref_col should be named character vector for all classes!")
    }
  }
  if (show_predictions & !is.null(pred) & !is.null(pred_col)) {
    if (!is.character(pred_col) | any(!classes %in% names(pred_col))) {
      stop("pred_col should be named character vector for all classes!")
    }
  }
  if (show_platform & !is.null(P) & !is.null(platform_col)) {
    if (!is.character(platform_col) | any(!P %in% names(platform_col))) {
      stop("platform_col should be named character vector for all platforms/studies!")
    }
  }
  xx_colors <- c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", 
                 "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", 
                 "#fabed4", "#469990", "#dcbeff", "#9A6324", "#fffac8", 
                 "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", 
                 "#a9a9a9")
  xx_colors2 <- c("#aaffc3", "#808000", "#ffd8b1", "#000075", 
                  "#a9a9a9", "#469990", "#dcbeff", "#9A6324", "#fffac8", 
                  "#800000", "#e6194B", "#3cb44b", "#ffe119", "#4363d8", 
                  "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", 
                  "#fabed4")
  if (is.null(ref_col) & !is.null(L)) {
    if (length(classes) < 20) {
      ref_col <- xx_colors[1:length(classes)]
    }
    else {
      ref_col <- sample(xx_colors, size = length(classes), 
                        replace = T)
    }
    names(ref_col) <- classes
  }
  if (is.null(pred_col) & !is.null(pred)) {
    if (length(classes) < 20) {
      pred_col <- xx_colors[1:length(classes)]
    }
    else {
      pred_col <- sample(xx_colors, size = length(classes), 
                         replace = T)
    }
    names(pred_col) <- classes
  }
  if (is.null(platform_col) & !is.null(P)) {
    if (length(platforms_ord) < 20) {
      platform_col <- xx_colors2[1:length(platforms_ord)]
    }
    else {
      platform_col <- sample(xx_colors, size = length(platforms_ord), 
                             replace = T)
    }
    names(platform_col) <- platforms_ord
  }
  sam_names <- colnames(D)
  if (top_anno == "ref") {
    lab <- L
    groups <- classes
    groups_col <- ref_col
  }
  if (top_anno == "prediction") {
    lab <- pred$max_score
    groups <- classes
    groups_col <- pred_col
  }
  if (top_anno == "platform") {
    lab <- P
    groups <- platforms_ord
    groups_col <- platform_col
  }
  anno_ord <- c("ref", "prediction", "platform")
  anno_ord <- anno_ord[c(!is.null(L) & show_ref, !is.null(pred) & 
                           show_predictions, !is.null(P) & show_platform)]
  anno_ord <- anno_ord[!anno_ord %in% top_anno]
  
  ## Check if the labels are the same as the prediction labels
  ### If so do as it before
  uniq_pred <- unique(lab)
  pred_eq_ref <- all(unlist(groups) %in% unlist(uniq_pred))
  if (pred_eq_ref) {
    if (cluster_cols & (top_anno %in% c("ref", "prediction"))) {
      tmp <- c()
      for (i in groups) {
        select_samples <- sam_names[lab == i]
        tmp_r <- C$classifiers[[i]]$TSPs
        tmp_binary <- D[tmp_r[, 1], select_samples, drop = FALSE] > 
          D[tmp_r[, 2], select_samples, drop = FALSE]
        print(i)
        d <- dist(t(tmp_binary[, select_samples, drop = FALSE]), 
                  method = "euclidean")
        fit <- hclust(d, method = "ward.D2")
        tmp <- c(tmp, fit$labels[fit$order])
      }
      sam_ord <- order(match(sam_names, tmp))[1:length(tmp)]
      rm(tmp)
    }
    if (cluster_cols & top_anno == "platform") {
      tmp_r <- data.frame(matrix(NA, ncol = 2, nrow = 0), 
                          stringsAsFactors = FALSE)
      for (cl in classes) {
        tmp_r <- rbind(tmp_r, C$classifiers[[cl]]$TSPs)
      }
      tmp <- c()
      for (i in groups) {
        select_samples <- sam_names[lab == i]
        tmp_binary <- D[tmp_r[, 1], select_samples, drop = FALSE] > 
          D[tmp_r[, 2], select_samples, drop = FALSE]
        d <- dist(t(tmp_binary[, select_samples, drop = FALSE]), 
                  method = "euclidean")
        fit <- hclust(d, method = "ward.D2")
        tmp <- c(tmp, fit$labels[fit$order])
      }
      sam_ord <- order(match(sam_names, tmp))[1:length(tmp)]
      rm(tmp)
    }
    if (!cluster_cols) {
      sam_ord <- order(match(lab, groups))
    }
  } else {
    sam_ord <- order(match(lab, groups))
  }
  
  D <- D[, sam_ord]
  L <- L[sam_ord]
  P <- P[sam_ord]
  pred <- pred[sam_ord, ]
  lab <- lab[sam_ord]
  sam_names <- sam_names[sam_ord]
  num_sam <- ncol(D)
  splits <- table(lab)[order(match(names(table(lab)), groups))]
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  {
    AreaStart <- 0.94
    SizeUnit <- anno_height
    Size <- SizeUnit * 1
    AreaEnd <- AreaStart - Size
    par(fig = c(0, 1, AreaEnd, AreaStart), mar = margin, 
        mgp = c(3, 0.5, 0), new = FALSE)
    plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", yaxs = "i", 
         xlab = "", ylab = "", main = "", xlim = c(0, num_sam), 
         ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n")
    text_positions <- cumsum(splits)[1]/2
    for (i in 1:(length(cumsum(splits)) - 1)) {
      text_positions <- c(text_positions, ((cumsum(splits)[i + 
                                                             1] - cumsum(splits)[i])/2 + cumsum(splits)[i]))
    }
    mtext(groups, side = 3, line = 0, outer = FALSE, at = text_positions, 
          adj = NA, padj = NA, cex = 0.8, col = groups_col, 
          font = NA)
    mtext(title, side = 3, line = -1, outer = TRUE, font = 2)
    axis(side = 2, at = 0.5, labels = c(ref = "Ref. labels", 
                                        prediction = "Predictions", platform = "Platform/Study")[top_anno], 
         las = 1, cex.axis = 0.7, tick = 0)
    for (f in groups) {
      for (g in which(lab == f)) {
        rect(g - 1, 0, g, 1, col = groups_col[f], border = NA)
      }
    }
    box(lwd = 1)
    li <- cumsum(splits)
    abline(v = li, lwd = 1.5, lty = 1, col = "black")
  }
  for (i in anno_ord) {
    {
      Gap <- 0
      AreaStart <- AreaStart - Size - Gap
      SizeUnit <- anno_height
      Size <- SizeUnit * 1
      AreaEnd <- AreaStart - Size
      par(fig = c(0, 1, AreaEnd, AreaStart), mar = margin, 
          mgp = c(3, 0.5, 0), new = TRUE)
      plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", yaxs = "i", 
           xlab = "", ylab = "", main = "", xlim = c(0, 
                                                     num_sam), ylim = c(0, 1), xaxt = "n", yaxt = "n", 
           bty = "n")
      axis(side = 2, at = 0.5, labels = c(ref = "Ref. labels", 
                                          prediction = "Predictions", platform = "Platform/Study")[i], 
           las = 1, cex.axis = 0.7, tick = 0)
      if (i == "ref") {
        tmp_color <- ref_col
        tmp_lab <- L
      }
      if (i == "prediction") {
        tmp_color <- pred_col
        tmp_lab <- pred$max_score
      }
      if (i == "platform") {
        tmp_color <- platform_col
        tmp_lab <- P
      }
      for (f in unique(tmp_lab)) {
        for (g in which(tmp_lab == f)) {
          rect(g - 1, 0, g, 1, col = tmp_color[f], border = NA)
        }
      }
      box(lwd = 1)
      li <- cumsum(splits)
      abline(v = li, lwd = 1.5, lty = 1, col = "black")
    }
  }
  if (top_anno == "platform") {
    groups <- classes
  }
  if (show_scores & !is.null(pred)) {
    score_matrix <- pred[, groups, drop = FALSE]
    Gap <- 0.01
    AreaStart <- AreaStart - Size - Gap
    SizeUnit <- score_height
    Size <- SizeUnit * 1
    AreaEnd <- AreaStart - Size
    for (class in colnames(score_matrix)) {
      par(fig = c(0, 1, AreaEnd, AreaStart), mar = margin, 
          mgp = c(3, 0.5, 0), new = TRUE)
      barplot(as.numeric(score_matrix[, class]), col = pred_col[class], 
              space = F, xaxs = "i", yaxs = "i", xlim = c(0, 
                                                          num_sam), border = NA, ylim = c(0, 1), xaxt = "n", 
              yaxt = "n", bty = "n", ylab = "", xlab = "")
      box(lwd = 1)
      axis(2, at = 0.5, labels = paste("Scores:", class), 
           las = 1, cex.axis = 0.7, tick = 0)
      axis(4, at = c(0.1, 0.5, 0.9), labels = c("0", "0.5", 
                                                "1"), las = 1, cex.axis = 0.4, tick = FALSE)
      li <- cumsum(splits)
      abline(v = li, lwd = 1.5, lty = 1, col = "black")
      if (class == colnames(score_matrix)[ncol(score_matrix)]) {
        next
      }
      Gap <- 0
      AreaStart <- AreaStart - Size - Gap
      SizeUnit <- score_height
      Size <- SizeUnit * 1
      AreaEnd <- AreaStart - Size
    }
  }
  Size <- (AreaEnd - (0.005 * length(groups)) - 0.08)/length(groups)
  for (o in groups) {
    tmp <- C$classifiers[[o]]$TSPs
    binary <- D[tmp[, 1], , drop = FALSE] > D[tmp[, 2], 
                                              , drop = FALSE]
    binary <- binary + 1
    #binary<- na.omit(binary)
    binary[is.na(binary)] <- 0
    if (cluster_rows & pred_eq_ref) {
      d <- dist(binary, method = "euclidean")
      fit <- hclust(d, method = "ward.D2")
      tmp_ord <- fit$labels[fit$order]
    }
    else {
      tmp_ord <- 1:nrow(binary)
    }
    binary <- binary[tmp_ord, ]
    Gap <- 0.005
    AreaStart <- AreaEnd - Gap
    AreaEnd <- AreaStart - Size
    par(fig = c(0, 1, AreaEnd, AreaStart), mar = margin, 
        mgp = c(3, 0.5, 0), new = TRUE)
    myplot <- plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", 
                   yaxs = "i", xlab = "", ylab = "", main = "", xlim = c(0, 
                                                                         num_sam), ylim = c(0, length(tmp_ord)), xaxt = "n", 
                   yaxt = "n", bty = "n")
    if (show_rule_name) {
      rule_names <- paste(tmp[, 1], tmp[, 2], sep = ">")
      axis(4, at = seq(0, (length(tmp_ord) - 1), 1) + 
             0.5, labels = rule_names, las = 1, cex.axis = 0.5, 
           tick = 0)
    }
    for (f in 1:ncol(binary)) {
      for (g in 1:nrow(binary)) {
        rect(f - 1, g, f, g - 1, col = binary_col[binary[g, 
                                                         f]], border = NA, lwd = 0)
      }
    }
    axis(2, at = nrow(binary)/2, labels = paste(o, "\n", 
                                                nrow(binary), "rules"), las = 1, cex.axis = 0.7, 
         tick = 0, col = groups_col[o])
    box(lwd = 1)
    li <- cumsum(splits)
    abline(v = li[-length(li)], lwd = 3, lty = 3, col = "red")
  }
  if (legend) {
    par(fig = c(0, 1, 0.02, (AreaEnd - 0.01)), mar = margin, 
        mgp = c(3, 0.5, 0), new = TRUE)
    plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", yaxs = "i", 
         xlab = "", ylab = "", main = "", xlim = c(0, num_sam), 
         ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n")
    if (!is.null(P) & show_platform) {
      legend(x = "topright", title = "Platform", ncol = length(platforms_ord), 
             cex = 0.5, legend = platforms_ord, fill = platform_col)
    }
    if (!is.null(L) & show_ref) {
      title <- "Ref"
      if (!is.null(pred) & show_predictions & length(ref_col) == 
          length(pred_col) & all(ref_col %in% pred_col)) {
        title <- "Classes"
      }
      legend(x = "topleft", title = title, ncol = length(classes), 
             cex = 0.5, legend = names(ref_col), fill = ref_col)
    }
    if (!is.null(pred) & show_predictions & (length(ref_col) != 
                                             length(pred_col) | any(!ref_col %in% pred_col))) {
      legend(x = "top", title = "Predictions", ncol = length(unique(pred$max_score)), 
             cex = 0.5, legend = names(pred_col), fill = pred_col)
    }
  }
}
