
SingleDimPlot_rast <- function (data, dims, col.by = NULL, cols = NULL, pt.size = NULL,
          shape.by = NULL, alpha.by = NULL, order = NULL, label = FALSE,
          repel = FALSE, label.size = 4, cells.highlight = NULL, cols.highlight = "#DE2D26",
          sizes.highlight = 1, na.value = "grey50")
{
  pt.size <- pt.size %||% Seurat:::AutoPointSize(data = data)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  }
  else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(cells.highlight = cells.highlight,
                                   cells.all = rownames(x = data), sizes.highlight = sizes.highlight %||%
                                     pt.size, cols.highlight = cols.highlight, col.base = cols[1] %||%
                                     "#C3C3C3", pt.size = pt.size)
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- "highlight"
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(data[, col.by]), ]
      }
    }
    else {
      order <- rev(x = c(order, setdiff(x = unique(x = data[,
                                                            col.by]), y = order)))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  }
  else {
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = "^\\d", x = col.by)) {
      col.by <- paste0("x", col.by)
    }
    else if (grepl(pattern = "-", x = col.by)) {
      col.by <- gsub(pattern = "-", replacement = ".",
                     x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
    warning("Cannot find alpha variable ", alpha.by,
            " in data, setting to NULL", call. = FALSE,
            immediate. = TRUE)
    alpha.by <- NULL
  }
  plot <- ggplot(data = data) + ggrastr::geom_point_rast(mapping = aes_string(x = dims[1],
                                                                y = dims[2], color = paste0("`", col.by, "`"),
                                                                shape = shape.by, alpha = alpha.by), size = pt.size) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(color = NULL)
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(plot = plot, id = col.by, repel = repel,
                          size = label.size)
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) ||
                                  cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    }
    else if (length(x = cols) == 1 && (cols %in% c("alphabet",
                                                   "alphabet2", "glasbey", "polychrome",
                                                   "stepped"))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])),
                                palette = cols)
      scale <- scale_color_manual(values = colors, na.value = na.value)
    }
    else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + theme_cowplot()
  return(plot)
}

FeaturePlot_rast <- function(object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
  c("lightgrey", "#ff0000", "#00ff00")} else {c("lightgrey", "blue")}, pt.size = NULL, order = FALSE, min.cutoff = NA, max.cutoff = NA,
reduction = NULL, split.by = NULL, shape.by = NULL, slot = "data",
blend = FALSE, blend.threshold = 0.5, label = FALSE, label.size = 4,
repel = FALSE, ncol = NULL, coord.fixed = FALSE, by.col = TRUE,
sort.cell = NULL, interactive = FALSE, combine = TRUE, rasterize=FALSE, xlab=NULL, xlab.size=8){
  if (!is.null(x = sort.cell)) {
    warning("The sort.cell parameter is being deprecated. Please use the order ",
            "parameter instead for equivalent functionality.",
            call. = FALSE, immediate. = TRUE)
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  if (interactive) {
    return(IFeaturePlot(object = object, feature = features[1],
                        dims = dims, reduction = reduction, slot = slot))
  }
  no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(),
                    axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold",
                                                                                           size = 14, margin = margin(r = 7)))
  reduction <- reduction %||% Seurat:::DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)),
                   `0` = {
                     warning("No colors provided, using default colors",
                             call. = FALSE, immediate. = TRUE)
                     default.colors
                   }, `1` = {
                     warning("Only one color provided, assuming specified is double-negative and augmenting with default colors",
                             call. = FALSE, immediate. = TRUE)
                     c(cols, default.colors[2:3])
                   }, `2` = {
                     warning("Only two colors provided, assuming specified are for features and agumenting with '",
                             default.colors[1], "' for double-negatives",
                             call. = FALSE, immediate. = TRUE)
                     c(default.colors[1], cols)
                   }, `3` = cols, {
                     warning("More than three colors provided, using only first three",
                             call. = FALSE, immediate. = TRUE)
                     cols[1:3]
                   })
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(dims, "ident",
                                              features), cells = cells, slot = slot)
  if (ncol(x = data) < 4) {
    stop("None of the requested features were found: ",
         paste(features, collapse = ", "), " in slot ",
         slot, call. = FALSE)
  }
  else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found",
         call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[,
                                                           feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[,
                                                           feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff,
                                              max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols,
  ]$maxcolors, no = length(x = cols))
  data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data),
                                     FUN = function(index) {
                                       data.feature <- as.vector(x = data[, index])
                                       min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index -
                                                                                    3], data.feature)
                                       max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index -
                                                                                    3], data.feature)
                                       data.feature[data.feature < min.use] <- min.use
                                       data.feature[data.feature > max.use] <- max.use
                                       if (brewer.gran == 2) {
                                         return(data.feature)
                                       }
                                       data.cut <- if (all(data.feature == 0)) {
                                         0
                                       }
                                       else {
                                         as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature),
                                                                          breaks = brewer.gran)))
                                       }
                                       return(data.cut)
                                     })
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  data$split <- if (is.null(x = split.by)) {
    Seurat:::RandomName()
  }
  else {
    switch(EXPR = split.by, ident = Idents(object = object)[cells,
                                                            drop = TRUE], object[[split.by, drop = TRUE]][cells,
                                                                                                          drop = TRUE])
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list", length = ifelse(test = blend,
                                                 yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[,
                                                                   dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[,
                                                               dims[2]])))
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold,
                                negative.color = cols[1])
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1], color.matrix[1, ],
                   as.vector(x = color.matrix))
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident,
                      , drop = FALSE]
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[,
                                                       features]) == 0]
      if (length(x = no.expression) != 0) {
        stop("The following features have no value: ",
             paste(no.expression, collapse = ", "),
             call. = FALSE)
      }
      data.plot <- cbind(data.plot[, c(dims, "ident")],
                         BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[,
                                                              feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      }
      else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, "ident",
                                   feature, shape.by)]

      if (rasterize==TRUE){

        plot <- SingleDimPlot_rast(data = data.single, dims = dims,
                              col.by = feature, order = order, pt.size = pt.size,
                              cols = cols.use, shape.by = shape.by, label = FALSE) +
          scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) +
          theme_cowplot() + theme(plot.title = element_text(hjust = 0.5),  plot.background=element_rect(colour="black"))
        if (!is.null(xlab)){
          plot <- plot +xlab(xlab)+theme(axis.title.x = element_text(size=xlab.size))}
      } else {
      plot <- Seurat:::SingleDimPlot(data = data.single, dims = dims,
                            col.by = feature, order = order, pt.size = pt.size,
                            cols = cols.use, shape.by = shape.by, label = FALSE) +
        scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) +
        theme_cowplot() + theme(plot.title = element_text(hjust = 0.5),  plot.background=element_rect(colour="black"))
      }
      if (label) {
        plot <- LabelClusters(plot = plot, id = "ident",
                              repel = repel, size = label.size)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA,
                                                         colour = "black"))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        }
        else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident),
                                                                    limits = ylims) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(axis.line.y = element_blank(),
                               axis.ticks.y = element_blank(), axis.text.y = element_blank(),
                               axis.title.y.left = element_blank())
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(axis.line.x = element_blank(),
                               axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                               axis.title.x = element_blank())
        }
      }
      else {
        plot <- plot + labs(title = feature)
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        }
        else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning("All cells have the same value (",
                    unique.feature.exp, ") of ", feature,
                    ".")
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            }
            else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad,
                                                                       guide = "colorbar"))
        }
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(x = plots,
                                              values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) >
                                                                                                                                  1, yes = levels(x = data$split)[ii], no = "")),
                                                                                              expand = c(0, 0)) + labs(x = features[1], y = features[2],
                                                                                                                       title = if (ii == 1) {
                                                                                                                         paste("Color threshold:", blend.threshold)
                                                                                                                       } else {
                                                                                                                         NULL
                                                                                                                       }) + no.right), after = 4 * ii - 1))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol,
                 no = length(x = features))
  legend <- if (blend) {
    "none"
  }
  else {
    split.by %iff% "none"
  }
  if (combine) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() +
                                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = ""),
                                                                   limits = ylims) + no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) +
                 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] +
                                         scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]),
                                                            limits = ylims) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots)%%length(x = features) ==
                      1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      }
      else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots),
                                                          f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    }
    else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff%
                            length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
  }
  return(plots)
}

`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}

DimPlot_rast <- function(object, dims = c(1, 2), cells = NULL, cols = NULL,
          pt.size = NULL, reduction = NULL, group.by = NULL, split.by = NULL,
          shape.by = NULL, order = NULL, shuffle = FALSE, seed = 1,
          label = FALSE, label.size = 4, label.color = "black",
          label.box = FALSE, repel = FALSE, cells.highlight = NULL,
          cols.highlight = "#DE2D26", sizes.highlight = 1, na.value = "grey50",
          ncol = NULL, combine = TRUE, rasterize=FALSE)
{
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  reduction <- reduction %||% Seurat:::DefaultDimReduc(object = object)
  cells <- cells %||% colnames(x = object)
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    cells <- sample(x = cells)
  }
  data <- Embeddings(object = object[[reduction]])[cells, dims]
  data <- as.data.frame(x = data)
  dims <- paste0(Key(object = object[[reduction]]), dims)
  object[["ident"]] <- Idents(object = object)
  orig.groups <- group.by
  group.by <- group.by %||% "ident"
  data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    data[, split.by] <- object[[split.by, drop = TRUE]]
  }
  plots <- lapply(X = group.by, FUN = function(x) {
    if (rasterize==T){
      plot <- SingleDimPlot_rast(data = data[, c(dims, x, split.by,
                                            shape.by)], dims = dims, col.by = x, cols = cols,
                            pt.size = pt.size, shape.by = shape.by, order = order,
                            label = FALSE, cells.highlight = cells.highlight,
                            cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,
                            na.value = na.value)
    } else {
    plot <- Seurat:::SingleDimPlot(data = data[, c(dims, x, split.by,
                                          shape.by)], dims = dims, col.by = x, cols = cols,
                          pt.size = pt.size, shape.by = shape.by, order = order,
                          label = FALSE, cells.highlight = cells.highlight,
                          cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,
                          na.value = na.value)}
    if (label) {
      plot <- LabelClusters(plot = plot, id = x, repel = repel,
                            size = label.size, split.by = split.by, box = label.box,
                            color = label.color)
    }
    if (!is.null(x = split.by)) {
      plot <- plot + FacetTheme() + facet_wrap(facets = vars(!!sym(x = split.by)),
                                               ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                                                 length(x = unique(x = data[, split.by]))
                                               }
                                               else {
                                                 ncol
                                               })
    }
    return(plot)
  })
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}

CurvePlot_rast <- function(object, sds = NULL, group.by = NULL, reduction = "umap",
          dims = 1:2, cols = NULL, label = T){
  object[["ident"]] <- Idents(object = object)
  group.by <- group.by %||% "ident"
  dims <- paste0(Key(object = object[[reduction]]), dims)
  curved <- bind_rows(lapply(names(slingCurves(sds)), function(x) {
    c <- slingCurves(sds)[[x]]
    d <- as.data.frame(c$s[c$ord, dims])
    d$curve <- x
    return(d)
  }))
  DimPlot_rast(object, cols = cols, label = label, group.by = group.by,
          reduction = reduction, rasterize=TRUE) + geom_path(aes_string(dims[1],
                                                        dims[2], linetype = "curve"), curved, size = 1)
}

`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

