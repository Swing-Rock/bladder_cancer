# plots_live.R — the fast, interactive plots shown while you move the sliders.
# (The figures saved on export come from each script's own make_plots().)

APP_REASON_COLOURS <- c(pass = "grey75", low_counts = "#E69F00", low_features = "#D55E00",
                        high_counts = "#56B4E9", high_features = "#0072B2", mito = "#CC79A7",
                        doublet = "#009E73", contamination = "#7F3C8D", small_area = "#8C564B",
                        large_area = "#17BECF", high_control = "#CC79A7")

reason_colours <- function(env) {
  cols <- if (exists("REASON_COLOURS", envir = env, inherits = FALSE)) env$REASON_COLOURS else APP_REASON_COLOURS
  c(cols, APP_REASON_COLOURS[setdiff(names(APP_REASON_COLOURS), names(cols))])
}

fmt_num <- function(v) {
  ifelse(!is.finite(v), "none",
         ifelse(abs(v) >= 100, formatC(round(v), format = "d", big.mark = ","), formatC(v, digits = 3, format = "g")))
}

theme_qc <- function() theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
                                                     plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9.5))

# One histogram per slider metric, with the removed zones shaded and the cutoffs dashed.
# th = NULL: QC not applied yet, so just the distributions.
plot_distributions <- function(meta, specs, th, sources, note = NULL) {
  panels <- lapply(specs, function(s) {
    x <- meta[[s$col]]
    ok <- is.finite(x)
    if (!any(ok)) return(NULL)
    lo <- if (!is.null(s$lo) && !is.null(th)) th[[s$lo]] else -Inf
    hi <- if (!is.null(s$hi) && !is.null(th)) th[[s$hi]] else Inf
    removed <- sum(ok & (x < lo | x > hi))
    tf <- if (s$scale == "log") function(v) log10(pmax(v, 1)) else identity
    df <- data.frame(v = tf(x[ok]))
    g <- ggplot(df, aes(v)) + geom_histogram(bins = 120, fill = "#5B7DB1", colour = NA)
    if (is.finite(lo)) g <- g + annotate("rect", xmin = -Inf, xmax = tf(lo), ymin = -Inf, ymax = Inf, fill = "#D62728", alpha = 0.12) +
      geom_vline(xintercept = tf(lo), colour = "#D62728", linetype = "dashed", linewidth = 0.7)
    if (is.finite(hi)) g <- g + annotate("rect", xmin = tf(hi), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#D62728", alpha = 0.12) +
      geom_vline(xintercept = tf(hi), colour = "#D62728", linetype = "dashed", linewidth = 0.7)
    src <- function(nm) if (!is.null(nm) && isTRUE(unname(sources[nm]) == "manual")) " (manual)" else ""
    cut_txt <- paste(c(if (!is.null(s$lo)) paste0("min ", fmt_num(lo), src(s$lo)),
                       if (!is.null(s$hi)) paste0("max ", fmt_num(hi), src(s$hi))), collapse = " · ")
    sub <- if (is.null(th)) sprintf("median %s%s", fmt_num(median(x[ok])), note %||% "") else
      sprintf("%s  —  %s removed (%.1f%%)", cut_txt, format(removed, big.mark = ","), 100 * removed / sum(ok))
    g <- g + labs(title = s$label, subtitle = sub,
                  x = NULL, y = "count") + theme_qc()
    if (s$scale == "log") g <- g + scale_x_continuous(labels = function(b) fmt_num(10^b))
    g
  })
  panels <- Filter(Negate(is.null), panels)
  wrap_plots(panels, ncol = 2)
}

subsample <- function(meta, n_max, seed = 1) {
  if (nrow(meta) <= n_max) return(meta)
  set.seed(seed)
  meta[sort(sample(nrow(meta), n_max)), , drop = FALSE]
}

# counts vs genes coloured by QC reason, plus counts vs the "max" metric (mito / control)
plot_scatter <- function(meta, specs, th, colours, n_max) {
  d <- subsample(meta, n_max)
  d$qc_reason <- factor(d$qc_reason, levels = intersect(names(colours), c("pass", unique(d$qc_reason))))
  d <- d[order(d$qc_reason == "pass", decreasing = TRUE), ]   # removed cells on top
  cs <- specs[[1]]; fs <- specs[[2]]
  g1 <- ggplot(d, aes(.data[[cs$col]], .data[[fs$col]], colour = qc_reason)) +
    geom_point(size = 0.6, alpha = 0.7, shape = 16) + scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = colours, drop = TRUE, name = "QC result") +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = sprintf("%s vs %s", cs$label, tolower(fs$label)), x = cs$label, y = fs$label) + theme_qc()
  ms <- Filter(function(s) s$scale == "linear", specs)
  if (!length(ms)) return(g1)
  ms <- ms[[1]]
  g2 <- ggplot(d, aes(.data[[cs$col]], .data[[ms$col]], colour = qc_reason)) +
    geom_point(size = 0.6, alpha = 0.7, shape = 16) + scale_x_log10() +
    scale_colour_manual(values = colours, drop = TRUE, guide = "none") +
    labs(title = sprintf("%s vs %s", cs$label, ms$label), x = cs$label, y = ms$label) + theme_qc()
  if (!is.null(th) && isTRUE(is.finite(th[[ms$hi]]))) g2 <- g2 + geom_hline(yintercept = th[[ms$hi]], colour = "#D62728", linetype = "dashed")
  if (!is.null(th) && isTRUE(is.finite(th[[cs$lo]]))) g2 <- g2 + geom_vline(xintercept = max(th[[cs$lo]], 1), colour = "#D62728", linetype = "dashed")
  g1 + g2 + plot_layout(guides = "collect")
}

plot_spatial <- function(meta, colour_by, colours, image, n_max, point_size, removed_only = FALSE) {
  d <- meta[is.finite(meta$x) & is.finite(meta$y), , drop = FALSE]
  if (!nrow(d)) return(ggplot() + annotate("text", 0, 0, label = "No spatial coordinates in this sample") + theme_void())
  if (removed_only) d <- d[!d$keep, , drop = FALSE]
  d <- subsample(d, n_max)
  g <- ggplot(d, aes(x, -y))
  if (!is.null(image)) g <- g + annotation_raster(image$raster, xmin = 0, xmax = image$width, ymin = -image$height, ymax = 0)
  if (colour_by == "qc_reason") {
    d$qc_reason <- factor(d$qc_reason, levels = intersect(names(colours), c("pass", unique(d$qc_reason))))
    d <- d[order(d$qc_reason == "pass", decreasing = TRUE), ]
    g <- g + geom_point(data = d, aes(colour = qc_reason), size = point_size, shape = 16, alpha = 0.9) +
      scale_colour_manual(values = colours, drop = TRUE, name = "QC result") +
      guides(colour = guide_legend(override.aes = list(size = 3)))
  } else {
    v <- d[[colour_by]]
    lab <- colour_by
    if (grepl("^nCount|^nFeature|area$", colour_by)) { v <- log10(pmax(v, 1)); lab <- paste0("log10 ", colour_by) }
    d$.v <- v
    g <- g + geom_point(data = d, aes(colour = .v), size = point_size, shape = 16) +
      scale_colour_viridis_c(name = lab, option = if (grepl("mt|control", colour_by)) "magma" else "viridis")
  }
  g + coord_fixed() + theme_void(base_size = 12) +
    theme(legend.position = "right", plot.background = element_rect(fill = "white", colour = NA))
}

plot_extra <- function(meta, th) {
  panels <- list()
  d <- meta[!is.na(meta$scDblFinder.score), ]
  if (nrow(d)) {
    g <- ggplot(d, aes(scDblFinder.score, fill = scDblFinder.class)) + geom_histogram(bins = 60) +
      scale_fill_manual(values = c(singlet = "grey70", doublet = "#009E73"), name = NULL) +
      labs(title = "scDblFinder doublet scores", subtitle = sprintf("%s doublets (%.1f%% of scored cells); cutoff %s",
           format(sum(d$scDblFinder.class == "doublet"), big.mark = ","), 100 * mean(d$scDblFinder.class == "doublet"),
           fmt_num(th$doublet_threshold)), x = "doublet score", y = "cells") + theme_qc()
    if (is.finite(th$doublet_threshold)) g <- g + geom_vline(xintercept = th$doublet_threshold, linetype = "dashed", colour = "#D62728")
    panels <- c(panels, list(g))
  }
  d <- meta[!is.na(meta$decontX_contamination), ]
  if (nrow(d)) {
    g <- ggplot(d, aes(decontX_contamination)) + geom_histogram(bins = 60, fill = "#4C9BE8") +
      labs(title = "DecontX ambient contamination", subtitle = sprintf("median %.1f%%", 100 * median(d$decontX_contamination)),
           x = "fraction of counts from ambient RNA", y = "cells") + theme_qc()
    panels <- c(panels, list(g))
  }
  if (!length(panels)) return(NULL)
  wrap_plots(panels, ncol = 2)
}
