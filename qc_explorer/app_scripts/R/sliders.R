# sliders.R — slider ranges and the rule that turns slider positions into cutoffs.
# Pure functions (no Shiny), so the same rule is used for the sample on screen and
# for samples saved in bulk.

# Slider settings for each metric, starting at the MAD cutoffs in `th`.
# UMIs/genes/area sliders work on a log10 scale; the app shows real values on them.
slider_config <- function(specs, meta, th, load_n = 0) {
  cfg <- lapply(specs, function(s) {
    s <- s[setdiff(names(s), c("min", "max", "step", "value"))]       # recomputed below (specs may be an old config)
    s$input_id <- s$input_id %||% sprintf("cut_%s_%d", s$id, load_n)   # fresh ID per load: no stale positions
    x <- meta[[s$col]]; x <- x[is.finite(x)]
    if (!length(x)) return(NULL)
    tf <- if (s$scale == "log") function(v) log10(pmax(v, 1)) else identity
    cuts <- c(if (!is.null(s$lo)) th[[s$lo]], if (!is.null(s$hi)) th[[s$hi]])
    cuts <- cuts[is.finite(cuts)]
    if (s$scale == "log") {
      mx <- ceiling(max(tf(c(x, cuts))) * 20) / 20 + 0.1; step <- 0.01; mn <- 0
    } else {
      mx <- max(ceiling(max(c(x, cuts))) + 1, 1); if (grepl("mt|mito", s$col)) mx <- min(mx, 100); step <- 0.1; mn <- 0
    }
    val <- c(if (!is.null(s$lo)) min(max(tf(th[[s$lo]]), mn), mx),
             if (!is.null(s$hi)) { h <- th[[s$hi]]; if (is.finite(h)) min(max(tf(h), mn), mx) else mx })
    c(s, list(min = mn, max = mx, step = step, value = round(val / step) * step))
  })
  cfg <- Filter(Negate(is.null), cfg)
  setNames(cfg, vapply(cfg, `[[`, "", "id"))
}

slider_sides <- function(c) if (!is.null(c$lo)) c(lo = 1, hi = 2) else c(hi = 1)

# A cutoff value -> slider position (log10 for log-scale sliders), kept inside the slider's range
slider_pos <- function(c, value) {
  tf <- if (c$scale == "log") log10(max(value, 1)) else value
  min(max(tf, c$min), c$max)
}

# Final cutoffs. A slider still where the MAD rule put it (`set`) uses the exact MAD value
# ("auto"); a slider that was moved uses its own value ("manual"). `exact` holds values typed
# into the number boxes (by cutoff name): the slider can only sit on its step grid, so while it
# is still where the typed value put it, the typed value is used as is.
resolve_cutoffs <- function(th, cfg, vals, set, exact = NULL) {
  used <- th; sources <- character(0)
  for (c in cfg) {
    v <- vals[[c$id]]; s0 <- set[[c$id]]
    sides <- slider_sides(c)
    for (side in names(sides)) {
      nm <- c[[side]]; i <- sides[[side]]
      ex <- exact[[nm]]
      if (!is.null(ex) && !is.null(v) && length(v) >= i && abs(slider_pos(c, ex) - v[i]) <= c$step / 2 + 1e-9) {
        if (isTRUE(all.equal(ex, th[[nm]]))) { sources[nm] <- "auto"; next }
        sources[nm] <- "manual"; used[[nm]] <- ex; next
      }
      if (is.null(v) || length(v) < i || is.null(s0) || abs(v[i] - s0[i]) <= c$step / 2 + 1e-9) { sources[nm] <- "auto"; next }
      sources[nm] <- "manual"
      val <- if (c$scale == "log") round(10^v[i]) else v[i]
      if (side == "hi" && v[i] >= c$max - c$step / 2) val <- Inf
      if (side == "lo" && v[i] <= c$min + c$step / 2) val <- 0
      used[[nm]] <- val
    }
  }
  list(th = used, auto = th, sources = sources)
}

# New slider settings for a sample whose MAD cutoffs may have changed since it was last
# on screen: sliders that were on "auto" move to the new MAD values, moved ones stay put.
restore_sliders <- function(specs, meta, th, load_n, old_cfg, old_vals, old_set) {
  cfg <- slider_config(specs, meta, th, load_n)
  set <- lapply(cfg, `[[`, "value")
  if (!is.null(old_vals)) for (id in names(cfg)) {
    c <- cfg[[id]]; v <- old_vals[[id]]; s0 <- old_set[[id]]; oc <- old_cfg[[id]]
    if (is.null(v) || is.null(s0) || is.null(oc)) next
    for (i in slider_sides(c)) if (length(v) >= i && abs(v[i] - s0[i]) > c$step / 2 + 1e-9)
      cfg[[id]]$value[i] <- min(max(v[i], c$min), c$max)
  }
  list(cfg = cfg, set = set)
}

# Cutoffs a sample-sheet row fixes (its cutoff_<name> columns): move the sliders there, as if typed
# into the number boxes (so they count as set by hand and keep their exact values).
# cut: named list of values by cutoff name. Returns list(cfg with the new positions, exact).
sheet_cutoff_sliders <- function(cfg, cut) {
  exact <- list()
  for (id in names(cfg)) {
    c <- cfg[[id]]; sides <- slider_sides(c)
    for (side in names(sides)) {
      nm <- c[[side]]; v <- cut[[nm]]
      if (is.null(v) || length(v) != 1 || is.na(v)) next
      cfg[[id]]$value[sides[[side]]] <- slider_pos(c, v)
      exact[[nm]] <- v
    }
  }
  list(cfg = cfg, exact = exact)
}
