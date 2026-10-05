# Output aggregation, publication tables and figures.

read_csv_shards <- function(path, pattern) {
  ff <- list.files(path, pattern = pattern, full.names = TRUE)
  if (!length(ff)) stop("No files found for pattern: ", pattern)
  do.call(rbind, lapply(ff, function(f) read.csv(f, stringsAsFactors = FALSE)))
}

dedupe_rows <- function(x, keys) {
  key <- do.call(paste, c(x[keys], sep = "|"))
  x[!duplicated(key, fromLast = TRUE), , drop = FALSE]
}

mean_sd_string <- function(x, digits = 3) {
  x <- x[is.finite(x)]
  if (!length(x)) return("--")
  if (length(x) == 1L) return(sprintf(paste0("%.", digits, "f"), x))
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(x), sd(x))
}

mean_se_values <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(Mean = NA_real_, SD = NA_real_, SE = NA_real_, N = 0))
  c(Mean = mean(x), SD = if (length(x) > 1L) sd(x) else NA_real_,
    SE = if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_, N = length(x))
}

summary_long <- function(data, group_cols, metrics, success_col = NULL) {
  if (!is.null(success_col)) data <- data[data[[success_col]] %in% TRUE, , drop = FALSE]
  keys <- unique(data[group_cols])
  ans <- list(); z <- 1L
  for (i in seq_len(nrow(keys))) {
    keep <- rep(TRUE, nrow(data))
    for (g in group_cols) keep <- keep & data[[g]] == keys[[g]][i]
    for (m in metrics) {
      ss <- mean_se_values(data[[m]][keep])
      ans[[z]] <- cbind(keys[i, , drop = FALSE], Metric = m,
                        Mean = ss[["Mean"]], SD = ss[["SD"]],
                        SE = ss[["SE"]], N = ss[["N"]])
      z <- z + 1L
    }
  }
  do.call(rbind, ans)
}

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x
}

write_latex_table <- function(df, path, caption, label, align = NULL) {
  if (is.null(align)) align <- paste0("l", paste(rep("c", ncol(df) - 1L), collapse = ""))
  con <- base::file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines("\\begin{table}[!ht]", con)
  writeLines("\\centering", con)
  writeLines(paste0("\\caption{", caption, "}"), con)
  writeLines(paste0("\\label{", label, "}"), con)
  writeLines(paste0("\\begin{tabular}{", align, "}"), con)
  writeLines("\\hline", con)
  header <- gsub("_", "\\\\_", names(df), fixed = TRUE)
  writeLines(paste(paste0(header, collapse = " & "), "\\\\"), con)
  writeLines("\\hline", con)
  for (i in seq_len(nrow(df))) {
    row <- vapply(df[i, , drop = FALSE], as.character, character(1))
    row <- gsub("_", "\\\\_", row, fixed = TRUE)
    row <- gsub("%", "\\\\%", row, fixed = TRUE)
    writeLines(paste(paste(row, collapse = " & "), "\\\\"), con)
  }
  writeLines("\\hline", con)
  writeLines("\\end{tabular}", con)
  writeLines("\\end{table}", con)
}

save_plot_both <- function(plot, stem, width = 9, height = 6, dpi = 300) {
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in")
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height,
                  units = "in", dpi = dpi)
}

long_metrics <- function(data, id_cols, metrics) {
  out <- lapply(metrics, function(m) {
    z <- data[id_cols]
    z$Metric <- m
    z$Value <- data[[m]]
    z
  })
  do.call(rbind, out)
}
