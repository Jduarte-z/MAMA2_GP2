#!/usr/bin/env Rscript
# Forest plots for MAMA loci: one row per population GWAS, plus GWAMA RE and/or
# MR-MEGA rows, with a side table (EAF, N, p).
#
# MR-MEGA BETA/SE are PC-axis coefficients, not a pooled effect, so the OR/CI on
# the MR-MEGA row comes from GWAMA RE; p, EAF and N come from MR-MEGA.
# GWAMA's OR_se is not the SE of log(OR), so the SE is recovered from the 95% CI
# and checked against Z.


# one list() per plot. snpid as CHR:POS:REF:ALT
# meta: any of "RE", "MR-MEGA", or character(0) for populations only
panels <- list(

  list(
    gene       = "UCHL1",
    snpid      = "4:41239695:C:T",    # PLACEHOLDER
    provenance = "Novel locus Random Effects",
    meta       = c("RE")
  ),

  list(
    gene       = "UCHL1",
    snpid      = "4:41239353:G:A",    # PLACEHOLDER
    provenance = "Novel locus MR-MEGA",
    meta       = c("RE")
  ),

  list(
    gene       = "MFSD6",
    snpid      = "2:190413615:G:A",   # PLACEHOLDER
    provenance = "Novel locus, random effects",
    meta       = c("RE")
  ),
  list(
    gene       = "GRM7",
    snpid      = "3:7024398:G:A",   # PLACEHOLDER
    provenance = "Fine-mapped singleton, 95% credible set",
    meta       = c("RE")
  ),
  list(
    gene       = "SCARB2",
    snpid      = "4:76213633:C:T",   # PLACEHOLDER
    provenance = "Fine-mapped singleton, 95% credible set",
    meta       = c("RE")
  ),
  list(
    gene       = "SLC18B1",
    snpid      = "6:132797077:A:G",   # PLACEHOLDER
    provenance = "Fine-mapped singleton, 95% credible set",
    meta       = c("RE")
  ),	
  list(
    gene       = "LRRK2",
    snpid      = "12:40227006:C:G",   # PLACEHOLDER
    provenance = "Fine-mapped singleton, 95% credible set",
    meta       = c("RE")
  )
	
)

  # list(
  #   gene       = "",
  #   snpid      = "",    # PLACEHOLDER
  #   provenance = "Known locus",
  #   meta       = c("RE")
  # )

#)

# inputs (both lists must name the same files)
mr_mega_input_list <- "../full_mr-mega_10pops_g2019sAJcond/MR-MEGA_input.txt"
gwama_input_list   <- "../full_random-effects_10pops_g2019sAJcond/gwama_input_list.txt"

# base dir for relative paths in the lists; NULL = dir of each list file
input_list_base_dir <- NULL

re_sumstats      <- "../full_random-effects_10pops_g2019sAJcond/r12_mama_randomEff.maf1.gwaslab.tsv.gz"
mr_mega_sumstats <- "../full_mr-mega_10pops_g2019sAJcond/r12_mama_mr_mega.maf1.gwaslab.tsv.gz"

# GP2 ancestry PCA palette; population is taken from the file name
population_colours <- c(
  AAC  = "#8E8E37",
  AFR  = "#7DC5E9",
  AJ   = "#D25A0A",
  AMR  = "#2C1F7A",
  CAH  = "#3BA69C",
  CAS  = "#7C214A",
  EAS  = "#D8C673",
  EUR  = "#147135",
  FIN  = "#EEE04D",
  MDE  = "#5B1309",
  SAS  = "#C45C6C",
  META = "#000000"
)

meta_row_labels <- c(
  "RE"      = "Random effects",
  "MR-MEGA" = "MR-MEGA*"
)

mr_mega_footnote <- paste(
  "* MR-MEGA row: OR and 95% CI from the GWAMA random-effects meta-analysis",
  "(MR-MEGA betas are PC-axis coefficients); EAF, N and p-value from MR-MEGA."
)

# appearance / output
missing_label  <- "--"
psignif        <- 0.05      # p >= psignif drawn hollow
show_legend    <- TRUE
table_rel_width <- 2.2      # forest:table = 5:table_rel_width
panel_width_in <- 11
height_base_in <- 1.9
height_row_in  <- 0.32
footnote_in    <- 0.35
dpi            <- 600
combined_ncol  <- 1
out_dir        <- "."


required_pkgs <- c(
  tidyverse    = "install.packages(\"tidyverse\")",
  patchwork    = "install.packages(\"patchwork\")",
  ggforestplot = "install.packages(\"remotes\"); remotes::install_github(\"NightingaleHealth/ggforestplot\")"
)

missing_pkgs <- names(required_pkgs)[
  !vapply(names(required_pkgs), requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_pkgs, collapse = ", "), "\n",
    "Install with:\n  ", paste(required_pkgs[missing_pkgs], collapse = "\n  "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggforestplot)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")


# check config
if (length(panels) == 0) stop("No panels configured.", call. = FALSE)

for (i in seq_along(panels)) {
  p <- panels[[i]]
  needed <- c("gene", "snpid", "provenance", "meta")
  miss <- setdiff(needed, names(p))
  if (length(miss) > 0) {
    stop(sprintf("Panel %d is missing field(s): %s", i, paste(miss, collapse = ", ")),
         call. = FALSE)
  }
  if (!is.character(p$snpid) || length(p$snpid) != 1 || !nzchar(p$snpid)) {
    stop(sprintf("Panel %d: snpid must be a single non-empty string.", i), call. = FALSE)
  }
  bad_meta <- setdiff(p$meta %||% character(0), names(meta_row_labels))
  if (length(bad_meta) > 0) {
    stop(sprintf("Panel %d: unknown meta value(s): %s. Allowed: %s", i,
                 paste(bad_meta, collapse = ", "),
                 paste(names(meta_row_labels), collapse = ", ")),
         call. = FALSE)
  }
}

for (f in c(mr_mega_input_list, gwama_input_list, re_sumstats, mr_mega_sumstats)) {
  if (!file.exists(f)) stop("File not found: ", f, call. = FALSE)
}

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


read_input_list <- function(list_file, base_dir = NULL) {
  base_dir <- base_dir %||% dirname(normalizePath(list_file))

  paths <- readLines(list_file, warn = FALSE) %>%
    str_trim() %>%
    discard(~ .x == "" || str_starts(.x, "#")) %>%
    path.expand()

  resolved <- ifelse(str_starts(paths, "/"), paths, file.path(base_dir, paths))

  not_found <- resolved[!file.exists(resolved)]
  if (length(not_found) > 0) {
    stop("Files listed in ", list_file, " not found:\n  ",
         paste(not_found, collapse = "\n  "),
         "\n(relative paths resolved against: ", base_dir, ")", call. = FALSE)
  }

  normalizePath(resolved)
}

# decompressor from magic bytes
decompress_cmd <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", n = 6)

  starts_with_bytes <- function(bytes) {
    length(magic) >= length(bytes) && all(magic[seq_along(bytes)] == as.raw(bytes))
  }

  cmd <- if (starts_with_bytes(c(0x1f, 0x8b))) {
    "gzip -cd"   # also bgzip
  } else if (starts_with_bytes(c(0x42, 0x5a, 0x68))) {
    "bzip2 -cd"
  } else if (starts_with_bytes(c(0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00))) {
    "xz -cd"
  } else if (starts_with_bytes(c(0x28, 0xb5, 0x2f, 0xfd))) {
    "zstd -dcq"
  } else {
    "cat"
  }

  tool <- word(cmd, 1)
  if (!nzchar(Sys.which(tool))) {
    stop("'", tool, "' is needed to read ", path, " but is not on PATH.", call. = FALSE)
  }
  cmd
}

# grep -w -F prefilter, then exact match on SNPID
fetch_rows <- function(path, ids_file, ids) {
  dcmp <- decompress_cmd(path)

  header <- system(sprintf("%s %s | head -n 1", dcmp, shQuote(path)), intern = TRUE)
  if (length(header) == 0) stop("Empty file: ", path, call. = FALSE)

  cmd <- sprintf("%s %s | LC_ALL=C grep -w -F -f %s",
                 dcmp, shQuote(path), shQuote(ids_file))
  hits <- suppressWarnings(system(cmd, intern = TRUE))

  status <- attr(hits, "status") %||% 0L
  if (status > 1) stop("grep failed (exit ", status, ") on ", path, call. = FALSE)

  df <- read_tsv(
    I(paste0(paste(c(header, hits), collapse = "\n"), "\n")),
    col_types = cols(.default = col_character()),
    na = c("", "NA", "NaN", "nan", "NULL", "."),
    progress = FALSE
  )

  if (!"SNPID" %in% names(df)) stop("No SNPID column in ", path, call. = FALSE)

  df <- filter(df, SNPID %in% ids)

  dups <- df %>% count(SNPID) %>% filter(n > 1)
  if (nrow(dups) > 0) {
    warning("Duplicated SNPID(s) in ", path, ": ",
            paste(dups$SNPID, collapse = ", "), " (keeping the first row).",
            call. = FALSE)
    df <- distinct(df, SNPID, .keep_all = TRUE)
  }
  df
}

require_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(what, " is missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(df)
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

fmt_or_missing <- function(x, formatter) {
  if_else(is.na(x), missing_label, formatter(x))
}


# input lists
mrmega_files <- read_input_list(mr_mega_input_list, input_list_base_dir)
gwama_files  <- read_input_list(gwama_input_list,   input_list_base_dir)

if (anyDuplicated(mrmega_files)) {
  stop("Duplicated file(s) in ", mr_mega_input_list, call. = FALSE)
}

if (!setequal(mrmega_files, gwama_files)) {
  stop(
    "MR-MEGA and GWAMA input lists differ.\n",
    "  Only in ", mr_mega_input_list, ":\n    ",
    paste(setdiff(mrmega_files, gwama_files), collapse = "\n    "), "\n",
    "  Only in ", gwama_input_list, ":\n    ",
    paste(setdiff(gwama_files, mrmega_files), collapse = "\n    "),
    call. = FALSE
  )
}

if (!identical(mrmega_files, gwama_files)) {
  warning("MR-MEGA and GWAMA input lists contain the same files in a different order.",
          call. = FALSE)
}

message("Input lists match (", length(mrmega_files), " studies).")

population_codes <- setdiff(names(population_colours), "META")

label_from_path <- function(path) {
  b <- toupper(basename(path))
  pattern <- str_c("(^|[^A-Z])", population_codes, "([^A-Z]|$)")
  hits <- population_codes[str_detect(b, pattern)]
  if (length(hits) != 1) {
    stop("Could not assign exactly one population to ", basename(path),
         " (matches: ", if (length(hits)) paste(hits, collapse = ", ") else "none",
         "). Add/adjust codes in population_colours.", call. = FALSE)
  }
  hits
}

study_manifest <- tibble(
  file  = mrmega_files,
  study = map_chr(mrmega_files, label_from_path)
) %>%
  arrange(study)

if (anyDuplicated(study_manifest$study)) {
  stop("Two input files map to the same population label: ",
       paste(unique(study_manifest$study[duplicated(study_manifest$study)]), collapse = ", "),
       call. = FALSE)
}

message("Populations: ", paste(study_manifest$study, collapse = ", "))


# read variants
all_ids  <- unique(map_chr(panels, "snpid"))
ids_file <- tempfile(fileext = ".txt")
writeLines(all_ids, ids_file)

# populations
study_cols <- c("SNPID", "EA", "NEA", "EAF", "BETA", "SE", "P", "N")

study_hits <- study_manifest %>%
  mutate(rows = map2(file, study, function(f, s) {
    message("Reading ", s, ": ", f)
    fetch_rows(f, ids_file, all_ids) %>%
      require_cols(study_cols, f) %>%
      select(all_of(study_cols))
  })) %>%
  select(study, rows) %>%
  unnest(rows) %>%
  mutate(across(c(EAF, BETA, SE, P, N), as.numeric))

# GWAMA RE
re_cols <- c("SNPID", "EA", "NEA", "EAF", "OR", "OR_95L", "OR_95U", "Z", "P", "N")
z_crit  <- qnorm(0.975)

message("Reading random effects: ", re_sumstats)
re_hits <- fetch_rows(re_sumstats, ids_file, all_ids) %>%
  require_cols(re_cols, re_sumstats) %>%
  select(all_of(re_cols)) %>%
  mutate(across(c(EAF, OR, OR_95L, OR_95U, Z, P, N), as.numeric)) %>%
  mutate(
    BETA = log(OR),
    SE   = (log(OR_95U) - log(OR_95L)) / (2 * z_crit)
  )

# SE from CI vs GWAMA Z, 5% tolerance
se_check <- re_hits %>%
  filter(!is.na(Z), abs(Z) > 1) %>%
  mutate(rel_diff = abs(BETA / SE - Z) / abs(Z)) %>%
  filter(rel_diff > 0.05)

if (nrow(se_check) > 0) {
  warning("CI-derived SE disagrees with GWAMA Z for: ",
          paste(se_check$SNPID, collapse = ", "),
          ". Check how the random-effects CI was computed.", call. = FALSE)
}

# MR-MEGA
mm_cols <- c("SNPID", "EA", "NEA", "EAF", "P", "N")

message("Reading MR-MEGA: ", mr_mega_sumstats)
mm_hits <- fetch_rows(mr_mega_sumstats, ids_file, all_ids) %>%
  require_cols(mm_cols, mr_mega_sumstats) %>%
  select(all_of(mm_cols)) %>%
  mutate(across(c(EAF, P, N), as.numeric))

# allele and missing-SNP checks
allele_conflicts <- bind_rows(
  select(study_hits, SNPID, EA, NEA),
  select(re_hits,    SNPID, EA, NEA),
  select(mm_hits,    SNPID, EA, NEA)
) %>%
  distinct() %>%
  count(SNPID) %>%
  filter(n > 1)

if (nrow(allele_conflicts) > 0) {
  warning("EA/NEA differ across files for: ",
          paste(allele_conflicts$SNPID, collapse = ", "), call. = FALSE)
}

not_found <- setdiff(all_ids, c(study_hits$SNPID, re_hits$SNPID, mm_hits$SNPID))
if (length(not_found) > 0) {
  warning("SNPID(s) not found in any file: ", paste(not_found, collapse = ", "),
          call. = FALSE)
}


build_panel_df <- function(snpid, meta) {

  studies <- study_manifest %>%
    select(study) %>%
    left_join(filter(study_hits, SNPID == snpid), by = "study") %>%
    transmute(
      row_label    = study,
      colour_group = study,
      row_type     = "Study",
      estimate     = BETA,
      se           = SE,
      pvalue       = P,
      n_ind        = N,
      eaf          = EAF
    )

  re_row <- filter(re_hits, SNPID == snpid)
  mm_row <- filter(mm_hits, SNPID == snpid)

  meta_rows <- map(meta, function(m) {
    src <- if (m == "RE") re_row else mm_row
    has_src <- nrow(src) == 1
    has_re  <- nrow(re_row) == 1
    tibble(
      row_label    = meta_row_labels[[m]],
      colour_group = "META",
      row_type     = "Meta",
      # effect always from GWAMA RE
      estimate     = if (has_src && has_re) re_row$BETA else NA_real_,
      se           = if (has_src && has_re) re_row$SE   else NA_real_,
      pvalue       = if (has_src) src$P   else NA_real_,
      n_ind        = if (has_src) src$N   else NA_real_,
      eaf          = if (has_src) src$EAF else NA_real_
    )
  })

  bind_rows(studies, list_rbind(meta_rows))
}

make_panel <- function(panel) {

  meta <- panel$meta %||% character(0)
  df   <- build_panel_df(panel$snpid, meta)

  if (all(is.na(df$estimate))) {
    warning(sprintf("No effect estimates for %s (%s) - panel skipped.",
                    panel$gene, panel$snpid), call. = FALSE)
    return(NULL)
  }

  missing_colours <- setdiff(unique(df$colour_group), names(population_colours))
  if (length(missing_colours) > 0) {
    stop("No colour defined for: ", paste(missing_colours, collapse = ", "), call. = FALSE)
  }

  df <- mutate(df, row_label = fct_inorder(row_label))

  show_footnote <- "MR-MEGA" %in% meta
  title_text <- sprintf("%s  %s  |  %s", panel$gene, panel$snpid, panel$provenance)

  p_forest <- ggforestplot::forestplot(
    df       = df,
    name     = row_label,
    estimate = estimate,
    se       = se,
    pvalue   = pvalue,
    psignif  = psignif,
    logodds  = TRUE,
    colour   = colour_group,
    shape    = row_type
  ) +
    scale_colour_manual(values = population_colours) +
    scale_shape_manual(values = c(Study = 21, Meta = 23), guide = "none") +
    guides(colour = guide_legend(reverse = FALSE, nrow = 2)) +
    labs(
      x       = "OR (95% CI)",
      colour  = "Population",
      caption = if (show_footnote) str_wrap(mr_mega_footnote, 110) else NULL
    ) +
    theme_bw() +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      axis.title.y    = element_blank(),
      plot.caption    = element_text(hjust = 0, size = 7.5, colour = "grey30")
    )

  # table, same row order as the forest
  df_table <- df %>%
    mutate(
      row_label = fct_rev(row_label),
      face      = if_else(row_type == "Meta", "bold", "plain"),
      EAF_lab   = fmt_or_missing(eaf,    function(x) sprintf("%.3f", x)),
      N_lab     = fmt_or_missing(n_ind,  function(x) formatC(round(x), format = "d", big.mark = ",")),
      P_lab     = fmt_or_missing(pvalue, function(x) sprintf("%.2e", x))
    )

  p_table <- ggplot(df_table, aes(y = row_label)) +
    ggforestplot::geom_stripes() +
    geom_text(aes(x = 0, label = EAF_lab, fontface = face), size = 3.2) +
    geom_text(aes(x = 1, label = N_lab,   fontface = face), size = 3.2) +
    geom_text(aes(x = 2, label = P_lab,   fontface = face), size = 3.2) +
    scale_x_continuous(
      limits   = c(-0.5, 2.5),
      breaks   = 0:2,
      labels   = c("EAF", "N", "p-value"),
      position = "top"
    ) +
    theme_void() +
    theme(
      axis.text.x.top = element_text(face = "bold", size = 9),
      plot.margin     = margin(5, 5, 5, 5)
    )

  p_title <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = title_text,
             fontface = "bold", size = 4.6) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_void()

  body <- (p_forest | p_table) + plot_layout(widths = c(5, table_rel_width))

  panel_plot <- p_title / body +
    plot_layout(heights = unit(c(0.45, 1), c("in", "null")))

  height <- height_base_in + height_row_in * nrow(df) + if (show_footnote) footnote_in else 0

  fname <- sprintf("forest_%s_%s_%s.png",
                   safe_name(panel$gene), safe_name(panel$snpid), run_stamp)

  ggsave(file.path(out_dir, fname), plot = panel_plot,
         width = panel_width_in, height = height, dpi = dpi, bg = "white")
  message("Wrote ", fname)

  list(plot = panel_plot, height = height)
}


built <- panels %>% map(make_panel) %>% compact()

if (length(built) >= 2) {
  n_panels <- length(built)
  n_rows   <- ceiling(n_panels / combined_ncol)
  heights  <- map_dbl(built, "height")

  row_heights <- map_dbl(seq_len(n_rows), function(r) {
    idx <- ((r - 1) * combined_ncol + 1):min(r * combined_ncol, n_panels)
    max(heights[idx])
  })

  combined <- wrap_plots(map(built, "plot"), ncol = combined_ncol, heights = row_heights)

  fname <- sprintf("forest_combined_%s.png", run_stamp)
  ggsave(file.path(out_dir, fname), plot = combined,
         width = panel_width_in * combined_ncol, height = sum(row_heights),
         dpi = dpi, bg = "white", limitsize = FALSE)
  message("Wrote ", fname)
} else {
  message("Fewer than two panels were drawn - combined figure not written.")
}

unlink(ids_file)
