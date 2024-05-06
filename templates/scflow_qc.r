#!/usr/bin/env Rscript
# Perform quality-control on a feature-barcode matrix with scflow
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

options(mc.cores = future::availableCores())

##  ............................................................................
##  Load packages                                                           ####
library(argparse)
library(scFlow)

##  ............................................................................
##  Parse command-line arguments                                            ####


#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse
parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- unlist(lapply(args_list, function(y) strsplit(y, ' +')), recursive = FALSE)

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}


args <- parse_args("!{params.options.args}")


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Pre-process args                                                        ####

options("scflow_species" = args$species)

args[startsWith(names(args), "drop_")] <-
  as.logical(args[startsWith(names(args), "drop_")])
args$max_library_size <- ifelse(
  args$max_library_size == "adaptive",
  args$max_library_size,
  as.numeric(as.character(args$max_library_size))
)
args$max_features <- ifelse(
  args$max_features == "adaptive",
  args$max_features,
  as.numeric(as.character(args$max_features))
)
args$max_mito <- ifelse(
  args$max_mito == "adaptive",
  args$max_mito,
  as.numeric(as.character(args$max_mito))
)
args$pK <- if (toupper(args$pK) == "NULL") NULL else {
  as.numeric(as.character(args$pK))
}
args$dpk <- if (toupper(args$dpk) == "NULL") NULL else {
  as.numeric(as.character(args$dpk))
}

if (toupper(args$retain) == "NULL") {
  args$retain <- NULL
} else if (toupper(args$retain) == "AUTO") {
  args$retain <- "auto"
} else {
  args$retain <- as.numeric(as.character(args$retain))
}

args$find_singlets <- as.logical(args$find_singlets)
args$vars_to_regress_out <- strsplit(args$vars_to_regress_out, ",")[[1]]
args <- purrr::map(args, function(x) {
  if (length(x) == 1) {
    if (toupper(x) == "TRUE") return(TRUE)
    if (toupper(x) == "FALSE") return(FALSE)
    if (toupper(x) == "NULL") return(NULL)
  }
  return(x)
})

if (!is.null(args$factor_vars)) {
  args$factor_vars <- strsplit(args$factor_vars, ",")[[1]]
  col_classes <- rep("factor", length(args$factor_vars))
  names(col_classes) <- args$factor_vars
} else {
  col_classes <- NA
}

args$nmads <- as.numeric(as.character(args$nmads))
args$min_library_size <- as.numeric(as.character(args$min_library_size))
args$min_features <- as.numeric(as.character(args$min_features))
args$min_ribo <- as.numeric(as.character(args$min_ribo))
args$max_ribo <- as.numeric(as.character(args$max_ribo))
args$min_counts <- as.numeric(as.character(args$min_counts))
args$min_cells <- as.numeric(as.character(args$min_cells))
args$pca_dims <- as.numeric(as.character(args$pca_dims))
args$var_features <- as.numeric(as.character(args$var_features))
args$doublet_rate <- as.numeric(as.character(args$doublet_rate))
args$lower <- as.numeric(as.character(args$lower))
args$alpha_cutoff <- as.numeric(as.character(args$alpha_cutoff))
args$niters <- as.numeric(as.character(args$niters))


##  ............................................................................
##  Start QC                                                                ####

cli::boxx(paste0("Analysing: ", "!{key}"), float = "center")

mat <- scFlow::read_sparse_matrix("!{mat_path}")

metadata <- read_metadata(
  unique_key = "!{key}",
  key_colname = args$key_colname,
  samplesheet_path = "!{input}",
  col_classes = col_classes
)

sce <- generate_sce(mat, metadata)

rm(mat, metadata)

if (args$find_cells) {
  sce <- find_cells(
    sce,
    lower = args$lower,
    retain = args$retain,
    alpha_cutoff = args$alpha_cutoff,
    niters = args$niters
  )
}

sce <- annotate_sce(
  sce = sce,
  min_library_size = args$min_library_size,
  max_library_size = args$max_library_size,
  min_features = args$min_features,
  max_features = args$max_features,
  max_mito = args$max_mito,
  min_ribo = args$min_ribo,
  max_ribo = args$max_ribo,
  min_counts = args$min_counts,
  min_cells = args$min_cells,
  drop_unmapped = args$drop_unmapped,
  drop_mito = args$drop_mito,
  drop_ribo = args$drop_ribo,
  nmads = args$nmads,
  annotate_genes = TRUE,
  annotate_cells = TRUE,
  ensembl_mapping_file = "!{ensembl_mappings}",
  species = args$species
)

sce <- filter_sce(
  sce,
  filter_genes = TRUE,
  filter_cells = TRUE
)

if (args$find_singlets) {
  sce <- find_singlets(
    sce = sce,
    singlet_find_method = args$singlets_method,
    vars_to_regress_out = args$vars_to_regress_out,
    pca_dims = args$pca_dims,
    var_features = args$var_features,
    doublet_rate = args$doublet_rate,
    dpk = args$dpk,
    pK = args$pK,
    num.cores = future::availableCores()
  )
  sce <- filter_sce(
    sce,
    filter_genes = TRUE,
    filter_cells = TRUE
  )
}


sce <- sce[ , sce$total_counts >= args$min_library_size]
sce <- sce[ , sce$total_features_by_counts >= args$min_features]


dir.create(file.path(getwd(), "qc_report"))

report_qc_sce(
  sce = sce,
  #report_folder_path = file.path(getwd(), "qc_report"),
  report_folder_path = file.path(getwd()),
  report_file = paste0("!{key}", "_scflow_qc_report")
)

print("Analysis complete, saving outputs..")

##  ............................................................................
##  Save Outputs                                                            ####

# Save SingleCellExperiment
write_sce(
  sce = sce,
  folder_path = file.path(getwd(), paste0("!{key}", "_sce"))
)

new_dirs <- c(
  "qc_plot_data",
  "qc_summary",
  "qc_plots")

#make dirs
purrr::walk(new_dirs, ~ dir.create(file.path(getwd(), .)))

# Save QC plots (tables)
for (df in names(sce@metadata$qc_plot_data)) {
  write.table(
    sce@metadata$qc_plot_data[[df]],
    file.path(getwd(), "qc_plot_data",
              paste0("!{key}", "_", df, ".tsv")),
    sep = "\t",
    col.names = TRUE, row.names = FALSE)
}

# Save QC summary table
write.table(
  cbind(sce@metadata$metadata, sce@metadata$qc_summary),
  file.path(getwd(), "qc_summary",
            paste0("!{key}", "_qc_summary.tsv")),
  sep = "\t",
  col.names = TRUE, row.names = FALSE)

# Save QC plots (images)
for (pname in names(sce@metadata$qc_plots)) {
  png(file.path(getwd(), "qc_plots",
                paste0("!{key}", "_", pname, ".png")),
      width = 247, height = 170, units = "mm", res = 600)
  print(sce@metadata$qc_plots[[pname]])
  dev.off()
}

# Save doublet finder plots, square
for (pname in names(sce@metadata$qc_plots$doublet_finder)) {
  png(file.path(getwd(), "qc_plots",
                paste0("!{key}", "_", pname, "_doublet_finder.png")),
      width = 170, height = 170, units = "mm", res = 600)
  print(sce@metadata$qc_plots$doublet_finder[[pname]])
  dev.off()
}

system("mkdir sce; mv !{key}_sce sce/")

##  ............................................................................
##  Clean up                                                                ####
