#' make_tidy_from_matrix
#'
#' @param wide_matrix A wide matrix of count data with samples as columns and genes as rows
#' @param value_name Name for value attribute, i.e. "count".
#' @param row_name Name for row attribute, i.e. "gene_id".
#' @param column_name Name for column attribute, i.e. "name".
#'
#' @return A tidy data.table representation of input wide matrix
#' @export
#'
#' @examples
#' mat = matrix(1:9, ncol = 3, byrow = TRUE)
#' rownames(mat) = paste("gene", seq(nrow(mat)), sep = "_")
#' colnames(mat) = paste("sample", seq(ncol(mat)), sep = "_")
#' dt = make_tidy_from_matrix(mat)
#' dt
make_tidy_from_matrix = function(wide_matrix, value_name = "count", row_name = "gene_id", column_name = "name"){
  dt = data.table::as.data.table(reshape2::melt(wide_matrix))
  colnames(dt) = c(row_name, column_name, value_name)
  # dt = .add_name_and_group.no_treat(dt)
  dt[]
}

guess_lib = function(dt_stats,
                     miss_cutoff = 1.2,
                     hit_cutoff = 1.2,
                     split_cutoff = 1.5){
  dt_stats = dt_stats[grepl("N_", gene_id)]
  tdt_stats = data.table::as.data.table(reshape2::melt(dt_stats, variable.name = "lib", value.name = "count", id.vars = "gene_id"))

  un_count = tdt_stats[gene_id == "N_noFeature" & lib == "unstranded"]$count
  fir_count = tdt_stats[gene_id == "N_noFeature" & lib == "first"]$count
  sec_count = tdt_stats[gene_id == "N_noFeature" & lib == "second"]$count

  fir_ratio = fir_count / un_count
  sec_ratio = sec_count / un_count

  lib = "unrecognized"

  if(sec_ratio > miss_cutoff & fir_ratio < hit_cutoff){
    lib = "first"
  }else if(fir_ratio > miss_cutoff & sec_ratio < hit_cutoff){
    lib = "second"
  }else if(fir_ratio > split_cutoff & sec_ratio > split_cutoff){#this could be wrong
    lib = "unstranded"
  }
  lib
}

#' guess_lib_from_file
#'
#' Determines that library strandedness based on ratio of first to unstranded and second to unstranded.
#'
#' @param f count file to guess library type from.
#' @param show_plots if TRUE, diagnostic plots will be output to graphics device.
#' @param miss_cutoff Max allowed ratio of strand sensitive to unstranded. If strand is over this, it can't be that strand.
#' @param hit_cutoff Min allowed ratio of strand sensitive to unstranded. If strand is under this, it can't be that strand.
#' @param split_cutoff If both are over this ratio, it's unstranded.
#'
#' @return One of "unrecognized", "first", "second", or "unstranded" indicating library strandedness.
#' @export
#'
#' @examples
#' count_files = setup_count_files(example_honeybee_output())
#' guess_lib_from_file(count_files$file[1])
guess_lib_from_file = function(f,
                               show_plots = FALSE,
                               miss_cutoff = 1.2,
                               hit_cutoff = 1.2,
                               split_cutoff = 1.5){
  dt = data.table::fread(f, col.names = c("gene_id", "unstranded", "first", "second"))
  dt_stats = dt[grepl("N_", gene_id)]
  tdt_stats = reshape2::melt(dt_stats, variable.name = "lib", value.name = "count", id.vars = "gene_id")

  p_stats = ggplot(tdt_stats, aes(x = lib, y = count)) +
    geom_bar(stat = "identity") +
    facet_wrap(~gene_id, scales = "free_y") +
    scale_y_continuous(labels = function(x)paste(x/1e6, "M")) +
    labs(title = "Problematic reads per strandedness library_type")

  dt_genes = dt[!grepl("N_", gene_id)]
  tdt_genes = reshape2::melt(dt_genes, variable.name = "lib", value.name = "count", id.vars = "gene_id")


  if(show_plots){
    p_genes = ggplot(tdt_genes, aes(x = lib, y = log10(count+1))) +
      geom_boxplot() +
      labs(title = "Reads mapped to genes per strandedness library_type")

    plot(
      cowplot::plot_grid(ncol = 1, rel_heights = c(1, 12),
                         ggplot() + theme_void() + cowplot::draw_text("If first looks like second -> unstranded\nElse, one of first or second, whichever has fewer noFeature and more mapped to genes."),
                         cowplot::plot_grid(p_genes, p_stats))
    )

  }

  guess_lib(
    dt_stats,
    miss_cutoff = miss_cutoff,
    hit_cutoff = hit_cutoff,
    split_cutoff = split_cutoff
  )
}

#' load_matrix_from_ReadsPerGene.out.tab
#'
#' @param files Count files output by STAR.
#' @param library_type Strandedness of library. Should be one of "guess", "unstranded", "first", or "second". Default is "guess" and will result in attempting to automatically detect and apply library type.
#'
#' @return A matrix containing raw read counts for each file in files.
#' @export
#'
#' @examples
#' tap_out = "~/R_workspace.combined/TAPhelpR.data/honeybee_TAP_output"
#' count_files = setup_count_files(tap_out, variable_map = c("day", "sex", "rep"))
#' mat = load_matrix_from_ReadsPerGene.out.tab(count_files$file)
#' head(mat)
load_matrix_from_ReadsPerGene.out.tab = function(files, library_type = "guess"){
  stopifnot(library_type %in% c("guess", "unstranded", "first", "second"))
  if(is.data.frame(files)){
    if(!all(c("file", "name") %in% colnames(files))){
      stop("`files` must be either a named character vector or a data.frame containing variables 'file' and 'name'")
    }
    tmp = files$file
    names(tmp) = files$name
    file = tmp
  }
  if(library_type == "guess"){
    test_files = head(files)
    guesses = unique(sapply(files, guess_lib_from_file))
    if(length(guesses) > 1){
      stop("could not uniquely guess library type, please determine using guess_lib_from_file and manually specify library_type.")
    }
    library_type = guesses
  }
  stopifnot(library_type %in% c("unstranded", "first", "second"))

  message("loading files...")
  dtl = pbmcapply::pbmclapply(files, mc.cores = 20, load_ReadsPerGene.out.tab, library_type = library_type)

  message("assembling matrix...")
  mat = matrix(0, nrow = nrow(dtl[[1]]), ncol = length(dtl))
  rownames(mat) = dtl[[1]]$gene_id
  colnames(mat) = sapply(dtl, function(x)colnames(x)[2])
  mat[,1] = dtl[[1]][[2]]
  if(length(dtl) > 1){
    for(i in seq(2, length(dtl))){
      mat[,i] = dtl[[i]][[2]]
    }
  }
  mat
}

load_ReadsPerGene.out.tab = function(f, library_type){
  stopifnot(library_type %in% c("unstranded", "first", "second"))
  obs_lib = guess_lib_from_file(f)
  if(obs_lib != library_type){
    warning("file ", f, " had library_type of ", obs_lib, ', not ', library_type, " as expected.")
  }

  dt = data.table::fread(f, col.names = c("gene_id", "unstranded", "first", "second"))

  dt_genes = dt[!grepl("N_", gene_id)]
  dt_genes = dt_genes[, c("gene_id", library_type), with = FALSE]
  colnames(dt_genes)[2] = sub(".ReadsPerGene.out.tab", "", basename(f))
  dt_genes
}

aggregate_by_gene_name = function(exp_cnt_mat, gtf_file, name_attribute = "gene_name"){
  exp_cnt_dt = data.table::as.data.table(exp_cnt_mat, keep.rownames = "gene_id")
  create_matrix_from_data.table(exp_cnt_dt, gtf_file = gtf_file, name_attribute = name_attribute)
}

create_matrix_from_data.table = function(exp_cnt_dt, gtf_file, name_attribute = "gene_name"){
  #### convert to gene_name ####
  if(is.character(gtf_file)){
    message("loading gtf...")
    ref_gr = .load_ref(gtf_file, "gene")
  }else{
    ref_gr = gtf_file
  }
  stopifnot(is(ref_gr, class2 = "GRanges"))
  names(ref_gr) = ref_gr$gene_id

  message("aggregate by gene_name...")
  if(!name_attribute %in% colnames(mcols(ref_gr))){
    stop("name_attribute: ", name_attribute, " not found in gtf file. Possibe valid values included:\n",
         paste(collapse = "\n",
           setdiff(colnames(mcols(ref_gr)),
                   c("gene_id")))
    )
  }
  ref_dt = unique(data.table::data.table(gene_id = ref_gr$gene_id, gene_name = GenomicRanges::mcols(ref_gr)[[name_attribute]]))
  exp_cnt_dt.tidy = reshape2::melt(exp_cnt_dt, id.vars = "gene_id", value.name = "count", variable.name = "sample_id")
  stopifnot(all(unique(exp_cnt_dt.tidy$gene_id) %in% ref_dt$gene_id))
  #aggregate by gene_name
  exp_cnt_dt.tidy = merge(exp_cnt_dt.tidy, ref_dt, by = "gene_id")
  exp_cnt_dt.tidy = exp_cnt_dt.tidy[, list(count = sum(count)), list(sample_id, gene_name)]

  #### create matrix ####
  message("reshaping to final matrix...")
  exp_cnt_dt = reshape2::dcast(exp_cnt_dt.tidy, sample_id~gene_name, value.var = "count")
  exp_cnt_mat = as.matrix(exp_cnt_dt[, -1])
  rownames(exp_cnt_mat) = exp_cnt_dt$sample_id

  message("done!")
  t(exp_cnt_mat)
}

#from A to B
# x; value_var = NULL; A_cols = NULL; B_cols = NULL; A_name = NULL; B_name = NULL
calc_log2_FC_per_gene = function(x, value_var = NULL, A_cols = NULL, B_cols = NULL, A_name = NULL, B_name = NULL){
  if(is.null(value_var)) value_var = "value"
  if(is.matrix(x)){
    x = reshape2::melt(x, value.name = "value")
    colnames(x) = c("gene_id", "column_name", value_var)
  }
  cn = unique(x$column_name)
  if(is.null(A_cols)) A_cols = cn[1]
  if(is.null(B_cols)) B_cols = cn[2]
  if(is.null(A_name)) A_name = paste(A_cols, collapse = "&")
  if(is.null(B_name)) B_name = paste(B_cols, collapse = "&")
  FC_name = paste0("from_", A_name, "_to_", B_name)

}

symlink.exists = function(f){
  test_val = Sys.readlink(f)
  if(is.na(test_val)) return(FALSE)
  test_val != ""
}
