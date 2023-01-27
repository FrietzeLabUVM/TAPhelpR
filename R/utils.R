#' make_tidy_from_matrix
#'
#' @param wide_matrix
#'
#' @return
#' @export
#'
#' @examples
make_tidy_from_matrix = function(wide_matrix){
  dt = data.table::as.data.table(reshape2::melt(wide_matrix))
  colnames(dt) = c("gene_id", "sample", "count")
  if(any(grepl("rep", dt$sample))){
    dt[, c("cell", "rep") := data.table::tstrsplit(sample, "[_ ]")]
  }else{
    dt[, c("cell") := data.table::tstrsplit(sample, "[_ ]")]
  }
  dt = .add_name_and_group.no_treat(dt)
  dt
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
#' @param f
#' @param show_plots
#' @param miss_cutoff
#' @param hit_cutoff
#' @param split_cutoff
#'
#' @return
#' @export
#'
#' @examples
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
    labs(title = "Problematic reads per strandedness lib_type")

  dt_genes = dt[!grepl("N_", gene_id)]
  tdt_genes = reshape2::melt(dt_genes, variable.name = "lib", value.name = "count", id.vars = "gene_id")


  if(show_plots){
    p_genes = ggplot(tdt_genes, aes(x = lib, y = log10(count+1))) +
      geom_boxplot() +
      labs(title = "Reads mapped to genes per strandedness lib_type")

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
#' @param files
#' @param lib_type
#'
#' @return
#' @export
#'
#' @examples
load_matrix_from_ReadsPerGene.out.tab = function(files, lib_type = "guess"){
  stopifnot(lib_type %in% c("guess", "unstranded", "first", "second"))
  if(lib_type == "guess"){
    test_files = head(files)
    guesses = unique(sapply(files, guess_lib_from_file))
    if(length(guesses) > 1){
      stop("could not uniquely guess library type, please determine using guess_lib_from_file and manually specify lib_type.")
    }
    lib_type = guesses
  }
  stopifnot(lib_type %in% c("unstranded", "first", "second"))

  message("loading files...")
  dtl = pbmcapply::pbmclapply(files, mc.cores = 20, load_ReadsPerGene.out.tab, lib_type = lib_type)

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

load_ReadsPerGene.out.tab = function(f, lib_type){
  stopifnot(lib_type %in% c("unstranded", "first", "second"))
  obs_lib = guess_lib_from_file(f)
  if(obs_lib != lib_type){
    warning("file ", f, " had lib_type of ", obs_lib, ', not ', lib_type, " as expected.")
  }

  dt = data.table::fread(f, col.names = c("gene_id", "unstranded", "first", "second"))

  dt_genes = dt[!grepl("N_", gene_id)]
  dt_genes = dt_genes[, c("gene_id", lib_type), with = FALSE]
  colnames(dt_genes)[2] = sub(".ReadsPerGene.out.tab", "", basename(f))
  dt_genes
}

aggregate_by_gene_name = function(exp_cnt_mat, gtf_file){
  exp_cnt_dt = data.table::as.data.table(exp_cnt_mat, keep.rownames = "gene_id")
  create_matrix_from_data.table(exp_cnt_dt, gtf_file = gtf_file)
}

create_matrix_from_data.table = function(exp_cnt_dt, gtf_file){
  #### convert to gene_name ####
  if(is.character(gtf_file)){
    message("loading gtf...")
    ref_gr = rtracklayer::import.gff(gtf_file, feature.type = "gene")
  }else{
    ref_gr = gtf_file
  }
  stopifnot(is(ref_gr, class2 = "GRanges"))
  names(ref_gr) = ref_gr$gene_id

  message("aggregate by gene_name...")
  ref_dt = unique(data.table::data.table(gene_id = ref_gr$gene_id, gene_name = ref_gr$gene_name))
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
