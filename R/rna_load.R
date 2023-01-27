#' get_mapped_reads
#'
#' @param bam_file A .bam file. Matching .bai file must exist.
#'
#' @return The number of mapped reads in bam file.
#' @export
#'
#' @examples
#' bam_file = system.file("extdata/MCF10A_CTCF_R1.100peaks.bam", package = "TAPhelpR")
#' get_mapped_reads(bam_file)
get_mapped_reads = function(bam_file){
  stopifnot(file.exists(bam_file))
  stopifnot(file.exists(paste0(bam_file, ".bai")))
  stats = Rsamtools::idxstatsBam(bam_file)
  sum(stats[,3])
}

#' locate_files
#'
#' @param wd
#' @param pattern
#'
#' @return
#' @export
#' @rdname tap
#'
#' @examples
#' wd = "/slipstream_old/home/joeboyd/R_workspace/SFtapfly.data"
#' locate_files(wd, "final.out$")
locate_files = function(wd, pattern){
  dir(wd, pattern = pattern, full.names = TRUE)
}

.setup_files = function(wd, pattern, var_map = NULL){
  files = locate_files(wd, pattern)
  dt = data.table::data.table(file = files)
  if(!is.null(var_map)){
    dt[, names(var_map) := data.table::tstrsplit(basename(file), "[_\\.]", keep = var_map)]
  }
  dt[]
}

#' Creates data.table/data.frame with file paths and associated metadata.
#'
#' @return data.table/data.frame
#' @export
#' @rdname tap
#'
#' @rawNamespace import(data.table, except = c(shift, first, second, last, melt))
#' @import ggplot2
#' @import GenomicRanges
#'
#' @examples
#' wd = "/slipstream_old/home/joeboyd/R_workspace/SFtapfly.data"
#' setup_bam_files(wd)
#' #usage of var_map
#' setup_bam_files(wd, var_map = c("cell" = 1, "temp" = 2, "rep" = 3))
setup_bam_files = function(wd, var_map = NULL){
  dt = .setup_files(wd, ".bam$", var_map)
  dt[, mapped_reads := get_mapped_reads(file), .(file)]
  dt[]
}

#' @export
#' @rdname tap
#'
#' @examples
#' wd = "/slipstream_old/home/joeboyd/R_workspace/"
#' setup_peak_files(wd)
#' setup_peak_files(wd, var_map = c("cell" = 1, "mark" = 2, "rep" = 3))
setup_peak_files = function(wd, var_map = NULL){
  .setup_files(wd, ".narrowPeak$", var_map)
}

#' Creates data.table/data.frame with file paths and associated metadata.
#'
#' @return data.table/data.frame
#' @export
#' @rdname tap
#'
#' @rawNamespace import(data.table, except = c(shift, first, second, last, melt))
#' @import ggplot2
#' @import GenomicRanges
#'
#' @examples
#' wd = "/slipstream_old/home/joeboyd/R_workspace/SFtapfly.data"
#' setup_count_files(wd)
#'
#' setup_count_files(wd, var_map = c("cell" = 1, "temp" = 2, "rep" = 3))
setup_count_files = function(wd, var_map = NULL){
  .setup_files(wd, ".ReadsPerGene.out.tab$", var_map)
}

#' Creates data.table/data.frame with file paths and associated metadata.
#'
#' @return wide matrix of RNA counts
#' @export
#' @rdname tap
#'
#' @rawNamespace import(data.table, except = c(shift, first, second, last, melt))
#' @import ggplot2
#' @import GenomicRanges
#'
#' @examples
#' wd = "/slipstream_old/home/joeboyd/R_workspace/SFtapfly.data"
#' load_counts(wd)
#' load_counts(wd, name_composition = c(1, 2, 3, 4))
load_counts = function(wd, name_composition = NULL, just_check_library_type = FALSE, use_gene_name = TRUE, gtf_file = NULL){

  if(is.null(name_composition)){
    setup_var_map = NULL
  }else{
    setup_var_map = name_composition
    var_names = paste0("var", seq_along(name_composition))
    names(setup_var_map) = var_names
  }
  dt = setup_count_files(wd, var_map = setup_var_map)
  if(is.null(name_composition)){
    dt$name = sub(".ReadsPerGene.out.tab", "", basename(dt$file))
  }else{
    setup_var_map = name_composition
    paste0("var", seq_along(name_composition))
    found_names = apply(dt[, var_names, with = FALSE], 1, function(x)paste(x, collapse = "_"))
    dt$name = found_names
  }
  # checked and all are second
  if(just_check_library_type){
    return(sapply(dt$file, guess_lib_from_file))
  }
  toload = dt$file
  names(toload) = basename(dt$file)
  if(any(duplicated(names(toload)))){
    names(toload) = paste0(names(toload), "_", table(names(toload)))
  }
  mat = load_matrix_from_ReadsPerGene.out.tab(dt$file, lib_type = "second")
  colnames(mat) = dt$name
  if(use_gene_name){
    mat = aggregate_by_gene_name(mat, gtf_file = gtf_file)
  }
  mat
}

#' xeno_rna.load_norm_counts
#'
#' @return wide matrix of RNA counts, RPM normalized after trimming .01 quantile
#' @export
#'
#' @examples
#' load_norm_counts()
load_norm_counts = function(wd, name_composition = NULL, use_gene_name = TRUE){
  raw_counts = load_counts(use_gene_name = use_gene_name, name_composition = name_composition)
  # apply(raw_counts, 2, quantile, probs = .99)
  rpm_counts = apply(
    raw_counts,
    2,
    function(x){
      length(x)
      x.counted = x[x > 0]
      length(x.counted)
      x.counted = x.counted[x.counted <= quantile(x.counted, .99)]
      length(x.counted)
      x/sum(x.counted)*1e6
    })
  # colSums(rpm_counts)/1e6
  rpm_counts
}

#' load_gene_reference
#'
#' @return
#' @export
#'
#' @examples
load_gene_reference = function(gtf_file){
  ref_gr = rtracklayer::import.gff(gtf_file, feature.type = "gene")
  names(ref_gr) = ref_gr$gene_id
  ref_gr
}

#' load_exon_reference
#'
#' @return
#' @export
#'
#' @examples
load_exon_reference = function(gtf_file){
  ref_gr = rtracklayer::import.gff(gtf_file, feature.type = "exon")
  names(ref_gr) = ref_gr$gene_id
  ref_gr
}
