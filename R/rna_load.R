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

#' setup_files
#'
#' Creates data.table/data.frame with file paths and associated metadata.
#'
#' @param wd Path to TAP output directory
#' @param pattern File pattern to match (.bam, .samon_quant, .ReadsPerGene.out.tab, etc)
#' @param var_map Variables that can be extracted from file names. File names are assumed to be delimitted by "_" and/or ".". See details.
#'
#' @return data.table/data.frame
#'
#' `var_map` encodes how metadata variables can be extracted from file names returned by matching the `pattern` in `wd`.
#'
#'  Let's use 6d_queen_A.Aligned.sortedByCoord.out.bam as an example file.
#'
#'  `var_map = c("day", "sex", "rep")` would extract 6d, queen, and A into columns "day", "sex", and "rep" respectively.
#'
#'  If you need to be more specific you could do:
#'
#'  `var_map = c("sex" = 2)` which would extract only the second item into column "sex".
#'
#'  Using this specification type we can get the same result as the first exmaple with:
#'
#'  `var_map = c("day" = 1, "sex" = 2, "rep" = 3)`
#'
#' @rawNamespace import(data.table, except = c(shift, first, second, last, melt))
#' @export
#' @rdname tapfiles
#' @import ggplot2
#' @import GenomicRanges
setup_files = function(wd, pattern, var_map = NULL){
  files = dir(wd, pattern = pattern, full.names = TRUE)
  dt = data.table::data.table(file = files)
  if(!is.null(var_map)){
    if(is.character(var_map)){
      tmp = seq_along(var_map)
      names(tmp) = var_map
      var_map = tmp
    }
    # dt[, names(var_map) := data.table::tstrsplit(basename(file), "[_\\.]", keep = var_map)]
    data.table::set(dt, j = names(var_map), value = data.table::tstrsplit(basename(dt$file), "[_\\.]", keep = var_map))
    # dt[, c(names(var_map)) := data.table::tstrsplit(basename(file), "[_\\.]", keep = var_map)]
    # dt[, `:=`(names(var_map), data.table::tstrsplit(basename(file), "[_\\.]", keep = var_map))]
  }
  dt[]
}



#' @export
#' @rdname tapfiles
#' @examples
#' wd = "/slipstream_old/home/joeboyd/R_workspace/SFtapfly.data"
#' setup_bam_files(wd)
#' #usage of var_map
#' setup_bam_files(wd, var_map = c("cell" = 1, "temp" = 2, "rep" = 3))
#' # general form
#' setup_bam_files(wd, var_map = c("cell" = 1, "temp" = 2, "rep" = 3), pattern = 'Aligned.sortedByCoord.out.bam$')
setup_bam_files = function(wd, var_map = NULL){
  dt = setup_files(wd, "Aligned.sortedByCoord.out.bam$", var_map)
  dt[, mapped_reads := get_mapped_reads(file), .(file)]
  dt[]
}

#' @export
#' @rdname tapfiles
setup_bam_files.transcriptome = function(wd, var_map = NULL){
  dt = setup_files(wd, "Aligned.toTranscriptome.out.bam$", var_map)
  dt[, mapped_reads := get_mapped_reads(file), .(file)]
  dt[]
}

#' @export
#' @rdname tapfiles
setup_peak_files = function(wd, var_map = NULL){
  setup_files(wd, ".narrowPeak$", var_map)
}

#' @export
#' @rdname tapfiles
setup_count_files = function(wd, var_map = NULL){
  setup_files(wd, ".ReadsPerGene.out.tab$", var_map)
}

#' Creates data.table/data.frame with file paths and associated metadata.
#'
#' @param wd
#' @param lib_type
#' @param name_composition
#' @param just_check_library_type
#' @param gtf_file
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
load_counts = function(wd, lib_type = NULL, name_composition = NULL, just_check_library_type = FALSE, gtf_file = NULL){
  use_gene_name = ifelse(is.null(gtf_file), FALSE, TRUE)
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
  if(is.null(lib_type)){
    lib_type = unique(sapply(dt$file, guess_lib_from_file))
    if(length(lib_type) > 1){
      stop("Multiple potential library types detected: ", paste(lib_type, collapse = ", "),
           "\nUse TAPhelpR::guess_lib_from_file to investigate and then run load_counts() again but manually specify lib_type.")
    }
  }
  message("Loading all files using library type: ", lib_type)
  mat = load_matrix_from_ReadsPerGene.out.tab(dt$file, lib_type = lib_type)
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
