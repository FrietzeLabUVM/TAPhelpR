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
#' @param work_dir Path to TAP output directory
#' @param pattern File pattern to match (.bam, .samon_quant, .ReadsPerGene.out.tab, etc)
#' @param variable_map Variables that can be extracted from file names. File names are assumed to be delimitted by "_" and/or ".". See details.
#'
#' @return data.table/data.frame
#'
#' `variable_map` encodes how metadata variables can be extracted from file names returned by matching the `pattern` in `work_dir`.
#'
#'  Let's use 6d_queen_A.Aligned.sortedByCoord.out.bam as an example file.
#'
#'  `variable_map = c("day", "sex", "rep")` would extract 6d, queen, and A into columns "day", "sex", and "rep" respectively.
#'
#'  If you need to be more specific you could do:
#'
#'  `variable_map = c("sex" = 2)` which would extract only the second item into column "sex".
#'
#'  Using this specification type we can get the same result as the first exmaple with:
#'
#'  `variable_map = c("day" = 1, "sex" = 2, "rep" = 3)`
#'
#' @rawNamespace import(data.table, except = c(shift, first, second, last, melt))
#' @export
#' @rdname tapfiles
#' @import ggplot2
#' @import GenomicRanges
setup_files = function(work_dir, pattern, variable_map = NULL){
  files = dir(work_dir, pattern = pattern, full.names = TRUE)
  dt = data.table::data.table(file = files)
  if(!is.null(variable_map)){
    if(is.character(variable_map)){
      tmp = seq_along(variable_map)
      names(tmp) = variable_map
      variable_map = tmp
    }
    data.table::set(dt, j = names(variable_map), value = data.table::tstrsplit(basename(dt$file), "[_\\.]", keep = variable_map))
    if(!"name" %in% names(variable_map)){
      data.table::set(dt, j = "name", value = apply(dt[, names(variable_map), with = FALSE], 1, paste, collapse = "_"))
    }
  }
  dt[]
}



#' @export
#' @rdname tapfiles
#' @examples
#' work_dir = "/slipstream_old/home/joeboyd/R_workspace/SFtapfly.data"
#' setup_bam_files(work_dir)
#' #usage of variable_map
#' setup_bam_files(work_dir, variable_map = c("cell" = 1, "temp" = 2, "rep" = 3))
#' # general form
#' setup_bam_files(work_dir, variable_map = c("cell" = 1, "temp" = 2, "rep" = 3), pattern = 'Aligned.sortedByCoord.out.bam$')
setup_bam_files = function(work_dir, variable_map = NULL){
  dt = setup_files(work_dir, "Aligned.sortedByCoord.out.bam$", variable_map)
  dt[, mapped_reads := get_mapped_reads(file), .(file)]
  dt[]
}

#' @export
#' @rdname tapfiles
setup_bam_files.transcriptome = function(work_dir, variable_map = NULL){
  dt = setup_files(work_dir, "Aligned.toTranscriptome.out.bam$", variable_map)
  dt[, mapped_reads := get_mapped_reads(file), .(file)]
  dt[]
}

#' @export
#' @rdname tapfiles
setup_peak_files = function(work_dir, variable_map = NULL){
  setup_files(work_dir, ".narrowPeak$", variable_map)
}

#' @export
#' @rdname tapfiles
setup_count_files = function(work_dir, variable_map = NULL){
  setup_files(work_dir, ".ReadsPerGene.out.tab$", variable_map)
}

#' Creates data.table/data.frame with file paths and associated metadata.
#'
#' @param work_dir Path to TAP output directory
#' @param library_type Strandedness of library. Should be one of "guess", "unstranded", "first", or "second". Default is "guess" and will result in attempting to automatically detect and apply library type.
#' @param name_composition Used to extract column names for final count matrix from basename of count files. Should be a vector of numbers defining positions in _ and/or . delimited basename. For instance, name_composition = c(1, 3) on a basename of donor1_august_liver.siteA would yield a column name of donor1_liver.
#' @param just_check_library_type If TRUE, counts will not be loaded into a matrix and instead a vector of library type guesses will be returned. For a more detailed look at individual libraries use [guess_lib_from_file] with show_plots = TRUE.
#' @param gtf_file Supply to convert gene_ids to gene_names (HGNC gene symbols). Either a path to a gtf file or GRanges object loaded from one with [load_gene_reference] or [rtracklayer::import.gff].  The "gene_id" attribute will be matched to "gene_name" and aggregated by sum if required.
#' @param name_attribute The attribute in the gtf file to which gene_ids should be aggregated. Default of "gene_name" mostly works but some gtfs use "gene".
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
#' work_dir = "/slipstream_old/home/joeboyd/R_workspace/SFtapfly.data"
#' load_counts(work_dir)
#' load_counts(work_dir, name_composition = c(1, 2, 3, 4))
load_counts = function(work_dir, library_type = NULL, name_composition = NULL, just_check_library_type = FALSE, gtf_file = NULL, name_attribute = "gene_name"){
  use_gene_name = ifelse(is.null(gtf_file), FALSE, TRUE)
  if(is.null(name_composition)){
    setup_var_map = NULL
  }else{
    setup_var_map = name_composition
    var_names = paste0("var", seq_along(name_composition))
    names(setup_var_map) = var_names
  }
  if(is.data.frame(work_dir)){
    if(!is.null(name_composition)){
      stop("name_composition is not allowed when work_dir is a data.frame.")
    }
    if(!all(c("file", "name") %in% colnames(work_dir))){
      dt = work_dir
    }
  }else{
    dt = setup_count_files(work_dir, variable_map = setup_var_map)
  }
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
  if(is.null(library_type)){
    library_type = "guess"
  }
  if(library_type == "guess"){
    library_type = unique(sapply(dt$file, guess_lib_from_file))
    if(length(library_type) > 1){
      stop("Multiple potential library types detected: ", paste(library_type, collapse = ", "),
           "\nUse TAPhelpR::guess_lib_from_file to investigate and then run load_counts() again but manually specify library_type.")
    }
  }
  message("Loading all files using library type: ", library_type)
  mat = load_matrix_from_ReadsPerGene.out.tab(dt$file, library_type = library_type)
  colnames(mat) = dt$name
  if(use_gene_name){
    mat = aggregate_by_gene_name(mat, gtf_file = gtf_file, name_attribute = name_attribute)
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
load_norm_counts = function(work_dir, library_type = NULL, name_composition = NULL, just_check_library_type = FALSE, gtf_file = NULL){
  raw_counts = load_counts(
    work_dir = work_dir,
    library_type = library_type,
    name_composition = name_composition,
    just_check_library_type = just_check_library_type,
    gtf_file = gtf_file
  )
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
#' Loads gene reference from the same .gtf file used to run TAP. The result of this function can be input for `gtf_file` to any function with that parameter to avoid loading the same reference repeatedly.
#'
#' Only the gene entries are loaded by this function. Use [load_exon_reference] if you require exon information.
#'
#' @return A GRanges object that can be supplied as `gtf_file` for other TAPhelpR functions.
#' @export
#'
#' @examples
#' ref_dir = example_honeybee_reference()
#' gtf_file = file.path(ref_dir, "GTF/current.gtf")
#' load_gene_reference(gtf_file)
#' # for simplicity you can just specify the base path to the reference used to run TAP.
#' load_gene_reference(ref_dir)
load_gene_reference = function(gtf_file){
  .load_ref(gtf_file, "gene")
}

.load_ref = function(gtf_file, feature_type){
  if(dir.exists(gtf_file)){
    if(file.exists(file.path(gtf_file, "GTF/current.gtf"))){
      gtf_file = file.path(gtf_file, "GTF/current.gtf")
    }
  }
  if(!file.exists(gtf_file)){
    stop("gtf_file could not be found. Use the same reference directory used to run TAP or the gtf located at GTF/current.gtf in that same reference.")
  }
  ref_gr = rtracklayer::import.gff(gtf_file, feature.type = feature_type)
  names(ref_gr) = ref_gr$gene_id
  ref_gr
}

#' load_exon_reference
#'
#' Loads exon reference from the same .gtf file used to run TAP.
#'
#' This loads all exon entries  which is the most detail possibled. Use [load_gene_reference] if you only need basic gene information.
#'
#' @return A GRanges object with exonic information.
#' @export
#'
#' @examples
#' ref_dir = example_honeybee_reference()
#' gtf_file = file.path(ref_dir, "GTF/current.gtf")
#' load_exon_reference(gtf_file)
#' # for simplicity you can just specify the base path to the reference used to run TAP.
#' load_exon_reference(ref_dir)
load_exon_reference = function(gtf_file){
  .load_ref(gtf_file, "exon")
}
