

#' stage_output_for_UCSC_tracks
#'
#' @param pipeline_outputs_dir
#'
#' @return
#' @export
#'
#' @examples
#' config_file = "/slipstream/home/joeboyd/lab_shared/data/SF_SE_data"
#' config_file = "/slipstream/home/joeboyd/R_workspace.combined/chiptsne_projects/zhang_IK_SE/config_zhang_IK_SE.csv"
#' pipeline_outputs_dir = "/slipstream/home/joeboyd/lab_shared/data/SF_SE_data"
#' pipeline_outputs_dir = "/slipstream/home/joeboyd/R_workspace.combined/chiptsne_projects/zhang_IK_SE/output"
#'
#'
#' c("fe", "nu", "ns", "ru", "rs")
#'
#' stage_output_for_UCSC_tracks("/slipstream/home/joeboyd/R_workspace.combined/chiptsne_projects/zhang_IK_SE/output")
#' stage_output_for_UCSC_tracks("/slipstream/home/joeboyd/lab_shared/data/SF_SE_data")
stage_output_for_UCSC_tracks = function(pipeline_outputs_dir, track_hosting_dir){
  bw_outs = dir(pipeline_outputs_dir, pattern = "bigwigs$", full.names = TRUE)
  bw_files.fe = dir(pipeline_outputs_dir, pattern = "FE.bw$", full.names = TRUE)
  bw_files.reads = dir(bw_outs, pattern = ".bw$", full.names = TRUE)

  derive_bw_group = function(f){
    bw_suffixes = sub(".+.Aligned.sortedByCoord.out_", "", basename(f))
    norms = sapply(strsplit(bw_suffixes, "[\\._]"), function(x)x[1])
    strands = sapply(strsplit(bw_suffixes, "[\\._]"), function(x)x[2])
    show_splice = sapply(strsplit(bw_suffixes, "[\\._]"), function(x)x[3])
    show_splice[show_splice == "bw"] = ""
    grp = paste(norms, strands, show_splice, sep = "_")
    sub("_$", "", grp)
  }

  bw_suffixes = sub(".+.Aligned.sortedByCoord.out_", "", basename(bw_files.reads))
  norms = sapply(strsplit(bw_suffixes, "[\\._]"), function(x)x[1])
  strands = sapply(strsplit(bw_suffixes, "[\\._]"), function(x)x[2])
  show_splice = sapply(strsplit(bw_suffixes, "[\\._]"), function(x)x[3])
  show_splice[show_splice == "bw"] = "noSplice"
  valid_groups = unique(paste(norms, strands, show_splice))
  if(length(bw_files.fe) > 0){
    valid_groups = c(valid_groups, "fold-enrichment")
  }
  valid_groups = gsub(" ", "_", valid_groups)
  valid_short = sapply(strsplit(valid_groups, "[ -_]"), function(x){
    paste(substr(x, 1, 1), collapse = "")
  })
  # message(paste(valid_groups, collapse = " "))
  # message(paste(valid_short, collapse = " "))
  bw_df = data.frame(file = c(bw_files.fe, bw_files.reads))
  bw_df$group = derive_bw_group(bw_df$file)
  for(i in seq_len(nrow(bw_df))){
    f = bw_df$file[i]
    group = bw_df$group[i]
    dest_dir = file.path(track_hosting_dir, group)
    dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
    if(!file.exists(file.path(dest_dir, basename(f)))){
      file.symlink(f, file.path(dest_dir, basename(f)))
    }
  }
}


stage_rnaseq_for_UCSC = function(){

}
