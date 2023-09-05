config_file = "/slipstream/home/joeboyd/lab_shared/data/SF_SE_data"
config_file = "/slipstream/home/joeboyd/R_workspace.combined/chiptsne_projects/zhang_IK_SE/config_zhang_IK_SE.csv"
pipeline_outputs_dir = "/slipstream/home/joeboyd/lab_shared/data/SF_SE_data"
pipeline_outputs_dir = "/slipstream/home/joeboyd/R_workspace.combined/chiptsne_projects/zhang_IK_SE/output"


c("fe", "nu", "ns", "ru", "rs")

stage_output_for_UCSC_tracks = function(pipeline_outputs_dir){
  bw_outs = dir(pipeline_outputs_dir, pattern = "bigwigs$", full.names = TRUE)
  bw_files.fe = dir(pipeline_outputs_dir, pattern = "FE.bw$", full.names = TRUE)
  bw_files.reads = dir(bw_outs, pattern = ".bw$", full.names = TRUE)
  norms = sapply(strsplit(basename(bw_files.reads), "[\\._]"), function(x)x[length(x)-2])
  strands = sapply(strsplit(basename(bw_files.reads), "[\\._]"), function(x)x[length(x)-1])
  valid_groups = unique(paste(norms, strands))
  if(length(bw_files.fe) > 0){
    valid_groups = c(valid_groups, "fold-enrichment")
  }
  valid_short = sapply(strsplit(valid_groups, "[ -]"), function(x){
    paste(substr(x, 1, 1), collapse = "")
  })
  message(paste(valid_groups, collapse = " "))
  message(paste(valid_short, collapse = " "))
}
stage_output_for_UCSC_tracks("/slipstream/home/joeboyd/R_workspace.combined/chiptsne_projects/zhang_IK_SE/output")
stage_output_for_UCSC_tracks("/slipstream/home/joeboyd/lab_shared/data/SF_SE_data")

stage_rnaseq_for_UCSC = function(){

}
