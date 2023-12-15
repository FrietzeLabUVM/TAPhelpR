#
# for(name in names(bam_df.grouped)){
#   message(name)
#   ssvSplicing::suppa_joinFiles(bam_df.grouped[[name]], output_name = name)
# }
#
# suppa_diffSplice.within_group
#
# supp_diffSplice.within_group = function(run_by){
#   bam_df.sp = split(bam_df, bam_df[[run_by]])
#   sub_wd = file.path(diff_out, run_by)
#   grp = names(bam_df.sp)[1]
#   for(grp in names(bam_df.sp)){
#     stage_dir = file.path(sub_wd, grp)
#     dir.create(stage_dir, showWarnings = FALSE, recursive = TRUE)
#     sel_bam_df = bam_df.sp[[grp]]
#     lapply(unique(sel_bam_df$group), function(nam){
#       match_files = dir(diff_out, pattern = paste0(nam, "\\."), full.names = TRUE)
#       file.symlink(normalizePath(match_files), stage_dir)
#       match_files
#     })
#     ssvSplicing::suppa_diffSplice(
#       wd = stage_dir,
#       ref_location = "/slipstream/home/joeboyd/indexes/honeybee/SUPPA2"
#     )
#   }
# }


library(TAPhelpR)
tap_out = "~/R_workspace.combined/TAPhelpR.data/honeybee_TAP_output"
bam_files = setup_bam_files(tap_out, var_map = c("day", "sex", "rep"))

suppa_joinFiles(bam_files, by = "sex")
suppa_diffSplice(
  ref_location = "~/lab_shared/indexes/honeybee",
  PSI_todo = TAP_SPLICE_EVENTS$SkippingExon)
debug(suppa_diffSplice.within_group)
suppa_diffSplice.within_group(
  input_files = bam_files,
  within_group = "day",
  between_group = "sex",
  ref_location = "~/lab_shared/indexes/honeybee",
  PSI_todo = TAP_SPLICE_EVENTS$SkippingExon)
suppa_clusterEvents(PSI_todo = TAP_SPLICE_EVENTS$SkippingExon)

