#' launch_UCSC_tracks_app
#'
#' @param track_hosting_dir The directory where track files are available via a
#'   public URL. To conveniently place pipeline outputs in a public location,
#'   use the [stage_output_for_UCSC_tracks] function.
#' @param base_host_path The portion of the `track_hosting_dir` to replace with
#'   `base_url` to generate a valid URL.
#' @param base_url The public URL pointing to `base_host_path` locally.
#'
#' @return A shiny application is launched. Nothing is returned.
#' @export
#'
#' @examples
#' work_dir = "/slipstream/home/joeboyd/public_files/honeybee"
#' track_hosting_dir = "/slipstream/home/joeboyd/public_files/honeybee"
#' base_host_path = "/slipstream/galaxy/production/galaxy-dist"
#' if(!grepl(base_host_path, track_hosting_dir)){
#'   sp = strsplit(track_hosting_dir, "/")[[1]]
#'   user = sp[which(sp == "home") + 1]
#'   track_hosting_dir = paste0("/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/", user, "_files")
#'   if(!dir.exists(track_hosting_dir)){
#'     "track_hosting_dir is not a sub-directory of base_host_path and could not convert to a UVM galaxy path."
#'   }
#' }
#' launch_UCSC_tracks_app(track_hosting_dir = track_hosting_dir,
#'                        base_host_path,
#'                        base_url = "https://galaxy.med.uvm.edu")
launch_UCSC_tracks_app = function(track_hosting_dir, base_host_path, base_url){
  if(!grepl(base_host_path, track_hosting_dir)){
    stop("track_hosting_dir must be a sub-directory (possibly indirectly) of base_host_path.")
  }
  shiny::runApp(list(ui = ucsc_track_ui(), server = ucsc_track_server(track_hosting_dir, base_host_path, base_url)))
  invisible()
}

# launch_UCSC_tracks_app(track_hosting_dir = "~/public_files/honeybee")
#"/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/joeboyd_files/honeybee"


#' launch_UCSC_tracks_app.UVM_galaxy
#'
#' This variant of launch_UCSC_tracks_app has convenient presets for users of the Galaxy server at UVM and should not be used by others.
#'
#' @param track_hosting_dir The directory where track files are available via a public URL. To conveniently place pipeline outputs in a public location, use the [stage_output_for_UCSC_tracks] function.
#'
#' @return A shiny application is launched. Nothing is returned.
#' @export
#'
#' @examples
#' launch_UCSC_tracks_app.UVM_galaxy("/slipstream/home/joeboyd/public_files/honeybee")
launch_UCSC_tracks_app.UVM_galaxy = function(track_hosting_dir){
  base_host_path = "/slipstream/galaxy/production/galaxy-dist"
  if(!grepl(base_host_path, track_hosting_dir)){
    track_hosting_dir = normalizePath(track_hosting_dir)
    if(!grepl(base_host_path, track_hosting_dir)){
      sp = strsplit(track_hosting_dir, "/")[[1]]
      user = sp[which(sp == "home") + 1]
      galaxy_track_dir = "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks"
      local_path = sub(".+public_files/?", "", track_hosting_dir)
      track_hosting_dir = file.path(galaxy_track_dir, paste0(user, "_files"), local_path)
      if(!dir.exists(track_hosting_dir)){
        stop("track_hosting_dir is not a sub-directory of base_host_path and could not convert to a UVM galaxy path.")
      }
    }
  }
  launch_UCSC_tracks_app(track_hosting_dir = track_hosting_dir,
                         base_host_path,
                         base_url = "https://galaxy.med.uvm.edu")
}



#' stage_output_for_UCSC_tracks
#'
#' Uses symbolic links to make bigwig files available at a publicly hosted
#' location for later configuration for viewing on UCSC. Uses
#' [launch_UCSC_tracks_app] with the same `track_hosting_dir` to configure a
#' UCSC browser session.
#'
#' @param pipeline_outputs_dir The directory where completed TAP outputs are
#'   located.
#' @param track_hosting_dir The directory where track files are available via a
#'   public URL. To conveniently place pipeline outputs in a public location,
#'   use the [stage_output_for_UCSC_tracks] function.
#'
#' @return `track_hosting_dir` is invisibly returned.
#' @export
#'
#' @examples
#' stage_output_for_UCSC_tracks("~/R_workspace.combined/TAPhelpR.data/honeybee_TAP_output", "~/public_files/honeybee_tracks")
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
  invisible(track_hosting_dir)
}



