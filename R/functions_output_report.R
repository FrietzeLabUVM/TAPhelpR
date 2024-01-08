tap_out = "~/R_workspace.combined/TAPhelpR.data/honeybee_TAP_output.rename"

#' report_completion
#'
#' @param tap_out TAP output directory, either in progress or completed
#'
#' @return
#' @export
#' @rdname report-output
#' @examples
#' tap_out = "~/R_workspace.combined/TAPhelpR.data/honeybee_TAP_output.rename"
#' report_completion(tap_out)
#' report_progress(tap_out)
#' report_errors(tap_out)
report_completion = function(tap_out){
  detect_files = function(suffix, df = NULL){
    match_files = dir(tap_out, pattern = paste0(suffix, "$"))
    suffix_name = sub("^\\.", "", suffix)
    if(is.null(df)){
      df = data.frame(name = sub(suffix, "", basename(match_files)))
      rownames(df) = df$name
    }
    df[[suffix_name]] = FALSE
    if(length(match_files) > 0){
      df[sub(suffix, "", match_files),][[suffix_name]] = TRUE
    }
    df
  }

  stat_df = detect_files(".start")
  stat_df = detect_files(".finish", stat_df)
  stat_df = detect_files(".complete", stat_df)

  message(nrow(stat_df), " samples detected")
  finish_perc = paste0(round(sum(stat_df$finish)/nrow(stat_df)*100, 2), "%")
  message(sum(stat_df$finish), " (", finish_perc, ") are no longer in process.")

  complete_perc = paste0(round(sum(stat_df$complete)/nrow(stat_df)*100, 2), "%")
  message(sum(stat_df$complete), " (", complete_perc, ") have completed successfully.")


  n_errors = sum(stat_df$finish) - sum(stat_df$complete)
  if(n_errors > 0){
    error_perc = paste0(round(n_errors / nrow(stat_df)*100, 2), "%")
    message(n_errors, " (", error_perc, ") have errors that should be investigated.")
  }

  stat_df
}


#' report_progress
#'
#' @param tap_out
#'
#' @return
#' @export
#'
#' @rdname report-output
report_progress = function(tap_out){
  log_dirs = dir(tap_out, pattern = "logs$", full.names = TRUE)
  log_files = dir(log_dirs, full.names = TRUE)
  echo_files = log_files[grepl("echo_sub", log_files)]
  echo_files = echo_files[grepl("out$", echo_files)]

  log_files = log_files[!grepl("echo_sub", log_files)]
  log_files = log_files[!grepl("finish", log_files)]
  log_files = log_files[!grepl("completion", log_files)]

  out_files = log_files[grepl("out$", log_files)]
  out_files = filter_files_for_most_recent(out_files)


  is_finished = sapply(out_files, function(x){
    df = read.table(x, sep = "\n")
    df[nrow(df),] == "FINISHED"
  })
  dt1 = data.table(log_file = out_files, finished = is_finished)
  dt1[, step := sub("\\..+", "", basename(log_file))]

  finish_files = sub(".logs$", ".finish", dirname(echo_files))
  complete_files = sub(".logs$", ".finish", dirname(echo_files))
  dt2 = rbind(
    data.table(log_file = echo_files, finished = file.exists(finish_files), step = "finish"),
    data.table(log_file = echo_files, finished = file.exists(complete_files), step = "complete")
  )

  dt = rbind(dt1, dt2)

  dt[, name := sub(".logs", "", basename(dirname(log_file)))]


  step_lev = c("STAR_align", "bsortindex", "salmon_quant", "suppa2", "exactSNP", "make_bigwigs", "finish", "complete")
  stopifnot(all(dt$step %in% step_lev))
  dt$step = factor(dt$step, levels = step_lev)

  finished_samples = dt[step == "finish"][finished == TRUE]$name
  complete_samples = dt[step == "complete"][finished == TRUE]$name
  finished_not_completed = setdiff(finished_samples, complete_samples)
  if(length(finished_not_completed > 0)){
    warning(paste0("Some samples did not complete successfully. Use report_errors for details on errors.\n",
                   "report_errors(", tap_out, ")"))
  }

  dt$step = factor(dt$step, levels = rev(step_lev))
  ggplot(dt, aes(x = name, y = step, fill = finished)) +
    geom_tile(color = "white", linewidth = 1.4) +
    scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "black")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), panel.background = element_blank(), panel.grid = element_blank())

  dt$step = factor(dt$step, levels = step_lev)
  mat = dcast(dt, step~name, value.var = "finished", fill = NA)
  mat
}

filter_files_for_most_recent = function(f){
  f_info = file.info(f)
  f_info = f_info[order(f_info$mtime, decreasing = TRUE),]
  f_info$sample_name = sub(".logs$", "", basename(dirname(rownames(f_info))))
  f_info$job_names = sapply(strsplit(basename(rownames(f_info)), "\\."), function(x){
    paste(x[-(length(x)-1):-length(x)], collapse = ".")
  })
  is_dupe = duplicated(paste(f_info$sample_name, f_info$job_names))
  # table(f_info$job_names[!is_dupe])
  # table(f_info$job_names[is_dupe])
  rownames(f_info)[!is_dupe]
}

#' report_errors
#'
#' @param tap_out
#'
#' @return
#' @export
#' @rdname report-output
report_errors = function(tap_out){
  log_dirs = dir(tap_out, pattern = "logs$", full.names = TRUE)
  log_files = dir(log_dirs, full.names = TRUE)

  echo_files = log_files[grepl("echo_sub", log_files)]
  echo_files = echo_files[grepl("out$", echo_files)]

  log_files = log_files[!grepl("echo_sub", log_files)]
  log_files = log_files[!grepl("finish", log_files)]
  log_files = log_files[!grepl("completion", log_files)]

  out_files = log_files[grepl("out$", log_files)]
  out_files = filter_files_for_most_recent(out_files)

  is_finished = sapply(out_files, function(x){
    df = read.table(x, sep = "\n")
    df[nrow(df),] == "FINISHED"
  })
  dt = data.table(out_files = out_files, finished = is_finished)
  dt[, error_files := sub(".out$", ".error", out_files)]
  dt[, sample_name := sub(".logs$", "", basename(dirname(out_files)))]
  dt$job_name = sapply(strsplit(basename(dt$out_files), "\\."), function(x){
    paste(x[-(length(x)-1):-length(x)], collapse = ".")
  })

  finish_files = sub(".logs$", ".finish", dirname(echo_files))
  dt_finished = data.table(finish_files = finish_files)
  dt_finished[, sample_name := sub(".finish", "", basename(finish_files))]

  #errors are present when job output is not FINISHED but sample is finished

  dt_finished = merge(dt, dt_finished, by = "sample_name")

  dt_error = dt_finished[finished == FALSE]
  if(nrow(dt_error) == 0){
    message("No errors were detected.")
  }else{
    job_lev = c("STAR_align", "bsortindex", "salmon_quant", "suppa2", "exactSNP", "make_bigwigs")
    dt_error$job_name = factor(dt_error$job_name, levels = job_lev)
    table(dt_error$job_name)
    dt_error = dt_error[order(sample_name)][order(job_name)]
    dt_error.report = dt_error[, .(sample_name, job_name)]
    message("Errors were detected at the following steps:")
    print(dt_error.report)
    message(paste(collapse = "\n", c(
      "See the following log files to diagnose these errors:",
      dt_error$error_files
    )))
    invisible(dt_error)
  }
  invisible(NULL)

}
