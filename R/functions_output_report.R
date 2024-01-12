#' report_completion
#'
#' @param tap_out TAP output directory, either in progress or completed
#'
#' @return data.frame reporting completion status per sample.
#' @export
#' @rdname report-output
#' @examples
#' tap_out = example_honeybee_output()
#' report_completion(tap_out)
#' report_progress(tap_out)
#' report_progress_plot(tap_out)
#' report_errors(tap_out)
#'
#' tap_out = example_honeybee_output.in_progress()
#' report_completion(tap_out)
#' report_progress(tap_out)
#' report_progress_plot(tap_out)
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
#' @return Report additional details for each sample.
#' @export
#'
#' @rdname report-output
report_progress = function(tap_out){
  dt_jobs = .job_report(tap_out)
  if(any(dt_jobs$job_status == "error")){
    warning(paste0("Some samples did not complete successfully. Use report_errors for details on errors.\n",
                   "report_errors(", tap_out, ")"))
  }
  mat = dcast(dt_jobs, step~sample_name, value.var = "job_status", fill = NA)
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

.sample_report = function(tap_out){
  start_files = dir(tap_out, pattern = ".start$", full.names = TRUE)
  finish_files = dir(tap_out, pattern = ".finish$", full.names = TRUE)
  complete_files = dir(tap_out, pattern = ".complete$", full.names = TRUE)
  dt_pipeline_status = data.table(file = c(start_files, finish_files, complete_files))
  dt_pipeline_status$job_status = sapply(strsplit(basename(dt_pipeline_status$file), "\\."), function(x)x[length(x)])
  dt_pipeline_status[, sample_name := sub(paste0(".", job_status), "", basename(file)), .(file)]
  dt_pipeline_status$present = TRUE
  dt_pipeline_status$job_status = factor(dt_pipeline_status$job_status, levels = c("start", "finish", "complete"))
  dt_pipeline_status = dcast(dt_pipeline_status, sample_name~job_status, value.var = "present", fill = FALSE, drop = FALSE)
  dt_pipeline_status$step = "sample_status"
  dt_pipeline_status$step = factor(dt_pipeline_status$step, levels = step_lev)
  dt_pipeline_status$sample_status = "none"
  dt_pipeline_status[start == TRUE, sample_status := "in progress"]
  dt_pipeline_status[start == TRUE & finish == TRUE, sample_status := "incomplete, errors"]
  dt_pipeline_status[start == TRUE & finish == TRUE & complete == TRUE, sample_status := "completed successfully"]
  dt_pipeline_status[]
}

# step_lev = c("STAR_align", "bsortindex", "salmon_quant", "suppa2", "exactSNP", "make_bigwigs", "finish", "complete", "sample_status")
step_lev = c("STAR_align", "bsortindex", "salmon_quant", "suppa2", "exactSNP", "make_bigwigs", "sample_status")
# step_lev.main = setdiff(step_lev, c("finish", "complete"))

.job_report = function(tap_out){
  start_files = dir(tap_out, pattern = ".start$", full.names = TRUE)
  finish_files = dir(tap_out, pattern = ".finish$", full.names = TRUE)
  complete_files = dir(tap_out, pattern = ".complete$", full.names = TRUE)
  dt_pipeline_status = data.table(file = c(start_files, finish_files, complete_files))
  dt_pipeline_status$job_status = sapply(strsplit(basename(dt_pipeline_status$file), "\\."), function(x)x[length(x)])
  dt_pipeline_status[, sample_name := sub(paste0(".", job_status), "", basename(file)), .(file)]
  dt_pipeline_status$present = TRUE
  dt_pipeline_status$job_status = factor(dt_pipeline_status$job_status, levels = c("start", "finish", "complete"))
  sample_lev = unique(dt_pipeline_status$sample_name)
  dt_pipeline_status = dcast(dt_pipeline_status, sample_name~job_status, value.var = "present", fill = FALSE, drop = FALSE)



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
  dt_sample_status = data.table(log_file = out_files, finished = is_finished)
  dt_sample_status[, step := sub("\\..+", "", basename(log_file))]
  dt_sample_status[, sample_name := sub(".logs", "", basename(dirname(log_file)))]
  dt_sample_status[, job_status := ifelse(finished, "complete", "in_progress")]

  dt_logs = dt_sample_status[, .(log_file, step, sample_name)]

  dt_sample_status = dcast(dt_sample_status, sample_name~step, value.var = "job_status", fill = "not_started")

  dt_jobs = merge(dt_sample_status, dt_pipeline_status, by = "sample_name")
  dt_jobs = melt(dt_jobs, id.vars = c("sample_name", "complete", "finish", 'start'), variable.name = "step", value.name = "job_status")
  dt_jobs[finish == TRUE & complete == FALSE & job_status == "in_progress", job_status := "error"]

  dt_jobs = merge(dt_jobs, dt_logs, by = c("step", "sample_name"), all.x = TRUE)

  dt_jobs$step = factor(dt_jobs$step, levels = step_lev)
  dt_jobs$sample_name = factor(dt_jobs$sample_name, levels = sample_lev)

  dt_jobs[]
}

.sample_report = TAPhelpR:::.sample_report
.job_report = TAPhelpR:::.job_report

#' report_progress_plot
#'
#' @return Plot of job progress.
#' @export
#' @rdname report-output
report_progress_plot = function(tap_out){
  dt_jobs = .job_report(tap_out)
  dt_samples = .sample_report(tap_out)
  dt_jobs$step = factor(dt_jobs$step, levels = rev(levels(dt_jobs$step)))
  dt_samples$step = factor(dt_samples$step, levels = rev(levels(dt_samples$step)))

  dt_jobs$job_status = factor(dt_jobs$job_status, levels = c("not_started", "error", "in_progress", "complete"))
  dt_samples$sample_status = factor(dt_samples$sample_status, levels = c("incomplete, errors", "in progress", "completed successfully"))

  p_size = 8
  p = ggplot(dt_jobs, aes(x = sample_name, y = step)) +
    geom_tile(aes(fill = job_status), color = "white", linewidth = 2) +
    scale_y_discrete(drop = FALSE) +
    geom_point(data = dt_samples, color = "black", size = 1.2*p_size, shape = 19) +
    geom_point(data = dt_samples, aes(color = sample_status), size = p_size, shape = 19) +
    scale_fill_manual(values = c("complete" = "forestgreen", "error" = "darkred", "in_progress" = "yellow2", "not_started" = "gray"), drop = FALSE) +
    scale_color_manual(values = c("in progress" = "yellow2", "incomplete, errors" = "darkred", "completed successfully" = "forestgreen"), drop = FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.background = element_blank(),
          panel.grid = element_blank()) +
    labs(x = "", y = "")
  p
}


#' report_errors
#'
#' @return Locate jobs that experienced errors and return relevant log files.
#' @export
#' @rdname report-output
report_errors = function(tap_out){
  dt_jobs = .job_report(tap_out)
  dt_error = dt_jobs[job_status == "error"]
  if(nrow(dt_error) == 0){
    message("No errors were detected.")
  }else{
    table(dt_error$step)
    dt_error = dt_error[order(sample_name)][order(step)]
    dt_error.report = dt_error[, .(sample_name, step)]
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
