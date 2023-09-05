


get_encode_table = function(request_organisms = NULL,
                            request_cell_lines = NULL,
                            request_target_chipseq = NULL,
                            request_file_types = NULL){
  dt = ENCODExplorer::queryEncode(
    organism = request_organisms,
    biosample_name = request_cell_lines,
    target = request_target_chipseq,
    file_format = request_file_types,
    fixed = TRUE)
  message("pay attention to output_type! you probably don't want all of them.")
  dt
}

#' fetch_encode_files
#'
#' @param dt Output of get_encode_table with file_name attribute set.
#' @param out_dir Optional output directory. Will use current directory as default.
#'
#' @return
#' @export
#'
#' @examples
#' dt = get_encode_table()
#' #just to avoid downloading anything too big
#' dt = dt[grepl("Mb", file_size)]
#' dt = dt[sample(seq_len(nrow(dt)), 5), ]
#' dt$file_format
#' dt[, file_name := paste0("test_", seq(.N), ".", file_format)]
#' dt[, .(file_name)]
#' fetch_encode_files(dt)
fetch_encode_files = function(dt, out_dir = getwd()){
  if(is.null(dt$file_name)){
    stop("You must set file_name attribute in input dt.")
  }
  if(is.null(dt$old_file)){
    dt[, old_file := paste0(file_accession, ".", file_format)]
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  pbmcapply::pbmclapply(
    seq_len(nrow(dt)), mc.cores = 20,
    function(i){
      old_file = dt$old_file[i]
      new_file = dt$file_name[i]
      if(!file.exists(file.path(out_dir, new_file))){
        ENCODExplorer::downloadEncode(dt$file_accession[i])
        file.rename(old_file, file.path(out_dir, new_file))
      }
    }
  )
}

