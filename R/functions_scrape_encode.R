


#' ENCODE_get_file_info
#'
#' @param request_organisms Organisms to select for.
#' @param request_cell_lines Cell lines to select for.
#' @param request_target_chipseq ChIPseq target to select for.
#' @param request_file_types File types to select for.
#'
#' @return A data.table output from [ENCODExplorer::queryEncode()]. Filter and then pass to [ENCODE_download_files()]
#' @export
#'
#' @importFrom ENCODExplorer queryEncode
#' @examples
#' dt = ENCODE_get_file_info()
#' dt
ENCODE_get_file_info = function(request_organisms = NULL,
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
  dt[]
}

#' ENCODE_download_files
#'
#' @param dt Output of ENCODE_get_file_info with file_name attribute set. Inspect this data.table carefully and filter out unwanted files.
#' @param output_location Optional output directory. Will use current directory as default.
#'
#' @return The output location is invisibly returned.
#' @export
#'
#' @importFrom ENCODExplorer downloadEncode
#' @examples
#' dt = ENCODE_get_file_info()
#' #just to avoid downloading anything too big
#' dt = dt[grepl("Mb", file_size)]
#' dt = dt[sample(seq_len(nrow(dt)), 5), ]
#' dt$file_format
#' dt[, file_name := paste0("test_", seq(.N), ".", file_format)]
#' dt[, .(file_name)]
#' ENCODE_download_files(dt)
ENCODE_download_files = function(dt, output_location = getwd()){
  if(is.null(dt$file_name)){
    stop("You must set file_name attribute in input dt.")
  }
  if(is.null(dt$old_file)){
    dt[, old_file := paste0(file_accession, ".", file_format)]
  }
  dir.create(output_location, showWarnings = FALSE, recursive = TRUE)
  pbmcapply::pbmclapply(
    seq_len(nrow(dt)), mc.cores = 20,
    function(i){
      old_file = dt$old_file[i]
      new_file = dt$file_name[i]
      if(!file.exists(file.path(output_location, new_file))){
        ENCODExplorer::downloadEncode(dt$file_accession[i])
        file.rename(old_file, file.path(output_location, new_file))
      }
    }
  )
  invisible(output_location)
}

