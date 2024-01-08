SPLICE_EVENTS = list(
  "SkippingExon" = "SE",
  "Alternative5Prime" = "A5",
  "Alternative3Prime" = "A3",
  "MutuallyExclusiveExon" = "MX",
  "RetainedIntron" = "RI",
  "AlternativeFirstExon" = "AF",
  "AlternativeLastExon" = "AL",
  "Isoform" = 'isoform'
)

SPLICE_EVENTS.DECODE = unlist(SPLICE_EVENTS)
SPLICE_EVENTS.REVERSE = names(SPLICE_EVENTS)
names(SPLICE_EVENTS.REVERSE) = unlist(SPLICE_EVENTS)

find_suppa = function(){
  suppressWarnings({
    suppa_path = system("which suppa.py", intern = TRUE)
  })
  if(length(suppa_path) == 0){
    if(file.exists("/slipstream/home/joeboyd/anaconda2/envs/suppa2_env/bin/suppa.py")){
      suppa_path = "/slipstream/home/joeboyd/anaconda2/envs/suppa2_env/bin/suppa.py"
    }else{
      warning("Unable to locate suppa.py on current PATH. suppa_* functions will not work until set_suppa_path() is called with path to suppa.py.")
      suppa_path = "suppa.py"
    }
  }
  suppa_path
}

SUPPA_PATH = find_suppa()

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching TAPhelpR version ",
                        packageDescription("TAPhelpR")$Version, ".")
  TAP_SPLICE_EVENTS <<- SPLICE_EVENTS
  TAP_SPLICE_EVENTS.DECODE <<- SPLICE_EVENTS.DECODE
  TAP_SPLICE_EVENTS.REVERSE <<- SPLICE_EVENTS.REVERSE
  TAP_SUPPA_PATH <<-  SUPPA_PATH
}


#' set_suppa_path
#'
#' @param suppa_path Path to suppa.py. Will be used for all suppa_* functions.
#'
#' @return Invisibly returns TRUE.
#' @export
#'
#' @examples
#' set_suppa_path("/path/to/suppa.py")
set_suppa_path = function(suppa_path){
  TAP_SUPPA_PATH <<- suppa_path
  invisible(TRUE)
}

example_honeybee_config = function(){
  "/slipstream_old/home/joeboyd/R_workspace.combined/TAPhelpR.data/honeybee_TAP_input/config_tap.csv"
}

#' example_honeybee_input
#'
#' @return Path to example TAP output
#' @export
#'
#' @examples
#' example_honeybee_input()
example_honeybee_input = function(){
  "/slipstream_old/home/joeboyd/R_workspace.combined/TAPhelpR.data/honeybee_TAP_input"
}

#' example_honeybee_config
#'
#' @return Path to example TAP output
#' @export
#'
#' @examples
#' example_honeybee_config()
example_honeybee_config = function(){
  "/slipstream_old/home/joeboyd/R_workspace.combined/TAPhelpR.data/honeybee_TAP_input"
}


#' example_honeybee_output
#'
#' @return Path to example TAP output
#' @export
#'
#' @examples
#' example_honeybee_output()
example_honeybee_output = function(){
  "/slipstream_old/home/joeboyd/R_workspace.combined/TAPhelpR.data/honeybee_TAP_output.rename"
}

#' exampple_honeybee_metadata
#'
#' @return data.frame with honeybee SRR metadata.
#' @export
#'
#' @examples
#' exampple_honeybee_metadata()
example_honeybee_metadata = function(){
  f = system.file(package = "TAPhelpR", "extdata/honeybee_meta.csv", mustWork = TRUE)
  df = read.table(f, sep = ",")
  colnames(df) = c("srr", "name")
  df
}

#' example_honeybee_reference
#'
#' @return Path to example TAP reference
#' @export
#'
#' @examples
#' example_honeybee_reference()
example_honeybee_reference = function(){
  "/slipstream_old/home/joeboyd/R_workspace.combined/TAPhelpR.data/honeybee_TAP_reference"
}
