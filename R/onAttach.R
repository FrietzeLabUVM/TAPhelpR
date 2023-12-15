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

