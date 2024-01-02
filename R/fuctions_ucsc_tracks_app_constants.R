
TEMPLATE = paste("track name=NAME description=DESCRIPTION",
                 "visibility=VISIBILITY autoScale=AUTOSCALE alwaysZero=ALWAYSZERO",
                 "maxHeightPixels=16:SIZE:128 viewLimits=MIN:MAX",
                 "color=COLOR bigDataUrl=URL type=TYPE",
                 "yLineOnOff=UCSC_SHOWYLINE yLineMark=UCSC_VALYLINE",
                 "gridDefault=UCSC_SHOWZERO smoothingWindow=UCSC_SMOOTHING",
                 "windowingFunction=UCSC_WINDOWFUN")
# FILE_ROOT = "/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks"
BASE_CFG_DIR = ".track_configs"
# CFG_DIR = paste0(FILE_ROOT, "/", CFG_DIR)
URL_ROOT = "https://galaxy.med.uvm.edu/static/UCSCtracks"
UCSC_VIS=c(hide = 0, dense = 1, full = 2, pack = 3, squish = 4)
UCSC_WIN=c("mean", "mean+whiskers", "maximum", "minimum")

