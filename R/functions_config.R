config_create = function(){

}

get_header_args = function(){
  arg_parse_str = '
-c|--config) echo ignoring config file specified in config file.; shift ;;
-f1s|--f1_suffix) F1_suff="$2"; shift ;;
-f2s|--f2_suffix) F2_suff="$2"; shift ;;
-i|--inDir) input="$2"; shift ;;
-o|--outDir) align_path="$2"; shift ;;
-j|--jobDir|--jobsDir) JOBS_PATH="$2"; shift ;;
-ref|--reference) ref="$2"; shift ;;
-idx|--starIndex) star_index="$2"; shift ;;
-s|--suppaRef) suppa_ref="$2"; shift ;;
-g|--gtf) gtf="$2"; shift ;;
-fa|--fasta) fasta="$2"; shift ;;
-rDNA|--rDNA_starIndex) rDNA_index="$2"; shift ;;
-PE|--PE) read_mode=PE ;;
-SE|--SE) read_mode=SE ;;
-noSub|--noSub) sub_mode=bash ;;
-p|--pipeline) pipeline="$2"; shift ;;
-sl|--scriptLocation) scripts="$2"; shift ;;
-docker|--docker) docker="$2"; shift ;;
-singularity|--singularity) singularity="$2"; shift ;;
*) echo "Unknown parameter passed: $1"; cat $SCRIPT_PATH/help_msg.txt; exit 1 ;;
'
  arg_parse_str = strsplit(arg_parse_str, "\n")[[1]]
  arg_parse_str = arg_parse_str[!grepl("Unknown parameter passed", arg_parse_str)]
  arg_parse_str = arg_parse_str[nchar(arg_parse_str) > 0]
  value_args = arg_parse_str[grepl("shift", arg_parse_str)]
  flag_args = arg_parse_str[!grepl("shift", arg_parse_str)]

  parse_arg_keys = function(args){
    args = sapply(
      strsplit(args, ")"),
      function(x)x[1]
    )
    args =strsplit(args, "\\|")
    names(args) = gsub("-", "", sapply(args, function(x)x[2]))
    args
  }

  value_args = parse_arg_keys(value_args)
  flag_args = parse_arg_keys(flag_args)
  list(value_args = value_args, flag_args = flag_args)
}

parse_head_str = function(head_str){
  head_args = strsplit(head_str, " +")[[1]]
  header_args = get_header_args()
  parsed_args = list()
  for(v_arg in names(header_args$value_args)){
    # message(v_arg)
    possible_keys = header_args$value_args[[v_arg]]
    for(v_key in possible_keys){
      arg_match = grepl(v_key, head_args)
      if(any(arg_match)){
        if(is.null(parsed_args[[v_arg]])){
          # value follows key by one position
          parsed_args[[v_arg]] = head_args[which(arg_match)+1]
        }else{
          #should not be present already
          stop("Duplicate key present for ", v_arg, " argument: redundant key is ", v_key)
        }

      }
    }
  }
  for(f_arg in names(header_args$flag_args)){
    # message(f_arg)
    possible_keys = header_args$flag_args[[f_arg]]
    for(v_key in possible_keys){
      arg_match = grepl(v_key, head_args)
      if(any(arg_match)){
        if(is.null(parsed_args[[f_arg]])){
          # value follows key by one position
          parsed_args[[f_arg]] = TRUE
        }else{
          #should not be present already
          stop("Duplicate key present for ", f_arg, " argument: redundant key is ", v_key)
        }

      }
    }
  }
  parsed_args
}

parse_header = function(config_file){
  head_df = read.table(config_file, sep = "\n", header = FALSE, comment.char = "")
  head_df = subset(head_df, grepl("^#CFG", V1))
  head_str = gsub("#CFG", "", paste(head_df$V1, collapse = " "))
  parse_head_str(head_str)
}

.config_validate = function(config_file, extra_args){
  body_df = read.table(config_file, comment.char = "#", sep = ",")
  cfg_params = parse_header(config_file)
  if(!is.null(cfg_params$config)){
    warning("--config should not be specified in your config file but this will not prevent TAP from running.")
  }
  if(!is.null(extra_args)){
    override_params = parse_head_str(extra_args)
    for(nam in names(override_params)){
      cfg_params[[nam]] = override_params[[nam]]
    }
  }

  must_exist = c("inDir", 'reference', "jobDir", "reference", "starIndex", "suppaRef", "gtf", "fasta", "rDNA_starIndex", "pipeline", "scriptLocation", "singularity")
  for(param in must_exist){
    if(!is.null(cfg_params[[param]])){
      if(!file.exists(cfg_params[[param]])){
        stop(param, ' : ', cfg_params[[param]], " must exist but was not found.")
      }
    }
  }
  #apply defaults
  if(is.null(cfg_params$inDir)){
    cfg_params$inDir= paste0(getwd(), " (default)")
  }
  if(is.null(cfg_params$outDir)){
    cfg_params$outDir= paste0(getwd(), "/TAP_output (default)")
  }
  if(is.null(cfg_params$PE) & is.null(cfg_params$SE)){
    #PE is default for RNAseq
    cfg_params$PE = "default"
  }
  if(is.null(cfg_params$f1_suffix)){
    cfg_params$f1_suffix = "_R1_001.fastq.gz (default)"
  }
  if(is.null(cfg_params$f2_suffix)){
    cfg_params$f2_suffix = "_R2_001.fastq.gz (default)"
  }
  cfg_params
  #check body (files exist, names are well formed)
  in_dir = sub(" \\(default\\)", "", cfg_params$inDir)
  f1_suffix = sub(" \\(default\\)", "", cfg_params$f1_suffix)
  f2_suffix = sub(" \\(default\\)", "", cfg_params$f2_suffix)
  f1_files = body_df[[1]]
  f1_files = unlist(strsplit(f1_files, "&"))
  fq_1_exists = file.exists(file.path(in_dir, f1_files))
  if(!all(fq_1_exists)){
    stop(
      paste(c("Some R1 fastqs were not found!:",
              f1_files[!fq_1_exists]), collapse = "\n")
    )
  }
  fq_1_has_suffix = grepl(f1_suffix, f1_files)
  if(!all(fq_1_has_suffix)){
    stop(
      paste(c("Some R1 fastqs missing R1 suffix!:",
              paste0("suffix:", f1_suffix),
              f1_files[!fq_1_has_suffix]), collapse = "\n")
    )
  }
  if(!is.null(cfg_params$PE)){
    f2_files = sub(f1_suffix, f2_suffix, f1_files)
    fq_2_exists = file.exists(file.path(in_dir, f2_files))
    if(!all(fq_2_exists)){
      stop(
        paste(c("Some R2 fastqs were not found!:",
                f2_files[!fq_2_exists]), collapse = "\n")
      )
    }
  }

  if(ncol(body_df) == 1){
    body_df$name = sub(f1_suffix, "", basename(f1_files))
  }
  colnames(body_df)[2] = "name"
  names.sp = strsplit(body_df$name, "[_\\.]")
  max_len = max(lengths(names.sp))
  if(!all(lengths(names.sp) == max_len)){
    warning("Some names are not well formed. They do not have the same number of _ and . delimitted elements.\nTAP can run but this will complicated some downstream operations.")
  }
  names.sp = lapply(names.sp, function(x){
    if(length(x) < max_len){
      x = c(x, rep(">EMPTY<", max_len - length(x)))
    }
    x
  })
  df_names = as.data.frame(matrix(unlist(names.sp), byrow = TRUE, ncol = max_len))
  rownames(df_names) = body_df[[1]]
  df_names
}


#' Validate config file for running TAP on RNAseq
#'
#' @param config_file A configuration file for TAP.
#'
#' @return A data.frame with the metadata elements parsed for each found fastq file.
#' @export
#'
#' @examples
#' cfgs = dir("~/lab_shared/scripts/TAP/testing/test_configs/", full.names = TRUE)
#' cfgs = cfgs[!grepl("chip", cfgs)]
#' names(cfgs) = basename(cfgs)
#' cfgs = as.list(cfgs)
#'
#' test_dir = "/slipstream/home/joeboyd/lab_shared/scripts/TAP/testing"
#' cmd_extra = paste("-i", file.path(test_dir, "test_data/fastq_rnaseq_PE"), "-ref", file.path(test_dir, "references/dm6"))
#' config_validate(cfgs$test_dm6_config.basic.csv, extra_args = cmd_extra)
#'
#' config_validate(cfgs$test_dm6_config.basic.csv, wd = "~/lab_shared/scripts/TAP/testing/test_data/fastq_rnaseq_PE/")
#'
#' config_validate(cfgs$test_dm6_config.params.csv, wd = "~/lab_shared/scripts/TAP")
#'
#' config_validate(cfgs$test_dm6_config.pool.csv, wd = "~/lab_shared/scripts/TAP/testing/test_data/fastq_rnaseq_PE/")
#'
#' config_validate(cfgs$test_dm6_config.rDNA_only.csv, wd = "~/lab_shared/scripts/TAP")
#'
#' config_validate(cfgs$test_dm6_config.rename.csv, wd = "~/lab_shared/scripts/TAP/testing/test_data/fastq_rnaseq_PE/")
#'
#' config_validate(cfgs$test_dm6_config.SE.csv, wd = "~/lab_shared/scripts/TAP")
#'
config_validate = function(config_file, extra_args = NULL, wd = NULL){
  old_wd = getwd()
  if(!is.null(wd)){
    setwd(wd)
  }
  #tryCatch so that wd gets reset regardless of error
  df_names = tryCatch(
    expr = {
      .config_validate(config_file, extra_args)
    },
    error = function(e){
      e
    },
    finally = {
      setwd(old_wd)
    })
  if(is(df_names, "error")){
    stop(df_names)
  }
  df_names
}
