#' config_create
#'
#' @param inDir Directory with all .fastq.gz files. Default is current directory.
#' @param f1_suffix Suffix of R1 files. Default is _R1_001.fastq.gz
#' @param f2_suffix Suffix of R2 files. Default is _R2_001.fastq.gz. Only required if data is PE.
#' @param outDir Output directory. Default is TAP_output in current directory.
#' @param reference Parent directory for all reference directories. The location with output from setup_new_reference.sh
#' @param starIndex Directory with STAR index files. Only required if `reference` is not provided or you wish to override.
#' @param suppaRef Directory with suppa2 index files. Only required if `reference` is not provided or you wish to override.
#' @param gtf File, gtf or gff, with gene annotation. Only required if `reference` is not provided or you wish to override.
#' @param fasta Organism genomic sequence fasta. Only required if `reference` is not provided or you wish to override.
#' @param rDNA_starIndex Directory with STAR index for rDNA. Only required if you wish to quantify rDNA read composition.
#' @param PE Set TRUE if data is PE. RNAseq will default to PE.
#' @param SE Set TRUE if data is SE. ChIPseq will default to SE.
#' @param noSub Set TRUE if bash should be used for all commands in serial instead of sbatch.
#' @param pipeline Pipeline script to use. Only required if submission script is not in the same directory as default of "rnaseq_pipeline.sh".
#' @param jobDir Directory with all job scripts. Only required if submission script is not in the same directory as default of "deployed_job_scripts".
#' @param scriptLocatio Deprecated. Don't use.
#' @param docker If you want to use docker, docker image to use. i.e. jrboyd/tap after using "setup_scripts/pull_docker_image.sh"
#' @param singularity If you want to use singularity, singularity .sif file to use. i.e. tap_latest.sif after using "setup_scripts/pull_singularity_image.sh"
#'
#' @return list of config lines and fastq lines
#' @export
#' @rdname config-write
#'
#' @examples
#' cfg = config_create(
#'   inDir = example_honeybee_input(),
#'   f1_suffix = "_1.fastq.gz", f2_suffix = "_2.fastq.gz",
#'   outDir = paste0(example_honeybee_output(), ".test2"),
#'   singularity = "~/lab_shared/scripts/TAP/tap_latest.sif",
#'   reference = example_honeybee_reference()
#' )
#' config_write(cfg, "test_config.csv")
#' config_validate("test_config.csv")
#'
#' cfg_rename = config_create(
#'   inDir = example_honeybee_input(),
#'   f1_suffix = "_1.fastq.gz", f2_suffix = "_2.fastq.gz",
#'   outDir = paste0(example_honeybee_output(), ".rename2"),
#'   singularity = "~/lab_shared/scripts/TAP/tap_latest.sif",
#'   reference = example_honeybee_reference()
#' )
#' cfg_fq = cfg_rename$fastq_lines
#' colnames(cfg_fq) = c("file", "srr")
#' meta_df = example_honeybee_metadata()
#' cfg_fq = merge(meta_df, cfg_fq, by = "srr")
#' cfg_fq = cfg_fq[, c("file", "name")]
#' cfg_rename$fastq_lines = cfg_fq
#'
#' config_write(cfg_rename, "rename_config.csv")
#' config_validate("rename_config.csv")
config_create = function(
    inDir = NULL,
    f1_suffix = NULL,
    f2_suffix = NULL,
    outDir = NULL,
    reference = NULL,
    starIndex = NULL,
    suppaRef = NULL,
    gtf = NULL,
    fasta = NULL,
    rDNA_starIndex = NULL,
    PE = NULL,
    SE = NULL,
    noSub = NULL,
    pipeline = NULL,
    jobDir = NULL,
    scriptLocation = NULL,
    docker = NULL,
    singularity = NULL
){
  argg <- c(as.list(environment()))
  #add a couple defaults that make sense
  if(is.null(argg$inDir)){
    argg$inDir = getwd()
  }
  if(is.null(argg$f1_suffix)){
    argg$f1_suffix = "_R1_001.fastq.gz"
  }
  if(is.null(argg$f2_suffix)){
    argg$f2_suffix = "_R2_001.fastq.gz"
  }

  found_fq = dir(argg$inDir, pattern = paste0(argg$f1_suffix, "$"))
  if(length(found_fq) < 1){
    warning("No fastq files found, check inDir and f1_suffix.\ninDir: ", argg$inDir,"\nf1_suffix: ", argg$f1_suffix)
  }

  argg
  valid_args = character()
  for(nam in names(argg)){
    if(!is.null(argg[[nam]])){
      #these 2 are simple flags
      if(nam %in% c("PE", "SE", "noSub")){
        valid_args = c(valid_args, paste0("--", nam))
      }else{
        valid_args = c(valid_args, paste0("--", nam, " ", argg[[nam]]))
      }
    }
  }
  cfg_lines = paste("#CFG", valid_args)
  fq_names = sub(argg$f1_suffix, "", basename(found_fq))
  # fq_lines = paste(found_fq, fq_names, sep = ",")
  fq_lines = data.frame(file = found_fq, name = fq_names)
  message(paste(cfg_lines, collapse = "\n"))
  list(config_lines = cfg_lines, fastq_lines = fq_lines)
}

#' config_write
#'
#' Use with output from [config_create]. Modify fastq names as necessary.
#'
#' @param config_list list returned from [config_create]
#' @param config_file file to write configuration to, a specialized .csv format.
#' @param overwrite If TRUE, existing config file will be overwritted. Default is FALSE.
#'
#' @return path to configuration file.
#' @export
#'
#' @rdname config-write
#'
config_write = function(config_list, config_file, overwrite = FALSE){
  stopifnot(setequal(names(config_list), c("config_lines", "fastq_lines")))
  if(!overwrite){
    stopifnot(!file.exists(config_file))
  }

  line1 = paste("# written by TAPhelpR version", packageDescription("TAPhelpR")$Version, "on", date())
  fq_header = "# fastq file, sample name"
  writeLines(c(line1, config_list$config_lines, fq_header), config_file)
  write.table(config_list$fastq_lines, config_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
  config_file
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
      # arg_match = grepl(v_key, head_args)
      arg_match = v_key == head_args
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
  # if primary reference location is not specified then all other references must be
  if(!"reference" %in% names(cfg_params)){
    if(!all(c("starIndex", "suppaRef", 'gtf', "fasta") %in% names(cfg_params))){
      stop("When reference is not specified then starIndex, suppaRef, gtf, and fasta must be specified.")
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
#' @rdname config-write
config_validate = function(config_file, extra_args = NULL, work_dir = NULL){
  old_wd = getwd()
  if(!is.null(work_dir)){
    setwd(work_dir)
  }
  #tryCatch so that work_dir gets reset regardless of error
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
