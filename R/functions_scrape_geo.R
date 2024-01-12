
get_right = function(key, txt){
  split = sapply(strsplit(txt, key), function(x)paste(x[-1], collapse = key))
  return(split)
}
get_left = function(key, txt){
  split = sapply(strsplit(txt, key), function(x)x[1])
  return(split)
}
get_mid = function(key_left, key_right, txt){
  right_split = get_right(key_left, txt)
  mid_split = get_left(key_right, right_split)
  return(mid_split)
}

get_attrib = function(gsm_lines, key, w = NULL, debug = FALSE){
  k = which(grepl(key, gsm_lines))
  attrib_line = gsm_lines[k+1]
  split = strsplit(attrib_line, split = "[\\\"><]")[[1]]
  if(debug) print(split)
  if(is.null(w)){
    split[which(split == "/td") - 1]
  }else if(is.character(w)){
    split[which(grepl(w, split))[1]]
  }else{
    split[w]
  }
}

scrape_gse = function(
    gsms
){
  base_url = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
  gse = (parallel::mclapply(gsms, function(g){
    url = paste0(base_url, g)
    con = curl::curl(url = url)
    gsm_lines = readLines(con)
    close(con)
    gse_line = gsm_lines[grepl("GSE", gsm_lines)]
    gse = regmatches(gse_line, regexpr("GSE[0-9]+", gse_line))
    gse
  }))
  names(gse) = gsms
  gse
}

readLines.recursive = function(url, wait = .5){
  success = FALSE
  tryCatch({
    con = curl::curl(url = url)
    lines = readLines(con)
    close(con)
    success = TRUE
  }, error = function(e){
    message(e)

  })
  if(!success){
    lines = readLines.recursive(url, wait)
  }
  lines
}

#' GEO_get_file_info
#'
#'
#' Characteristics
#'   cell type: stromal adherent pre-B cell progenitors
#'   genotype: IKDN CD2 Cre
#'   passages: 3 to 8
#'   antibody: Abnova/H00001879-M01
#'
#' @return A data.table containing information on .fastq files available via SRA associated with GEO GSE series accession. Should used to call [GEO_download_files].
#' @export
#' @import curl
#' @examples
#' my_gse = "GSE86897"
#' gse_dt = GEO_get_file_info(my_gse)
#' gse_dt
#' stopifnot(!any(duplicated(gse_dt$title)))
#' dump_script = system.file("extdata/fasterq_dump_wrapper.sh", package = "TAPhelpR", mustWork = TRUE)
#' srr_tofetch =
#' library(data.table)
#'
#' gse_dt[, title := sub("Badh", "_Bad", title)]
#' gse_dt[, title := sub("InputControl_", "input_", title)]
#' gse_dt[, title := sub("pe_", "pre_", title)]
#' gse_dt[, title := sub("pre_", "pre", title)]
#' gse_dt[, title := sub("Pre", "pre", title)]
#' gse_dt[, title := sub("pre_", "pre", title)]
#' gse_dt[, title := sub("pre", "_pre", title)]
#' gse_dt[, title := gsub("_+", "_", title)]
#' gse_dt[, title := sub("H3K27Ac", "H3K27ac", title)]
#' gse_dt[, c("mark", "treatment", "cell", "method") := tstrsplit(title, "_", keep = 1:4)]
#'
#' table(gse_dt$mark)
#' table(gse_dt$treatment)
#' table(gse_dt$cell)
#' gse_dt$cell = "preBad"
#' gse_dt[, treatment := sub("IkDN", "IKDN", treatment)]
#'
#' table(gse_dt$mark)
#' table(gse_dt$treatment)
#' table(gse_dt$cell)
#'
#' gse_dt[, prefix := paste("mouse", cell, mark, treatment, sep = "_")]
#'
#' fastq_prefixes = gse_dt$prefix
#' out_dir = "zhang_IK_SE/fastqs"
#' fetch_cmds = submit_sra_fetch(gse_dt$srr, gse_dt$prefix, out_dir)#, docker = "jrboyd/tap")
#' writeLines(
#' sapply(fetch_cmds, system),
#' "submit_fetch_zhange.sh")
GEO_get_file_info = function(GSE_id,
                             gsm_override = NULL,
                             gsm_todo = NULL,
                             gsm_skip = NULL,
                             do_srr = TRUE,
                             debug = FALSE,
                             chara = c("cancer type", "cell type", "genotype", "passages", "antibody"),
                             key_str = c("Title", rep("Characteristics", length(chara))),
                             key_idx = append(list(5), chara)){
  stopifnot(length(key_str) == length(key_idx))



  if(!is.null(gsm_override)){
    gsms = gsm_override
  }else{
    gse_url = paste0("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GSE_id)
    gse_lines = readLines.recursive(gse_url)

    keep = which(grepl("GSM", gse_lines))
    gsm_lines = gse_lines[keep]
    # gsms = get_mid(key_left = "geoaxema_recenter)\">", key_right = "</a></t", gsm_lines)
    if(!is.null(gsm_todo)){
      gsms = gsm_todo
    }else{
      gsms = sapply(gsm_lines, function(x){
        m = regexpr(pattern = "GSM[0-9]+", text = x)
        regmatches(x = x, m = m)
      })
    }
    gsms = setdiff(gsms, gsm_skip)
  }
  names(gsms) = NULL

  mat = matrix("", nrow = 0, ncol = 2+length(key_str) + sum(do_srr))
  # rownames(mat) = gsms
  base_url = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
  for(g in gsms){
    print(g)
    url = paste0(base_url, g)
    message(url)

    gsm_lines = readLines.recursive(url)

    if(debug) browser()

    vals = sapply(seq_along(key_str), function(ks){
      get_attrib(gsm_lines, key_str[ks], key_idx[[ks]], debug = debug)
    })
    if(do_srr){
      srx_line = gsm_lines[grepl("href", gsm_lines) & grepl("SRX", gsm_lines)][1]
      if(is.na(srx_line)){
        srr = "na"
      }else{
        m = regexpr(pattern = "http.+term=SRX[0-9]+", text = srx_line)
        srx_url = regmatches(x = srx_line, m = m)

        srx_lines = readLines.recursive(srx_url)
        m = gregexpr(pattern = "SRR[0-9]+", text = srx_lines)
        srr = unique(unlist(regmatches(x = srx_lines, m = m)))
      }
      newL = cbind(matrix(rep(c(GSE_id, g, vals), length(srr)), nrow = length(srr), byrow = TRUE), srr)
      print(newL)
      mat = rbind(mat, newL)
      #mat[g, ] = newL
    }else{
      mat = rbind(mat, vals)
      # mat[g, ] = vals
    }
    #   close(srx_con)
    #
    #   m = regexpr(pattern = "SRX[0-9].+", text = ftp_line)
    #   ftp = sub("\"", "", regmatches(x = ftp_line, m = m))
    #   ftp_con = curl::curl(url = ftp)
    #   readLines(ftp_con)
  }

  dt = data.table::as.data.table(mat)
  data.table::setnames(dt, paste0("V", seq(3+length(chara))), c("gse", "gsm", "title", unlist(chara)))
  for(ch in chara){
    dt[[ch]] = sub(".+: ", "", dt[[ch]])
  }
  return(dt)
}

#' GEO_download_files
#'
#' Uses [GEO_get_file_info] to get the required information for `srr_tofetch` and `fastq_prefixes`.
#'
#' @param srr_tofetch SRR accessions. SRR#####.
#' @param fastq_prefixes Prefixes for final fastq files. "_R1_001.fastq.gz" or "_R2_001.fastq.gz" will be appended as appropriate.
#' @param out_dir Location where fastq files will be downloaded to. Will be created if necessary.
#' @param docker Optional name of docker image, REPOSITORY:TAG. docker will be used if provided.
#' @param singularity Optional path to singularity .sif file. singularity will be used if provided.
#' @param bash_or_sbatch How the script will be called. "bash" or "sbatch" are the most likely values. Additional flags can be passed to sbatch or similar here. For example, if using ncores > 1 you need to tell SLURM about this; i.e. bash_or_sbatch = "sbatch -c 4"
#' @param return_commands If TRUE, commands are not executed via [system] but returned as a character vector.
#'
#' @return Invisibly returns the output directory where fastqs have been downloaded.
#' @export
#'
#' @examples
#' # you need the GSE accession
#' gse_df = GEO_get_file_info("GSE152028")
#' # you typically need to do some cleanup from the title to generate usable file prefixes.
#' gse_df = dplyr::mutate(
#'   gse_df,
#'   name = sub(" \\+ ", "&", title)
#' )
#' gse_df = dplyr::mutate(
#'   gse_df,
#'   name = sub(" - ", " ", name)
#' )
#' gse_df = dplyr::mutate(
#'   gse_df,
#'   name = sub("input_DNA", "input", name)
#' )
#' gse_df$name = gsub(" +", "_", gse_df$name)
#' gse_df$name = paste0(gse_df$name, ".", gse_df$srr)
#' GEO_download_files(gse_df$srr, fastq_prefixes = gse_df$name, singularity = "tap_latest.sif", bash_or_sbatch = "bash")
GEO_download_files = function(srr_tofetch,
                              fastq_prefixes,
                              out_dir = getwd(),
                              docker = NULL,
                              singularity = NULL,
                              bash_or_sbatch = "sbatch",
                              return_commands = FALSE,
                              ncores = 1){
  SRA_download_files(srr_tofetch = srr_tofetch,
                     fastq_prefixes = fastq_prefixes,
                     out_dir = out_dir,
                     docker = docker,
                     singularity = singularity,
                     bash_or_sbatch = bash_or_sbatch,
                     return_commands = return_commands,
                     ncores = ncores)
}

.create_sra_fetch_cmds = function(srr_tofetch, fastq_prefixes, out_dir = getwd(), docker = NULL, singularity = NULL, ncores = 1){
  dump_script = system.file("extdata/fasterq_dump_wrapper.sh", package = "TAPhelpR", mustWork = TRUE)
  stopifnot(length(srr_tofetch) == length(fastq_prefixes))
  stopifnot(!any(duplicated(srr_tofetch)))
  stopifnot(!any(duplicated(fastq_prefixes)))
  if(!is.null(docker) && !is.null(singularity)){
    stop("only one of docker or singularity is allowed")
  }
  dir.create(out_dir, showWarnings = FALSE)
  i = 1
  all_cmds = character()
  for(i in seq_along(srr_tofetch)){
    srr = srr_tofetch[i]
    pre = fastq_prefixes[i]
    cmd_args = paste("-s", srr, "-p", pre, "-o", out_dir, "-n", ncores)
    if(!is.null(docker)){
      cmd_args = paste(cmd_args, "--docker", docker)

    }
    if(!is.null(singularity)){
      cmd_args = paste(cmd_args, "--singularity", singularity)
    }
    cmd = paste(dump_script, cmd_args)
    all_cmds = c(all_cmds, cmd)
  }
  all_cmds
}

#' SRA_download_files
#'
#' Uses [GEO_get_file_info] to get the required information for `srr_tofetch` and `fastq_prefixes`.
#'
#' @param srr_tofetch SRR accessions. SRR#####.
#' @param fastq_prefixes Prefixes for final fastq files. "_R1_001.fastq.gz" or "_R2_001.fastq.gz" will be appended as appropriate.
#' @param out_dir Location where fastq files will be downloaded to. Will be created if necessary.
#' @param docker Optional name of docker image, REPOSITORY:TAG. docker will be used if provided.
#' @param singularity Optional path to singularity .sif file. singularity will be used if provided.
#' @param bash_or_sbatch How the script will be called. "bash" or "sbatch" are the most likely values. Additional flags can be passed to sbatch or similar here. For example, if using ncores > 1 you need to tell SLURM about this; i.e. bash_or_sbatch = "sbatch -c 4"
#' @param return_commands If TRUE, commands are not executed via [system] but returned as a character vector.
#'
#' @return Invisibly returns the output directory where fastqs have been downloaded.
#' @export
#'
#' @examples
#' # you need the GSE accession
#' gse_df = GEO_get_file_info("GSE152028")
#' # you typically need to do some cleanup from the title to generate usable file prefixes.
#' gse_df = dplyr::mutate(
#'   gse_df,
#'   name = sub(" \\+ ", "&", title)
#' )
#' gse_df = dplyr::mutate(
#'   gse_df,
#'   name = sub(" - ", " ", name)
#' )
#' gse_df = dplyr::mutate(
#'   gse_df,
#'   name = sub("input_DNA", "input", name)
#' )
#' gse_df$name = gsub(" +", "_", gse_df$name)
#' gse_df$name = paste0(gse_df$name, ".", gse_df$srr)
#' SRA_download_files(gse_df$srr, fastq_prefixes = gse_df$name, singularity = "tap_latest.sif", bash_or_sbatch = "bash")
SRA_download_files = function(srr_tofetch,
                              fastq_prefixes,
                              out_dir = getwd(),
                              docker = NULL,
                              singularity = NULL,
                              bash_or_sbatch = "sbatch",
                              return_commands = FALSE,
                              ncores = 1){
  all_cmds = .create_sra_fetch_cmds(srr_tofetch, fastq_prefixes, out_dir, docker, singularity, ncores)
  if(return_commands){
    sapply(all_cmds, function(cmd){
      paste(bash_or_sbatch, cmd)
    })
  }else{
    sapply(all_cmds, function(cmd){
      system(paste(bash_or_sbatch, cmd))
    })
  }
  invisible(out_dir)
}
