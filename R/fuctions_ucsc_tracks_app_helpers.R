grab_files = function(root_dir, pattern = "\\.bw"){
  dir(root_dir, recursive = T, pattern = pattern, full.names = T)
}

grab_files.test = function(test_file){
  read.table(test_file, stringsAsFactors = FALSE)[,1]
}

clean_files_names = function(file_paths){
  # keep = !duplicated(basename(file_paths))
  # file_paths = file_paths[keep]
  file_df = data.frame(files = file_paths, directory = basename(dirname(file_paths)), raw_names = basename(file_paths))
  file_df$clean_names = file_df$raw_names
  #do a bunch of stuff to cleanup raw names
  file_df$clean_names = sub(".bw$", "", file_df$clean_names)
  file_df$clean_names = sub("_bigwig$", "", file_df$clean_names)
  file_df$clean_names = sub("].bigwig$", "", file_df$clean_names)

  file_df$clean_names = sub(".bb$", "", file_df$clean_names)
  file_df$clean_names = sub(".bigbed$", "", file_df$clean_names)
  file_df$clean_names = sub(".bigBed$", "", file_df$clean_names)

  file_df$clean_names = sub("Galaxy[0-9]+-\\[", "", file_df$clean_names)

  file_df$clean_names = sub("TM_", "", file_df$clean_names)
  file_df$clean_names = sub("_Control_", "_ctrl_", file_df$clean_names)
  file_df$clean_names = sub("MDA-MB-231", "MDA231", file_df$clean_names)
  #cells
  file_df$clean_names = sub("MCF10CA1", "MCF10A-CA1", file_df$clean_names)
  file_df$clean_names = sub("MCF10AT1", "MCF10A-AT1", file_df$clean_names)

  #drugs
  file_df$clean_names = sub("_BZA_E2_", "_e2bza_", file_df$clean_names)
  file_df$clean_names = sub("_BZAE2_", "_e2bza_", file_df$clean_names)
  file_df$clean_names = sub("_BZA_GC10_", "_gc10bza_", file_df$clean_names)
  file_df$clean_names = sub("_BZA_", "_bza_", file_df$clean_names)
  file_df$clean_names = sub("_E2_", "_e2_", file_df$clean_names)
  file_df$clean_names = sub("_GC10_", "_gc10_", file_df$clean_names)

  file_df$clean_names = sub("_rep", "_R", file_df$clean_names)

  file_df$clean_names = gsub("\\.", "_", file_df$clean_names)
  file_df
}

split_file_names = function(file_df){
  raw_split = strsplit(as.character(file_df$clean_names), "_")
  MAX = max(sapply(raw_split, length))
  mat_split = t(sapply(raw_split, function(x){
    return(c(x, rep("", MAX - length(x))))
  }))
  mat_split = as.data.frame(mat_split)

  bw_final_mat = cbind(file_df, mat_split)
  bw_final_mat$clean_names = NULL
  bw_final_mat$raw_names = NULL
  bw_final_mat = bw_final_mat[, c(seq_len(ncol(bw_final_mat))[-1], 1)]
  bw_final_mat
}

makeTrack = function(base_host_path, base_url, file, type, name, description, size, view_min, view_max,
                     color, visibility, autoscale, alwaysZero, showZero, showYline,
                     valYline, smoothingWindow, windowFun){
  url = sub(base_host_path, base_url, file)
  new_track = TEMPLATE
  new_track = sub("TYPE", type, new_track)
  new_track = sub("URL", url, new_track)
  new_track = sub("NAME", name, new_track)
  new_track = sub("DESCRIPTION", description, new_track)
  new_track = sub("SIZE", size, new_track)
  new_track = sub("MIN", view_min, new_track)
  new_track = sub("MAX", view_max, new_track)
  new_track = sub("COLOR", color, new_track)
  new_track = sub("VISIBILITY", UCSC_VIS[visibility], new_track)
  new_track = sub("AUTOSCALE", ifelse(autoscale, "on", "off"), new_track)
  new_track = sub("ALWAYSZERO", ifelse(alwaysZero, "on", "off"), new_track)
  new_track = sub("UCSC_SHOWYLINE", ifelse(showYline, "on", "off"), new_track)
  new_track = sub("UCSC_VALYLINE", valYline, new_track)
  new_track = sub("UCSC_SHOWZERO", ifelse(showZero, "on", "off"), new_track)
  new_track = sub("UCSC_SMOOTHING", ifelse(smoothingWindow > 1, smoothingWindow, "off"), new_track)
  new_track = sub("UCSC_WINDOWFUN", windowFun, new_track)
  return(new_track)
}

safeBrew = seqsetvis::safeBrew

hex2ucsc_rgb = function(colorChoices){
  apply(col2rgb(colorChoices),2, function(x)paste(x, collapse = ","))
}

gen_color_picker_ui = function(grps, brew_name, is_free_color){
  possibleColors = unique(safeBrew(50, brew_name))
  colorChoices = safeBrew(length(grps), brew_name)
  #special case for gradients where all color choices aren't in possible colors
  if(!all(colorChoices %in% possibleColors)){
    clen = length(possibleColors)
    glen = length(grps) + 1
    gi = round((1:glen-1) * ((clen-1) / (glen - 1)) + 1)
    colorChoices = possibleColors[gi]
  }

  names(colorChoices) = grps
  inputs = character(length(grps))

  for(i in seq_along(grps)){
    g = as.character(grps[i])
    g = gsub("/", "_", g)

    if(is_free_color){
      ginput = colourpicker::colourInput(inputId = paste0("color_", g),
                                         label = g,
                                         value = colorChoices[i])
    }else{
      ginput = colourpicker::colourInput(inputId = paste0("color_", g),
                                         label = g,
                                         value = colorChoices[i],
                                         palette = "limited",
                                         allowedCols = possibleColors)
    }

    inputs[i] = as.character(ginput)
  }
  HTML(paste(inputs, collapse = "</br>"))
}

refresh_bw = function(root_dir){
  split_file_names(clean_files_names(grab_files(root_dir = root_dir, pattern = "\\.bw|\\.bigwig|\\.bigWig")))
}
refresh_bb = function(root_dir){
  split_file_names(clean_files_names(grab_files(root_dir = root_dir, pattern = "\\.bb|\\.bigbed|\\.bigBed")))
}
