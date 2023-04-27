#' make_make_contrasts
#'
#' @param config_df Data.frame containing metadata with grouping variables for
#'   differential analysis.
#' @param groups_var Defines column with grouping assignments within which
#'   comparisons will be made.
#' @param vary_var Defines column that identifies treatment variable between
#'   which comparisons will be made.
#'
#' @return A list with 2 named elements. A "contrast" matrix produced by limma::makeContrasts and the "design" matrix that assigns items in config to groups.
#' @export
#'
#' @examples
#' cfg_file = system.file(package = "TAPhelpR", "extdata/test_diff_config_cnr.csv", mustWork = TRUE)
#' cfg_df = read.table(cfg_file, sep = ",", header = TRUE)
#' cfg_df = subset(cfg_df, mark == "H3K27ac")
#' cfg_df
#' make_make_contrasts(cfg_df, groups_var = "background", vary_var = "treatment")
make_make_contrasts = function(config_df, groups_var = "cell", vary_var = "treatment"){
  stopifnot(is.data.frame(config_df))
  config_df = as.data.frame(config_df)
  vars = c(groups_var, vary_var)

  group = factor(apply(config_df[, vars], 1, function(x)paste(x, collapse = "_")))
  design = model.matrix(~0+group)

  cont = character()
  cont_names = character()
  vary_values = unique(config_df[[vary_var]])
  group_values = unique(config_df[[groups_var]])

  num_groups = length(unique(config_df[[groups_var]]))
  num_vary = length(vary_values)
  for(v_i in seq(1, num_vary - 1)){
    for(v_j in seq(v_i + 1, num_vary)){
      if(num_groups == 1){
        cont_names = c(cont_names, paste(sep = "_", vary_values[v_j], "vs", vary_values[v_i], "in", group_values[1]))
        cont = c(cont, paste0("group", group_values[1], "_", vary_values[v_j], "-", "group", group_values[1], "_", vary_values[v_i]))
      }else if(num_groups == 2){

        cont_names = c(cont_names, paste(sep = "_", vary_values[v_j], "vs", vary_values[v_i], "in", "all"))
        cont = c(cont, paste0("(group", group_values[2], "_", vary_values[v_j], "+group", group_values[1], "_", vary_values[v_j], ")/2",
                              "-",
                              "(group", group_values[2], "_", vary_values[v_i], "+group", group_values[1], "_", vary_values[v_i], ")/2")
        )

        cont_names = c(cont_names, paste(sep = "_", vary_values[v_j], "vs", vary_values[v_i], "in", group_values[2]))
        cont = c(cont, paste0("group", group_values[2], "_", vary_values[v_j], "-", "group", group_values[2], "_", vary_values[v_i]))

        cont_names = c(cont_names, paste(sep = "_", vary_values[v_j], "vs", vary_values[v_i], "in", group_values[1]))
        cont = c(cont, paste0("group", group_values[1], "_", vary_values[v_j], "-", "group", group_values[1], "_", vary_values[v_i]))
      }
    }
  }


  cbind(cont_names, cont)
  names(cont) = cont_names
  cbind(colnames(design))

  cont
  colnames(design)

  all_contrasts = limma::makeContrasts(contrasts = cont, levels = colnames(design))
  colnames(all_contrasts) = names(cont)
  list(contrasts = all_contrasts, design = design)
}
