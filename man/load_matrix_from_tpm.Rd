% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_salmon_tpm.R
\name{load_matrix_from_tpm}
\alias{load_matrix_from_tpm}
\title{load_matrix_from_tpm}
\usage{
load_matrix_from_tpm(files)
}
\arguments{
\item{files}{paths to tpm.txt files generated by TAP in the *salmon_quant directories.}
}
\value{
matrix of TPM values generated by salmon.
}
\description{
load_matrix_from_tpm
}
\examples{
tap_out = "~/R_workspace.combined/TAPhelpR.data/honeybee_TAP_output"
tpm_files = setup_tpm_files(tap_out, variable_map = c("day", "sex", "rep"))
mat = load_matrix_from_tpm(tpm_files$file)
head(mat)
}
