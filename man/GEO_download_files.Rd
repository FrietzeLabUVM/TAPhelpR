% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_scrape_geo.R
\name{GEO_download_files}
\alias{GEO_download_files}
\title{GEO_download_files}
\usage{
GEO_download_files(
  srr_tofetch,
  fastq_prefixes,
  out_dir = getwd(),
  docker = NULL,
  singularity = NULL,
  bash_or_sbatch = "sbatch",
  return_commands = FALSE,
  ncores = 1
)
}
\arguments{
\item{srr_tofetch}{SRR accessions. SRR#####.}

\item{fastq_prefixes}{Prefixes for final fastq files. "_R1_001.fastq.gz" or "_R2_001.fastq.gz" will be appended as appropriate.}

\item{out_dir}{Location where fastq files will be downloaded to. Will be created if necessary.}

\item{docker}{Optional name of docker image, REPOSITORY:TAG. docker will be used if provided.}

\item{singularity}{Optional path to singularity .sif file. singularity will be used if provided.}

\item{bash_or_sbatch}{How the script will be called. "bash" or "sbatch" are the most likely values. Additional flags can be passed to sbatch or similar here. For example, if using ncores > 1 you need to tell SLURM about this; i.e. bash_or_sbatch = "sbatch -c 4"}

\item{return_commands}{If TRUE, commands are not executed via \link{system} but returned as a character vector.}
}
\value{
Invisibly returns the output directory where fastqs have been downloaded.
}
\description{
Uses \link{GEO_get_file_info} to get the required information for \code{srr_tofetch} and \code{fastq_prefixes}.
}
\examples{
# you need the GSE accession
gse_df = GEO_get_file_info("GSE152028")
# you typically need to do some cleanup from the title to generate usable file prefixes.
gse_df = dplyr::mutate(
  gse_df,
  name = sub(" \\\\+ ", "&", title)
)
gse_df = dplyr::mutate(
  gse_df,
  name = sub(" - ", " ", name)
)
gse_df = dplyr::mutate(
  gse_df,
  name = sub("input_DNA", "input", name)
)
gse_df$name = gsub(" +", "_", gse_df$name)
gse_df$name = paste0(gse_df$name, ".", gse_df$srr)
GEO_download_files(gse_df$srr, fastq_prefixes = gse_df$name, singularity = "tap_latest.sif", bash_or_sbatch = "bash")
}
