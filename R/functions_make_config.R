fq_files = dir("zhang_IK_SE/fastqs/", pattern = "fastq.gz", full.names = TRUE)
fq_files
cfg_lines = c(
  "#CFG -ref ~/lab_shared/indexes/MM10",
  "#CFG -i zhang_IK_SE/fastqs -o zhang_IK_SE/output"
)
fq_files
df = data.frame(fq = fq_files)
library(tidyverse)
df = df %>%
  mutate(name = basename(fq)) %>%
  separate(name, sep = "[_\\.]", into = c("org", "cell", "mark", "treatment"), extra = "drop")
df$rep = "rep1"
df = df %>%
  mutate(fq_name = basename(fq)) %>%
  mutate(rep_name = paste(org, cell, mark, treatment, rep, sep = "_")) %>%
  mutate(pool_name = paste(org, cell, mark, treatment, "pooled", sep = "_")) %>%
  mutate(input_name = paste(org, cell, "input", treatment, "pooled", sep = "_"))
df = df %>% select(fq_name, rep_name, pool_name, input_name)

stopifnot(nrow(df %>% filter(pool_name == input_name)) > 0)

body_lines = df %>% mutate(line = paste(fq_name, rep_name, pool_name, input_name, sep = ","), .keep = "none")
body_lines = body_lines[[1]]

cfg_f = "config_zhang_IK_SE.csv"
writeLines(c(cfg_lines, body_lines), cfg_f)
