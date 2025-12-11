### Miscellaneous utility functions for the pipeline ###


# Extracts the data for a given trial from a dataframe

get_trial <- function(dat, p_num, s_num, b_num, t_num) {
  subset(dat, id == p_num & session == s_num & block == b_num & trial == t_num)
}


# Prints a progress message to the console when running scripts with 'source()'

progress_msg <- function(msg, header = FALSE) {
  if (header) {
    cat("\n=== ", msg, " ===\n\n")
  } else {
    cat(" * ", msg, "...\n")
  }
}


# Creates a new folder, removing it first if it already exists

mkdir_replace <- function(outdir) {
  if (outdir != ".") {
    if (dir.exists(outdir)) {
      unlink(outdir, recursive = TRUE)
    }
    dir.create(outdir, recursive = TRUE)
  }
}
