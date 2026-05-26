install_github_package <- function(repo, commit, pkg_name) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  url <- paste0(repo, "@", commit)
  remotes::install_github(url, upgrade = "never")
}

install_github_package("xuranw/MuSiC", "f21fe67f5670d5e9fca0ad7550abaae3423eb59c", "MuSiC")
