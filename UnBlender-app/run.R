# === Parse Arguments ===
args <- commandArgs(TRUE)
wd <- args[1]  # Shiny app directory
debug <- tolower(args[2]) == "true"  # Debug flag

rd <- dirname(wd)  # Root directory

# === Logging Function ===
log_msg <- function(...) {
  if (debug) message(...)
}

# === Abort if no working directory ===
if (is.na(wd)) {
  stop("Please use run.bat to launch the application.")
}

# === Set CRAN Mirror ===
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# === Determine Library Path ===
lib_path <- if (
  file.exists(file.path(rd, "R", "bin", "R.exe")) &&
    .libPaths()[1] != file.path(rd, "R", "library")
) {
  file.path(rd, "R", "library")
} else {
  .libPaths()[1]
}
log_msg("Using library path: ", lib_path)

# === Install and Load CRAN Packages ===
install_and_load_cran_packages <- function(req_file, lib_path) {
  if (!file.exists(req_file)) {
    log_msg("No req.txt found.")
    return()
  }

  req <- readLines(req_file)
  log_msg("Required packages: ", paste(req, collapse = ", "))

  # Check which packages are already installed
  installed <- installed.packages(lib.loc = lib_path)[, "Package"]
  missing <- req[!(req %in% installed)]
  
  # Install what is missing
  if (length(missing) > 0) {
    log_msg("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, lib = lib_path, clean = TRUE, type = "source")
  }

  suppressPackageStartupMessages(
    invisible(lapply(req, function(pkg) {
      log_msg("Loading package: ", pkg)
      library(pkg, character.only = TRUE)
    }))
  )
}

# === Install Bioconductor Packages ===
install_bioc_packages <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    # log_msg("Installing BiocManager...")
    install.packages("BiocManager")
  }

  missing_bioc <- pkgs[!(pkgs %in% installed.packages(lib.loc = lib_path)[, "Package"])]
  if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
  }
  
}

install_github_package <- function(repo, commit, pkg_name) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  installed <- tolower(installed.packages(lib.loc = lib_path)[, "Package"])
  if (tolower(pkg_name) %in% installed) {
    log_msg("GitHub package already installed: ", pkg_name)
    return()
  }

  url <- paste0(repo, "@", commit)
  log_msg("Installing GitHub package: ", url)
  remotes::install_github(url, lib = lib_path, upgrade = "never")
}

# === Launch Shiny App ===
launch_app <- function(app_dir, browser_path = NULL) {
  if (!is.null(browser_path) && file.exists(browser_path)) {
    log_msg("Launching app in embedded browser...")
    shiny::runApp(
      appDir = app_dir,
      launch.browser = function(url) {
        system(paste0("\"", browser_path, "\" --app=", url, " -incognito --no-sandbox"), wait = FALSE)
      }
    )
  } else {
    log_msg("Launching app in default browser...")
    shiny::runApp(appDir = app_dir, launch.browser = TRUE)
  }
}

# === Main Execution ===
install_and_load_cran_packages(file.path(rd, "req.txt"), lib_path)

install_bioc_packages(c("SingleCellExperiment", "TOAST"))

install_github_package("xuranw/MuSiC", "f21fe67f5670d5e9fca0ad7550abaae3423eb59c", "MuSiC")


launch_app(
  app_dir = wd,
  browser_path = file.path(rd, "chrome", "chrome.exe")
)
