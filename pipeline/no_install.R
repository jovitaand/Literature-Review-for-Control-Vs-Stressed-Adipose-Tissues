# Loaded by R before the notebook runs when CRAN is not reachable: turns install.packages() into a clear message instead of a hang.
local({ r <- getOption("repos"); options(repos = c(CRAN = "http://127.0.0.1:1")); options(timeout = 5) })
utils::assignInNamespace("install.packages", function(pkgs, ...) { message("[pipeline] install.packages(", paste(pkgs, collapse = ", "), ") skipped: package source is offline. Run pipeline/setup_env.sh"); invisible(NULL) }, "utils")
