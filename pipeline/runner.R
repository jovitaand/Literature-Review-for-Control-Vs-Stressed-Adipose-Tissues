# Runs notebook cells (cells/cell_NNN.R) one after another in ONE R session, with logging, timing, memory and checkpoints.
# Usage: Rscript runner.R RUN_DIR START STOP_AFTER CHECKPOINT_MODE CELL_TIMEOUT_SECONDS
# All runner variables start with ".pl_" and are never saved in checkpoints, so restoring a checkpoint cannot overwrite them.
.pl_args <- commandArgs(TRUE); .pl_rd <- .pl_args[1]; .pl_start <- as.integer(.pl_args[2]); .pl_stop <- as.integer(.pl_args[3]); .pl_ckpt <- .pl_args[4]; .pl_tmo <- as.numeric(.pl_args[5])
.pl_cells <- sort(list.files(file.path(.pl_rd, "cells"), "^cell_[0-9]+\\.R$", full.names = TRUE)); .pl_n <- length(.pl_cells)
dir.create(file.path(.pl_rd, "checkpoints"), showWarnings = FALSE); dir.create(file.path(.pl_rd, "figures"), showWarnings = FALSE)
.pl_peak <- function() { x <- tryCatch(readLines("/proc/self/status"), error = function(e) character()); v <- grep("^VmHWM", x, value = TRUE); if (length(v)) round(as.numeric(gsub("[^0-9]", "", v)) / 1024) else NA }
.pl_log <- list(); .pl_state <- "running"; .pl_failed <- "null"; .pl_last <- .pl_start - 1
.pl_write <- function() {
  cs <- vapply(.pl_log, function(c) sprintf('{"cell": %d, "seconds": %.1f, "peak_mem_mb": %s, "ok": %s}', c$cell, c$seconds, ifelse(is.na(c$mem), "null", c$mem), tolower(c$ok)), "")
  f <- file.path(.pl_rd, "status.json"); old <- character()      # keep the cell history of earlier (resumed) segments
  if (file.exists(f) && .pl_start > 1) { t <- paste(readLines(f), collapse = ""); old <- regmatches(t, gregexpr('\\{"cell": [0-9]+[^}]*\\}', t))[[1]]; old <- old[as.integer(sub('\\{"cell": ([0-9]+).*', "\\1", old)) < .pl_start] }
  writeLines(sprintf('{"total_cells": %d, "last_completed": %d, "state": "%s", "failed_cell": %s, "cells": [%s]}', .pl_n, .pl_last, .pl_state, .pl_failed, paste(c(old, cs), collapse = ", ")), f)
}
OUT_DIR <- file.path(normalizePath(.pl_rd), "out"); RUN_DIR <- normalizePath(.pl_rd)   # defaults; a notebook can reassign them
pdf(file.path(.pl_rd, "figures", "default_device.pdf"))                 # plots printed without an explicit png() end up here
if (.pl_start > 1) {                                                    # resume: restore the workspace and the attached packages
  .pl_ck <- file.path(.pl_rd, "checkpoints", sprintf("after_%03d.RData", .pl_start - 1)); if (!file.exists(.pl_ck)) stop("No checkpoint ", .pl_ck, "; use --from-cell 1 or checkpoint mode 'each'")
  load(.pl_ck, envir = globalenv()); .pl_pk <- sub("RData$", "pkgs", .pl_ck)
  if (file.exists(.pl_pk)) for (.pl_p in readLines(.pl_pk)) suppressPackageStartupMessages(library(.pl_p, character.only = TRUE))
  cat(sprintf("[runner] resumed from checkpoint after cell %d\n", .pl_start - 1))
}
for (.pl_i in .pl_start:min(.pl_n, .pl_stop)) {
  cat(sprintf("\n[runner] ===== cell %d of %d : %s =====\n", .pl_i, .pl_n, format(Sys.time(), "%H:%M:%S"))); .pl_t0 <- Sys.time()
  .pl_ok <- tryCatch({ setTimeLimit(elapsed = .pl_tmo, transient = TRUE)
    source(.pl_cells[.pl_i], local = globalenv(), echo = TRUE, max.deparse.length = Inf, keep.source = TRUE, print.eval = TRUE); TRUE },
    error = function(e) { cat("[runner] ERROR in cell", .pl_i, ":", conditionMessage(e), "\n"); FALSE })
  setTimeLimit(elapsed = Inf, transient = TRUE)
  .pl_sec <- as.numeric(difftime(Sys.time(), .pl_t0, units = "secs")); .pl_mem <- .pl_peak()
  .pl_log[[length(.pl_log) + 1]] <- list(cell = .pl_i, seconds = .pl_sec, mem = .pl_mem, ok = .pl_ok)
  cat(sprintf("[runner] cell %d %s in %.1f s (peak memory %s MB)\n", .pl_i, ifelse(.pl_ok, "done", "FAILED"), .pl_sec, .pl_mem))
  if (!.pl_ok) { .pl_state <- "failed"; .pl_failed <- as.character(.pl_i); .pl_write(); quit(status = 1) }
  .pl_last <- .pl_i
  if (.pl_ckpt == "each" || (grepl("^[0-9]+$", .pl_ckpt) && .pl_i %% as.integer(.pl_ckpt) == 0)) {
    .pl_keep <- setdiff(ls(globalenv(), all.names = TRUE), ls(globalenv(), all.names = TRUE, pattern = "^\\.pl_"))
    save(list = .pl_keep, file = file.path(.pl_rd, "checkpoints", sprintf("after_%03d.RData", .pl_i)))
    writeLines(setdiff(.packages(), c("stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")), file.path(.pl_rd, "checkpoints", sprintf("after_%03d.pkgs", .pl_i)))
    .pl_old <- list.files(file.path(.pl_rd, "checkpoints"), "^after_[0-9]+\\.(RData|pkgs)$", full.names = TRUE); unlink(.pl_old[!grepl(sprintf("after_%03d", .pl_i), .pl_old)])   # keep only the newest checkpoint
  }
  .pl_write()
}
.pl_state <- if (.pl_last >= .pl_n) "finished" else "paused"; .pl_write(); cat("[runner] state:", .pl_state, "\n")
