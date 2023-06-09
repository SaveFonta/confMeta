# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(hMean)

get_old_FUN <- function(path, fun_name) {
  new_name <- paste0(fun_name, "_old")
  def <- readLines(path)
  comments <- grepl("^\\s*#", def)
  def <- def[!comments]
  def <- paste0(def, collapse = "\n")
  regex <- paste0(fun_name, "\\s*<-\\s*function")
  replace <- paste0(new_name, " <- function")
  def <- sub(regex, replace, def)
  eval(parse(text = def))
  get(new_name)
}

test_check("hMean")