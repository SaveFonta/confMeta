# get the old function
get_old_FUN <- function(path, fun_name) {
  # Source the utils file as well
  source(
    "https://raw.githubusercontent.com/felix-hof/hMean/main/R/utils.R",
    local = TRUE
  )
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
