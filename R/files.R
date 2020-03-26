##' Access extdata files whether they're in inst or not
##'
##' To build a package files should live in inst/extdata,
##' whereas to refer to them in the package, they should have the
##' inst in their path. devtools::system.file deals with this.
##' This function is a helper for the devtools-modified system.file but
##' with the added benefit that you can put whatever you want into "extdata" and
##' whatever you ever use will also get linked to from "inst/extdata" so that
##' extra files in extdata that aren't loaded by the package can easily be excluded.
##' 
##' @title Intelligent use of extdata files
##' @param ... The path elements, without any need for the 'extdata' part
##' @return The file-path, and a side-effect of creating a link in inst/extdata if necessary
##' @author Gavin Kelly
extdata <- function(...) {
  fname <- system.file("extdata", ..., package=getPackageName())
  if (!grepl(paste0(normalizePath(file.path("inst", "extdata", ...)), "$"), fname)) {
    file.symlink(fname, sub("/extdata/", "/inst/extdata/", fname))
  }
  fname
}

##' Cache data/ recover cached data
##'
##' Wrap this around an assignment.  If the variable is already in cache, recover it, otherwise
##' recalculate it and put it in cache.  Non-interactive sessions always recalculate from fresh.
##' The cache is stored in the 'data' directory, so as to provide ready-made data objects in the
##' eventual library.  They're stored with an underscore at the beginning of the file extension
##' as devtools::load_all doesn't lazy-load data, which can slow things down.
##' @title Intelligent library-targetted caching
##' @param expr The assignment to be cached/recovered
##' @return The assignment is made in the correct environment
##' @author Gavin Kelly
dcache <- function (expr) 
{
    pexpr <- parse(text = deparse(substitute(expr)))
    pexpr <- as.list(pexpr[[1]])
    name <- as.character(pexpr[[2]])
    RHS <- pexpr[[3]]
    cachefile <- file.path("data", paste0(name, "._RData"))
    if (file.exists(cachefile) & interactive()) {
      message("Loading from cache")
      load(cachefile)
      assign(name, get(name), envir = parent.frame())
    }
    else {
        assign(name, eval(RHS, envir = parent.frame()), envir = parent.frame())
        save(list = name, file = cachefile, envir = parent.frame())
    }
    invisible(name)
}

