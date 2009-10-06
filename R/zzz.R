#.First.lib <- function(lib, pkg) library.dynam("sde", pkg, lib) 

#.noGenerics <- TRUE

.onLoad <- function(libname, pkgname)
{
 require(methods)
 
 require(zoo)
 cat("\n\nThis is the S4 version of YUIMA Project package.\nIt is intended for internal use only!\n\n")
# require(KernSmooth, quietly=TRUE)
# library.dynam("yuima", pkgname, libname) 
}

