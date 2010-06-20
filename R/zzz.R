#.First.lib <- function(lib, pkg) library.dynam("sde", pkg, lib) 

#.noGenerics <- TRUE

.onLoad <- function(libname, pkgname)
{
 require(methods)
 
 require(zoo)
 cat("\n\nThis is YUIMA Project package.\nNon stable release. Use at your own risk!!!\nFunctions and documentation may change in the final release\n\n")
# require(KernSmooth, quietly=TRUE)
# library.dynam("yuima", pkgname, libname) 
}

