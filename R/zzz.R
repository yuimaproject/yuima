#.First.lib <- function(lib, pkg) library.dynam("sde", pkg, lib) 

#.noGenerics <- TRUE

.onLoad <- function(libname, pkgname)
{
 require(methods)
 
 require(zoo)
 cat(sprintf("\n\n%s\nThis is YUIMA Project package.\nNon stable release. Use at your own risk!!!\nFunctions and documentation may change in the final release\n%s\n\n",paste(rep("#",60),collapse=""),paste(rep("#",60),collapse="")))
# require(KernSmooth, quietly=TRUE)
# library.dynam("yuima", pkgname, libname) 
}

