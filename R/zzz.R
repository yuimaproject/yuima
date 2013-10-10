#.First.lib <- function(lib, pkg) library.dynam("sde", pkg, lib) 

#.noGenerics <- TRUE

.onAttach <- function(libname, pkgname)
{
    # require(methods)
 
    # require(zoo)
 packageStartupMessage(rep("#",60))
 packageStartupMessage("This is YUIMA Project package.")  
 packageStartupMessage("Non stable release. Use at your own risk!!!")
 packageStartupMessage("Functions and documentation may change in the final release")
 packageStartupMessage(rep("#",60))   
    
# require(KernSmooth, quietly=TRUE)
# library.dynam("yuima", pkgname, libname) 
}

