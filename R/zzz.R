#.First.lib <- function(lib, pkg) library.dynam("sde", pkg, lib) 

#.noGenerics <- TRUE

.onAttach <- function(libname, pkgname)
{
    # require(methods)
 
    # require(zoo)
  Pver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
    fields="Version")
 packageStartupMessage(rep("#",40))
 packageStartupMessage(sprintf("This is YUIMA Project package v.%s", Pver))
 packageStartupMessage("Why don't you try yuimaGUI package?")
 packageStartupMessage("Visit: http://www.yuima-project.com")
 packageStartupMessage(rep("#",40))
    
# require(KernSmooth, quietly=TRUE)
# library.dynam("yuima", pkgname, libname) 
}

