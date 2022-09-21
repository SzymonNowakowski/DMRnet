.onAttach=function(libname,pkgname){
  packageStartupMessage("Loaded DMRnet version ", as.character(packageDescription("DMRnet")[["Version"]]))
}
