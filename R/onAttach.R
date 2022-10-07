.onAttach=function(libname,pkgname){
  packageStartupMessage("Loaded DMRnet version ", as.character(utils::packageDescription("DMRnet")[["Version"]]))
}
