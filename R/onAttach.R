.onAttach=function(libname, pkgname){
  packageStartupMessage("Loaded DMRnet version ", as.character(utils::packageDescription("DMRnet")[["Version"]]), "\n")

  dmrnet_citation <- utils::citation("DMRnet")
  packageStartupMessage(dmrnet_citation$header)
  packageStartupMessage(paste(dmrnet_citation$textVersion))
}
