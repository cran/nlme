### Copyright 1999-2003 Douglas M. Bates <bates@stat.wisc.edu>,
###                     Saikat DebRoy <saikat@stat.wisc.edu>

.First.lib <- function(lib, pkg) {
  library.dynam(pkg, pkg, lib )
  require(lattice)
  require(nls)
  autoload("dist", "mva")
}

.onLoad <- function(lib, pkg) {
  require(lattice)
  require(nls)
}
