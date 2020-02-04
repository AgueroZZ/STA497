libplace = "~/lib"
install.packages("trustOptim",lib = libplace)
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE,lib = libplace)
install.packages("mvQuad",lib = libplace)
