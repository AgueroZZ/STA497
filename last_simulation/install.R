libplace = "~/lib"
install.packages("trustOptim",lib.loc = libplace)
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE,lib = libplace)
install.packages("mvQuad",lib.loc = libplace)