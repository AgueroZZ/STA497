# LOAD LIBRARIES FOR CASE-CROSSOVER DATA ANALYSIS ###
# ALSO INCLUDE HELPER FUNCTIONS
libplace = "~/lib"
# Packages required for computation
library(parallel)
library(Matrix)
library(trustOptim,lib.loc = libplace)
library(matrixStats) # for LogSumExp
library(purrr)
library(INLA,lib.loc = libplace)
library(survival)
library(mgcv)
library(mvQuad,lib.loc = libplace)
# tidyverse packages- required for analysis
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(lubridate)

