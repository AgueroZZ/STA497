# LOAD LIBRARIES FOR CASE-CROSSOVER DATA ANALYSIS ###
# ALSO INCLUDE HELPER FUNCTIONS

# Packages required for computation
suppressWarnings(suppressMessages(require(ipoptr)))
library(parallel)
library(coxed)
library(Matrix)
library(trustOptim)
library(matrixStats) # for LogSumExp
library(purrr)
library(rstanarm)
library(INLA)

# tidyverse packages- required for analysis
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(lubridate)

# other
suppressWarnings(suppressMessages(require(rstanarm)))

### HELPER FUNCTIONS ###
